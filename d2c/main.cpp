/*
 * File: main.cpp
 * Created Data: 2020-7-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "d2c.h"
#include <taskflow/taskflow.hpp>

#include <ctime>

#include <chrono>
#include <iostream>
#include <string>
#include <thread>
using namespace std;

#include <filesystem>
namespace fs = std::filesystem;

#include <CLI11.hpp>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

#include "timer.h"

//#define DEVEL
constexpr auto APP_NAME    = "D2C";
constexpr auto APP_VERSION = "1.0.0";

int main(int argc, char** argv)
{
    Timer timer;
    // Parse the command line parameters.
    CLI::App app{ string(APP_NAME) + ": Drop to Cell." };
    app.footer(string(APP_NAME) + " version: " + APP_VERSION);
    app.get_formatter()->column_width(40);

    // Required parameters
    string input_bam, output_path;
    app.add_option("-i", input_bam, "Input bam filename")->check(CLI::ExistingFile)->required();
    app.add_option("-o", output_path, "Output result path")->required();
    string barcode_list;
    app.add_option("-b", barcode_list, "Barcode list file")->check(CLI::ExistingFile)->required();

    // Optional parameters
    string barcode_tag = "XB";
    app.add_option("--bt", barcode_tag, "Barcode tag in bam file, default 'XB'");
    int mapq = 30;
    app.add_option("--mapq", mapq, "Filter thrshold of mapping quality, default 30");
    int cores = 0;
    app.add_option("-c", cores, "CPU core number, default detect");
    string run_name;
    app.add_option("-n", run_name, "Name for the all output files, default prefix of input bam file");

    bool tn5 = false;
    app.add_flag("--tn5", tn5, "Process data knowing that the barcodes were generated with a barcoded Tn5");
    double min_barcode_frags = 0.0;
    app.add_option("--bf", min_barcode_frags, "Minimum number of fragments to be thresholded for doublet merging");
    double min_jaccard_index = 0.0;
    app.add_option("--ji", min_jaccard_index, "Minimum jaccard index for collapsing bead barcodes to cell barcodes");

    // Model organism
    string ref = "hg19";
    app.add_option("-r", ref, "Specify supported reference genome, default hg19");

    // Non-model organism
    string mito_chr, bed_genome_file, blacklist_file, trans_file;
    app.add_option("--mc", mito_chr, "Name of the mitochondrial chromosome");
    app.add_option("--bg", bed_genome_file, "Bedtools genome file")->check(CLI::ExistingFile);
    app.add_option("--bl", blacklist_file, "Blacklist bed file")->check(CLI::ExistingFile);
    app.add_option("--ts", trans_file, "Path bed file of transcription start sites")->check(CLI::ExistingFile);

    // Specific parameters
    double barcode_threshold = 0.1;
    app.add_option("--bp", barcode_threshold,
                   "Percentage of minimum number of fragments to be thresholded for doublet merging")
        ->check(CLI::Range(0.0, 1.0));
    double jaccard_threshold = 0.3;
    app.add_option("--jp", jaccard_threshold,
                   "Percentage of minimum jaccard index for collapsing bead barcodes to cell barcodes")
        ->check(CLI::Range(0.0, 1.0));

    bool saturation_on = false;
    app.add_flag("--sat", saturation_on, "Output sequencing saturation file, default False");

    string barcode_runname_list = "";
    app.add_option("--br", barcode_runname_list, "Barcode runname list file, default detect")->check(CLI::ExistingFile);

    CLI11_PARSE(app, argc, argv);

    // Make sure the parameters are valid
    if (cores <= 0)
        cores = std::thread::hardware_concurrency();
    if (run_name.empty())
        run_name = fs::path(input_bam).stem().string();
    fs::path exe_path      = argv[0];
    exe_path               = fs::absolute(exe_path).parent_path();
    fs::path      ref_path = exe_path / "anno";
    set< string > supported_genomes;
    for (auto& p : fs::directory_iterator(ref_path / "bedtools"))
    {
        fs::path filename = p.path().filename();
        if (filename.extension() != ".sizes")
            continue;
        string name = filename.stem().string();
        if (name.substr(0, 6) == "chrom_")
            supported_genomes.insert(name.substr(6));
    }
    // Check jaccard index is valid
    if (min_jaccard_index > 1 || min_jaccard_index < 0)
    {
        cout << "User specified jaccard index > 1 or < 0:" << min_jaccard_index << endl;
    }
    // Handle reference genome
    if (supported_genomes.count(ref) != 0)
    {
        cout << "Found designated reference genome:" << ref << endl;
        if (trans_file.empty())
            trans_file = ref_path / "TSS" / (ref + ".refGene.TSS.bed");
        if (blacklist_file.empty())
            blacklist_file = ref_path / "blacklist" / (ref + ".full.blacklist.bed");
        if (bed_genome_file.empty())
            bed_genome_file = ref_path / "bedtools" / ("chrom_" + ref + ".sizes");
    }
    else
    {
        cout << "Could not identify this reference genome:" << ref << endl;
        cout << "Attempting to infer necessary input files from user specification" << endl;
        if (bed_genome_file.empty() || blacklist_file.empty() || trans_file.empty())
        {
            cout << "Invalid parameters:--bg --bl --ts" << endl;
            exit(1);
        }
    }

    // Make sure mito chr is valid
    if (mito_chr.empty())
    {
        if (ref == "hg19" || ref == "mm10" || ref == "hg38")
            mito_chr = "chrM";
        else if (ref == "GRCh37" || ref == "GRCh38" || ref == "GRCm37" || ref == "GRCm38")
            mito_chr = "MT";
        else if (ref == "hg19_mm10_c")
            mito_chr = "humanM";
        else
            mito_chr = "hg19_chrM";
    }
    // Check the output path is valid
    if (!fs::exists(output_path))
    {
        if (!fs::create_directories(fs::path(output_path)))
        {
            cout << "Failed to create directory: " << output_path << endl;
            exit(1);
        }
    }

    // Devel
#ifdef DEVEL
    for (auto& g : supported_genomes)
        cout << g << endl;
    cout << "cpu cores:" << cores << endl;
    cout << "run name:" << run_name << endl;
    cout << "bedtools:" << bed_genome_file << endl;
    cout << "blacklist:" << blacklist_file << endl;
    cout << "TSS:" << trans_file << endl;
    cout << "mito chr:" << mito_chr << endl;
#endif

    // Figure out if the specified reference genome is a species mix
    set< string > mix_species{ "hg19-mm10", "hg19_mm10_c", "hg19-mm10_nochr" };
    bool          species_mix = false;
    if (mix_species.count(ref) != 0)
        species_mix = true;

    // Set the default logger to file logger.
    std::time_t        t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream ostr;
    ostr << APP_NAME << "_" << std::put_time(std::localtime(&t), "%Y%m%d_%H%M%S") << ".log";
    try
    {
        // auto file_logger = spdlog::basic_logger_mt("main", "logs/" + ostr.str());
        auto file_sink      = std::make_shared< spdlog::sinks::basic_file_sink_mt >(exe_path / "logs" / ostr.str());
        auto main_logger    = std::make_shared< spdlog::logger >("main", file_sink);
        auto process_logger = std::make_shared< spdlog::logger >("process", file_sink);
        spdlog::register_logger(main_logger);
        spdlog::register_logger(process_logger);
        spdlog::set_default_logger(process_logger);
    }
    catch (const spdlog::spdlog_ex& ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }
    spdlog::set_level(spdlog::level::debug);  // Set global log level.
    spdlog::flush_on(spdlog::level::debug);
    spdlog::set_pattern("%Y-%m-%d %H:%M:%S.%e %L %n: %v");

    spdlog::get("main")->info("{} input_bam:{} output_path:{} barcode_tag:{} "
                              "mapq:{} cores:{} run_name:{} tn5:{} min_barcode_frags:{} min_jaccard_index:{} "
                              "ref:{} mito_chr:{} bed_genome_file:{} blacklist_file:{} trans_file:{} "
                              "species_mix:{} barcode_threshold:{} jaccard_threshold:{} saturation_on:{} "
                              "barcode_list:{} barcode_runname_list:{}",
                              argv[0], input_bam, output_path, barcode_tag, mapq, cores, run_name, tn5,
                              min_barcode_frags, min_jaccard_index, ref, mito_chr, bed_genome_file, blacklist_file,
                              trans_file, species_mix, barcode_threshold, jaccard_threshold, saturation_on,
                              barcode_list, barcode_runname_list);

    D2C d2c = D2C(input_bam, output_path, barcode_tag, mapq, cores, run_name, tn5, min_barcode_frags, min_jaccard_index,
                  ref, mito_chr, bed_genome_file, blacklist_file, trans_file, species_mix, exe_path.string(),
                  barcode_threshold, jaccard_threshold, saturation_on, barcode_list, barcode_runname_list);
    try
    {
        d2c.run();
    }
    catch (std::exception& e)
    {
        spdlog::error("Error: {}", e.what());
    }
    catch (...)
    {
        spdlog::error("Unknown error");
    }

    cout << "Finish process." << endl;
    spdlog::get("main")->info("{} process done. Elapsed time(s):{:.2f}", APP_NAME, timer.toc(1000));

    return 0;
}