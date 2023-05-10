/*
 * File: main.cpp
 * Created Data: 2020-7-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "d2c.h"
#include "reannotate.h"
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
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include "timer.h"

//#define DEVEL
constexpr auto APP_NAME    = "D2C";
constexpr auto APP_VERSION = "1.5.1";

int main(int argc, char** argv)
{
    Timer    timer;
    fs::path exe_path = argv[0];
    exe_path          = fs::absolute(exe_path).parent_path();

    // Parse the command line parameters.
    CLI::App app{ string(APP_NAME) + ": Drop to Cell." };
    app.footer(string(APP_NAME) + " version: " + APP_VERSION);
    app.get_formatter()->column_width(40);

    // Require and only require single subcommand
    app.require_subcommand(1);

    auto sub_count  = app.add_subcommand("merge", "Drop barcode to Cell barcode");
    auto sub_reanno = app.add_subcommand("transid", "Reannotate bam file using barcode translate file");
    sub_count->fallthrough();
    sub_reanno->fallthrough();

    // Required parameters
    string input_bam, output_path;
    app.add_option("-i", input_bam, "Input bam or bed filename")->check(CLI::ExistingFile)->required();
    app.add_option("-o", output_path, "Output result path")->required();

    // Optional parameters
    string barcode_in_tag = "XB";
    app.add_option("--bt1", barcode_in_tag, "Barcode tag in input bam file, default 'XB'");
    string barcode_out_tag = "DB";
    app.add_option("--bt2", barcode_out_tag, "Barcode tag in output bam file, default 'DB'");
    string log_path = "logs";
    app.add_option("--log", log_path, "Set logging path, default is './logs'");
    string run_name;
    app.add_option("-n", run_name, "Name for the all output files, default prefix of input bam file");

    string barcode_list = (exe_path / "barcode.list").string();
    sub_count->add_option("-b", barcode_list, "Barcode list file")->check(CLI::ExistingFile);
    string tn5_list = (exe_path / "tn5.list").string();
    int    mapq     = 30;
    sub_count->add_option("--mapq", mapq, "Filter thrshold of mapping quality, default 30");
    int cores = std::thread::hardware_concurrency();
    sub_count->add_option("-c", cores, "CPU core number, default detect")->check(CLI::PositiveNumber);
    int tn5 = 0;
    sub_count
        ->add_option("--tn5", tn5,
                     "Process data knowing that the barcodes were generated with a barcoded Tn5, default 0")
        ->check(CLI::NonNegativeNumber);
    double min_barcode_frags = 0.0;
    sub_count->add_option("--bf", min_barcode_frags,
                          "Minimum number of fragments to be thresholded for doublet merging");
    double min_jaccard_index = 0.0;
    sub_count->add_option("--ji", min_jaccard_index,
                          "Minimum jaccard index for collapsing bead barcodes to cell barcodes");

    // Model organism
    string ref = "";
    sub_count->add_option("-r", ref, "Specify supported reference genome, default None");

    // Non-model organism
    string mito_chr, bed_genome_file, blacklist_file, trans_file;
    sub_count->add_option("--mc", mito_chr, "Name of the mitochondrial chromosome, seperated by ',' for mixed species");
    sub_count->add_option("--bg", bed_genome_file, "Bedtools genome file")->check(CLI::ExistingFile);
    sub_count->add_option("--bl", blacklist_file, "Blacklist bed file")->check(CLI::ExistingFile);
    sub_count->add_option("--ts", trans_file, "Path bed file of transcription start sites")->check(CLI::ExistingFile);

    // Specific parameters
    int barcode_threshold = 0;
    sub_count->add_option("--bp", barcode_threshold, "Top N number of fragments to be thresholded for doublet merging")
        ->check(CLI::PositiveNumber);

    int beads_force = 0;
    sub_count
        ->add_option("--fb", beads_force,
                     "Top N number of fragments to be thresholded when calculated threshold is too large")
        ->check(CLI::NonNegativeNumber);

    int jaccard_threshold = 0;
    sub_count
        ->add_option("--jp", jaccard_threshold,
                     "Top N number of jaccard index for collapsing bead barcodes to cell barcodes")
        ->check(CLI::PositiveNumber);

    bool saturation_on = false;
    sub_count->add_flag("--sat", saturation_on, "Output sequencing saturation file, default False");

    bool species_mix = false;
    sub_count->add_flag("--mix-species", species_mix, "Set species mixed, default False");

    string rank = "inflection";
    sub_count->add_option("--rank", rank, "Rank method, support knee/inflection/kde, default knee")
        ->check(CLI::IsMember({ "knee", "inflection" }));

    string barcode_runname_list = "";
    // app.add_option("--br", barcode_runname_list, "Barcode runname list file, default
    // detect")->check(CLI::ExistingFile);

    // Reannotate specific parameters
    string barcode_translate_file;
    sub_reanno->add_option("-t", barcode_translate_file, "Translate file from drop barcode to cell barcode")
        ->check(CLI::ExistingFile)
        ->required();

    CLI11_PARSE(app, argc, argv);

    // Set the default logger to file logger.
    std::time_t        t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::ostringstream ostr;
    ostr << APP_NAME << "_" << std::put_time(std::localtime(&t), "%Y%m%d_%H%M%S") << ".log";
    try
    {
        // auto file_logger = spdlog::basic_logger_mt("main", "logs/" + ostr.str());
        fs::path fs_log_path(log_path);
        auto     file_sink    = std::make_shared< spdlog::sinks::basic_file_sink_mt >(fs_log_path / ostr.str());
        auto     console_sink = std::make_shared< spdlog::sinks::stdout_color_sink_mt >();
        std::vector< spdlog::sink_ptr > sinks{ console_sink, file_sink };
        auto                            logger = std::make_shared< spdlog::logger >("d2c", sinks.begin(), sinks.end());
        spdlog::register_logger(logger);
        spdlog::set_default_logger(logger);
    }
    catch (const spdlog::spdlog_ex& ex)
    {
        std::cout << "Log init failed: " << ex.what() << std::endl;
    }
    spdlog::set_level(spdlog::level::debug);  // Set global log level.
    spdlog::flush_on(spdlog::level::debug);
    spdlog::set_pattern("%Y-%m-%d %H:%M:%S.%e %L %n: %v");

    // Get runname from the prefix of input bam if it is empty
    if (run_name.empty())
        run_name = fs::path(input_bam).stem().string();
    // Check the output path is valid
    if (!fs::exists(output_path))
    {
        if (!fs::create_directories(fs::path(output_path)))
        {
            cout << "Failed to create directory: " << output_path << endl;
            exit(1);
        }
    }

    if (sub_count->parsed())
    {
        // Make sure the parameters are valid
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
            spdlog::info("User specified jaccard index > 1 or < 0: {}", min_jaccard_index);
        }
        // Handle reference genome
        if (supported_genomes.count(ref) != 0)
        {
            spdlog::info("Found designated reference genome: {}", ref);
            if (trans_file.empty())
                trans_file = ref_path / "TSS" / (ref + ".refGene.TSS.bed");
            if (blacklist_file.empty())
                blacklist_file = ref_path / "blacklist" / (ref + ".full.blacklist.bed");
            if (bed_genome_file.empty())
                bed_genome_file = ref_path / "bedtools" / ("chrom_" + ref + ".sizes");
        }
        else
        {
            // Allow no clear reference parameter
            // cout << "Could not identify this reference genome:" << ref << endl;
            // cout << "Attempting to infer necessary input files from user specification" << endl;
            // if (bed_genome_file.empty() || blacklist_file.empty() || trans_file.empty())
            // {
            //     cout << "Invalid parameters:--bg --bl --ts" << endl;
            //     exit(1);
            // }
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
            // Remove the default mito chr
            // else
            //     mito_chr = "hg19_chrM";
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
        if (mix_species.count(ref) != 0)
            species_mix = true;

        spdlog::info("{} input_bam:{} output_path:{} barcode_in_tag:{} barcode_out_tag:{} "
                     "mapq:{} cores:{} run_name:{} tn5:{} min_barcode_frags:{} min_jaccard_index:{} "
                     "ref:{} mito_chr:{} bed_genome_file:{} blacklist_file:{} trans_file:{} "
                     "species_mix:{} barcode_threshold:{} jaccard_threshold:{} saturation_on:{} "
                     "barcode_list:{} barcode_runname_list:{} beads_force:{} tn5_list:{} rank:{}",
                     argv[0], input_bam, output_path, barcode_in_tag, barcode_out_tag, mapq, cores, run_name, tn5,
                     min_barcode_frags, min_jaccard_index, ref, mito_chr, bed_genome_file, blacklist_file, trans_file,
                     species_mix, barcode_threshold, jaccard_threshold, saturation_on, barcode_list,
                     barcode_runname_list, beads_force, tn5_list, rank);

        D2C d2c = D2C(input_bam, output_path, barcode_in_tag, barcode_out_tag, mapq, cores, run_name, tn5,
                      min_barcode_frags, min_jaccard_index, ref, mito_chr, bed_genome_file, blacklist_file, trans_file,
                      species_mix, exe_path.string(), barcode_threshold, jaccard_threshold, saturation_on, barcode_list,
                      barcode_runname_list, beads_force, tn5_list, rank);
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
    }
    else if (sub_reanno->parsed())
    {
        spdlog::info(
            "{} input_bam:{} output_path:{} barcode_in_tag:{} barcode_out_tag:{} runname:{} barcode_translate_file:{}",
            argv[0], input_bam, output_path, barcode_in_tag, barcode_out_tag, run_name, barcode_translate_file);
        fs::path output_bam(output_path);
        // Add suffix so the input bam and output bam can be in the same directory
        output_bam /= run_name + "_transid.bam";
        bool ret = reannotate(input_bam, barcode_translate_file, output_bam.string(), barcode_in_tag, barcode_out_tag);
        if (!ret)
            spdlog::error("Transform failed!");
    }
    else
    {
        cerr << "Undefined subcommand." << endl;
        return -1;
    }

    // cout << "Finish process." << endl;
    spdlog::info("{} process done. Elapsed time(s):{:.2f}", APP_NAME, timer.toc(1000));

    return 0;
}
