/*
 * File: bap.h
 * Created Data: 2020-7-23
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
#include <map>
#include <mutex>
#include <set>
using namespace std;

#include "samReader.h"

#include <filesystem>
namespace fs = std::filesystem;


struct Bedpe
{
    int start;
    int end;
    string qname;
    string barcode;
};


class Bap
{
public:
    Bap(string input_bam, string output_path, string barcode_tag, int mapq, int cores, string run_name, bool tn5,
        double min_barcode_frags, double min_jaccard_index, string ref, string mito_chr, string bed_genome_file,
        string blacklist_file, string trans_file, bool species_mix, string bin_path);
    ~Bap() {};
    int run();
    int taskflow();
    int splitBamByChr(int chr_id);

private:
    map<string, string> parseChrsFromBedFile();
    void extractBedPE(const BamRecord b1, const BamRecord b2, vector<Bedpe>& bedpes);
    int determineHQBeads();
    double parseBeadThreshold(string filename);
    int computeStatByChr(int chr_id);
    int determineBarcodeMerge();
    int reannotateFragByChr(int chr_id);
    int annotateBamByChr(int chr_id);

private:
    // Input parameters
    string input_bam;
    fs::path output_path;
    string barcode_tag;
    int mapq, cores;
    string run_name;
    bool tn5;
    double min_barcode_frags, min_jaccard_index;
    string ref;
    string mito_chr;
    string bed_genome_file, blacklist_file, trans_file;
    bool species_mix;
    fs::path bin_path;

    // Specific parameters;
    int nc_threshold;   // Number of barcodes that a paired-end read must be observed for the read to be filtered
    int regularize_threshold; // Minimum number of inserts two barcodes must share to be counted per-chromosome
    bool one_to_one;    // Enforce that each bead barcode maps to one unique drop barcode (cancels the merging)
    fs::path temp_bam_path;
    string drop_tag;

    // Shared data
    vector<vector<Bedpe>> _bedpes_by_chr;
    map<string, int> _total_bead_quant;
    std::mutex _merge_chr_mutex;
    set<string> _hq_beads;
    vector<map<int, int>> _total_nc_cnts;
    vector<map<string, int>> _total_bead_cnts;
    map<string, string> _drop_barcodes;
    vector<string> _contig_names;
    vector<vector<string>> _dup_frags;
    vector<set<string>> _keep_qnames;
};