/*
 * File: d2c.h
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
#include <unordered_set>
#include <unordered_map>
using namespace std;

#include "samReader.h"
#include "saturation.h"
#include "types.h"

#include <filesystem>
namespace fs = std::filesystem;




struct SumStat
{
    int nuclear_total;
    int nuclear_uniq;
    int mito_total;
    int mito_uniq;
    int library_size;
    float dup_proportion;
    string drop_barcode;
};

class D2C
{
public:
    D2C(string input_bam, string output_path, string barcode_tag, int mapq, int cores, string run_name, bool tn5,
        double min_barcode_frags, double min_jaccard_index, string ref, string mito_chr, string bed_genome_file,
        string blacklist_file, string trans_file, bool species_mix, string bin_path, double barcode_threshold,
        double jaccard_threshold, bool saturation_on, string barcode_list);
    ~D2C() {};
    int run();
    int taskflow();
    int splitBamByChr(int chr_id);

private:
    map<string, string> parseChrsFromBedFile();
    void extractBedPE(const BamRecord b1, const BamRecord b2, vector<Bedpe>& bedpes, int l1, int l2);
    int determineHQBeads();
    pair<double, double> parseBeadThreshold(string filename);
    int computeStatByChr(int chr_id);
    int determineBarcodeMerge();
    int reannotateFragByChr(int chr_id);
    int annotateBamByChr(int chr_id);
    int finalQC();
    int plot();
    bool checkTn5(string s);
    bool checkTn5(size_t l);
    bool parseBarcodeList();
    inline string int2Barcode(int i);

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
    double barcode_threshold, jaccard_threshold;
    bool saturation_on;
    string barcode_list;

    // Specific parameters;
    int nc_threshold;   // Number of barcodes that a paired-end read must be observed for the read to be filtered
    int regularize_threshold; // Minimum number of inserts two barcodes must share to be counted per-chromosome
    bool one_to_one;    // Enforce that each bead barcode maps to one unique drop barcode (cancels the merging)
    fs::path temp_bam_path;
    string drop_tag;
    string peak_file;   // If supplied, compute FRIP (in QC stats) and generate Summarized Experiment

    // Shared data
    vector<vector<Bedpe>> _bedpes_by_chr;
    unordered_map<int, int> _total_bead_quant;
    //vector<string> _total_bead_order;
    std::mutex _merge_chr_mutex;
    unordered_set<int> _hq_beads;
    vector<map<int, int>> _total_nc_cnts;
    vector<map<size_t, int>> _total_bead_cnts;
    unordered_map<int, string> _drop_barcodes;
    vector<string> _contig_names;
    vector<vector<string>> _dup_frags;
    vector<string> _final_frags;
    vector<unordered_set<int>> _keep_qnames;
    vector<map<string, pair<int, int>>> _frag_stats;
    vector<SumStat> _sum_stats;

    // Instance of sequencing saturation
    Saturation saturation;

    // Barcode list data, for encoding and decodeing barcode strings
    vector<string> _barcode_names;
    unordered_map<string, int> _barcode2int;
};