/*
 * File: d2c.h
 * Created Data: 2020-7-23
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
using namespace std;

#include "samReader.h"
#include "saturation.h"
#include "types.h"

#include <filesystem>
namespace fs = std::filesystem;

struct SumStat
{
    int    nuclear_total;
    int    nuclear_uniq;
    int    mito_total;
    int    mito_uniq;
    int    human_mito_total;
    int    human_mito_uniq;
    int    mouse_mito_total;
    int    mouse_mito_uniq;
    int    human_total;
    int    human_uniq;
    int    mouse_total;
    int    mouse_uniq;
    int    library_size;
    float  dup_proportion;
    string drop_barcode;
};

class D2C
{
public:
    D2C(string input_bam, string output_path, string barcode_in_tag, string barcode_out_tag, int mapq, int cores,
        string run_name, int tn5, double min_barcode_frags, double min_jaccard_index, string ref, string mito_chr,
        string bed_genome_file, string blacklist_file, string trans_file, bool species_mix, string bin_path,
        int barcode_threshold, int jaccard_threshold, bool saturation_on, string barcode_list,
        string barcode_runname_list, int beads_force, string tn5_list);
    ~D2C(){};
    int run();
    int taskflow();
    int splitBamByChr(int chr_id);

private:
    map< string, string > parseChrsFromBedFile(string filename);
    void                  extractBedPE(const BamRecord b1, const BamRecord b2, vector< Bedpe >& bedpes, int l1, int l2);
    int                   determineHQBeads();
    pair< double, double > parseBeadThreshold(string filename);
    int                    computeStatByChr(int chr_id);
    int                    determineBarcodeMerge();
    int                    reannotateFragByChr(int chr_id);
    int                    annotateBamByChr(int chr_id);
    int                    simpleQC(vector<int>& used_chrs);
    int                    finalQC();
    int                    plot();
    bool                   checkTn5(string s);
    bool                   checkTn5(size_t l);
    bool                   parseBarcodeList();
    inline string          int2Barcode(int i);
    inline string          int2Tn5(int i);
    bool                   parseRunnameList();
    bool    parseTn5List();

private:
    // Input parameters
    string   input_bam;
    fs::path output_path;
    string   barcode_tag;
    string   drop_tag;
    int      mapq, cores;
    string   run_name;
    int     tn5;
    double   min_barcode_frags, min_jaccard_index;
    string   ref;
    string   bed_genome_file, blacklist_file, trans_file;
    bool     species_mix;
    fs::path bin_path;
    int      barcode_threshold, jaccard_threshold;
    bool     saturation_on;
    string   barcode_list;
    string   barcode_runname_list;
    int  beads_force;
    string tn5_list;

    // for mixed species
    string   single_mc, human_mc, mouse_mc;
    set<string> mito_chrs;

    // Specific parameters;
    int      nc_threshold;  // Number of barcodes that a paired-end read must be observed for the read to be filtered
    int      regularize_threshold;  // Minimum number of inserts two barcodes must share to be counted per-chromosome
    bool     one_to_one;  // Enforce that each bead barcode maps to one unique drop barcode (cancels the merging)
    fs::path temp_bam_path;
    string   peak_file;  // If supplied, compute FRIP (in QC stats) and generate Summarized Experiment

    // Shared data
    vector< vector< Bedpe > > _bedpes_by_chr;
    unordered_map< int, int > _total_bead_quant;
    // vector<string> _total_bead_order;
    std::mutex                                _merge_chr_mutex;
    unordered_set< int >                      _hq_beads;
    vector< map< int, int > >                 _total_nc_cnts;
    vector< unordered_map< size_t, int > >    _total_bead_cnts;
    unordered_map< int, int >              _drop_barcodes;
    vector< string >                          _contig_names;
    vector< vector< AnnotateFragment > >                _dup_frags;
    vector< AnnotateFragment >                          _final_frags;
    vector< unordered_set< int > >            _keep_qnames;
    vector< map< string, pair< unsigned long, unsigned long > > > _frag_stats; // ulong encode two uint32, for mix species
    vector< SumStat >                         _sum_stats;

    // Instance of sequencing saturation
    Saturation saturation;

    // Barcode list data, for encoding and decodeing barcode strings
    vector< string >             _barcode_names;
    unordered_map< string, int > _barcode2int;

    // Barcode runnames, just keep run name in bead barcodes
    // Note: runnames starts with '-', and the string in first position is empty string
    // vector< string >             _runnames;
    // unordered_map< string, int > _runname2int;
    string _runname;

    // Index to drop barcodes, encoding barcodes to int
    vector<string> _idx2drop;

    // Check input file type
    bool _is_bed;

    // Encode tn5 sequence to integer
    vector<string> _tail_names;
    unordered_map<string, int> _tail2int;
};
