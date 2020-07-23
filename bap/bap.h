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
using namespace std;

class Bap
{
public:
    Bap(string input_bam, string output_path, string barcode_tag, int mapq, int cores, string run_name, bool tn5,
        double min_barcode_frags, double min_jaccard_index, string ref, string mito_chr, string bed_genome_file,
        string blacklist_file, string trans_file, bool species_mix);
    ~Bap() {};
    int taskflow();

private:
    string input_bam;
    string output_path;
    string barcode_tag;
    int mapq, cores;
    string run_name;
    bool tn5;
    double min_barcode_frags, min_jaccard_index;
    string ref;
    string mito_chr;
    string bed_genome_file, blacklist_file, trans_file;
    bool species_mix;
};