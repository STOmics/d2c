/*
 * File: bap.cpp
 * Created Data: 2020-7-23
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "bap.h"

Bap::Bap(string input_bam, string output_path, string barcode_tag, int mapq, int cores, string run_name, bool tn5,
        double min_barcode_frags, double min_jaccard_index, string ref, string mito_chr, string bed_genome_file,
        string blacklist_file, string trans_file, bool species_mix) :
        input_bam(input_bam),
        output_path(output_path),
        barcode_tag(barcode_tag),
        mapq(mapq),
        cores(cores),
        run_name(run_name),
        tn5(tn5),
        min_barcode_frags(min_barcode_frags),
        min_jaccard_index(min_jaccard_index),
        ref(ref),
        mito_chr(mito_chr),
        bed_genome_file(bed_genome_file),
        blacklist_file(blacklist_file),
        trans_file(trans_file),
        species_mix(species_mix) {}
    
int Bap::taskflow()
{
    return 0;
}

