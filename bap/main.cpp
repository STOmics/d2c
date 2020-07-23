/*
 * File: main.cpp
 * Created Data: 2020-7-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "assemble_fragments.h"
#include "annotate_fragments.h"
#include "compute_bap_stat_chr.h"
#include "determine_barcode_merges.h"
#include "reannotate_fragments_get_summary_stats_chr.h"
#include "reannotate_bam_file.h"
#include "final_qc_se.h"

#include <iostream>
#include <string>
using namespace std;

int main(int argc, char** argv)
{
    // if (argc != 3)
    // {
    //     cout<<"PROG: <BAM_FILE> <FRAG_FILE>"<<endl;
    //     return -1;
    // }

    // 01.assemble_fragments
    // string bam_file = argv[1];
    // string frag_file = argv[2];
    // assemble_fragments(bam_file, frag_file);

    // 02.annotate_fragments
    // string blacklist_file = argv[1];
    // string frag_file = argv[2];
    // string read_bead_file = argv[3];
    // string annotated_file = argv[4];
    // string stat_file = argv[5];
    // annotate_fragments(blacklist_file, frag_file, read_bead_file, annotated_file, stat_file);

    // 04.compute_bap_stat_chr
    // int nc_thre = atoi(argv[1]); 
    // string HQbeads_file = argv[2];
    // string anno_bedpe_file = argv[3];
    // string out_frag_file = argv[4];
    // string out_nc_count_file = argv[5];
    // int regularize_thre = atoi(argv[6]);
    // compute_bap_stat_chr(nc_thre, HQbeads_file, anno_bedpe_file, out_frag_file, out_nc_count_file, regularize_thre);

    if (argc > 10) cout<<"Too much paras"<<endl;

    // 05
    // string dir = argv[1];
    // string for_knee_file = argv[2];
    // string HQbead_file = argv[3];
    // string out_csv_file = argv[4];
    // string name = argv[5];
    // float minimum_jaccard_index = 0.000440958395334479f;
    // bool one_to_one = false;
    // bool barcoded_tn5 = false;

    // determine_barcode_merges(dir, for_knee_file,
    //     HQbead_file, out_csv_file, minimum_jaccard_index,
    //     name, one_to_one, barcoded_tn5);

    // 06
    // string anno_bedpe_file = argv[1];
    // string bc_translate_file = argv[2];
    // string clean_bedpe_file = argv[3];
    // string frag_ss_file = argv[4];
    // reannotate_fragments_get_summary_stats_chr(anno_bedpe_file, bc_translate_file,
    //     clean_bedpe_file, frag_ss_file);

    // 07
    // string input_bam_file = argv[1];
    // string output_bam_file = argv[2];
    // string bead_barcode_file = argv[3];
    // string HQbead_file = argv[4];
    // string bead_barcode = "CB", drop_barcode = "DB";
    // reannotate_bam_file(input_bam_file, output_bam_file,
    //     bead_barcode, drop_barcode,
    //     bead_barcode_file, HQbead_file);

    // 11
    string frag_file = argv[1];
    string tss_file = argv[2];
    string peak_file = "";
    string qc_basic_file = argv[3];
    string bc_translate_file = argv[4];
    bool species_mix = false;
    bool one_to_one = false;
    final_qc_se(frag_file, tss_file, peak_file, qc_basic_file, bc_translate_file, species_mix, one_to_one);


    return 0;
}