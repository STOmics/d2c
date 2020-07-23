/*
 * File: reannotate_fragments_get_summary_stats_chr.h
 * Created Data: 2020-7-15
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
using namespace std;

int reannotate_fragments_get_summary_stats_chr(
    string anno_bedpe_file, string bc_translate_file,
    string clean_bedpe_file, string frag_ss_file
);