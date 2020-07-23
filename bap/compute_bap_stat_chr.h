/*
 * File: compute_bap_stat_chr.h
 * Created Data: 2020-7-13
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
using namespace std;

int compute_bap_stat_chr(int nc_thre, string HQbeads_file, string anno_bedpe_file,
    string out_frag_file, string out_nc_count_file, int regularize_thre);