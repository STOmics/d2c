/*
 * File: final_qc_se.h
 * Created Data: 2020-7-16
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
using namespace std;

int final_qc_se(string frag_file, string tss_file, string peak_file,
    string qc_basic_file, string bc_barcode_file,
    bool species_mix, bool one_to_one);