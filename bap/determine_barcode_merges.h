/*
 * File: determine_barcode_merges.h
 * Created Data: 2020-7-14
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
using namespace std;

int determine_barcode_merges(string dir, string for_knee_file,
    string HQbead_file, string out_csv_file, float minimum_jaccard_index,
    string name, bool one_to_one, bool barcoded_tn5);