/*
 * File: barcodeRank.h
 * Created Data: 2020-9-30
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <vector>

using namespace std;

enum CURVE_DATA_TYPE
{
    BEAD,
    JACCARD
};

enum INFLECTION_KERNEL_TYPE
{
    BAP2,
    DROPLETUTILS
};

double barcode_rank(vector< double >& input, INFLECTION_KERNEL_TYPE inflection_kernel_type, CURVE_DATA_TYPE curve_data_type);