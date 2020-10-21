/*
 * File: saturation.h
 * Created Data: 2020-6-16
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <mutex>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
using namespace std;

#include "types.h"

struct FragData
{
    int barcode;
    int start;
    int end;
    int chr;
};

class Saturation
{
public:
    Saturation();
    ~Saturation();

    // Parse raw data to vectors
    int addData(unordered_map< string, unordered_map< string, int > >& dups_per_cell);

    // Calculate sequencing saturation
    int calculateSaturation(string out_file);

private:
    vector< size_t > _keys;  // Unique name, composed of cell barcode and uniq id
    size_t           _nreads;
    char             _sep;
    vector< float >  _samples;

    int                          _barcode_index;
    int                          _frags_index;
    unordered_map< string, int > _uniq_barcodes;
};