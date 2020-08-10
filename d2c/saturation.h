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
    int addData(vector<Bedpe>& frags);

    // Calculate sequencing saturation
    int calculateSaturation(string out_file);

private:
    vector< FragData >                 _keys;  // Unique name, composed of barcode, gene and umi
    int                              _base2i[128];
    size_t _nreads;
    char _sep;
    vector< float > _samples;

    int _chr_num;
};