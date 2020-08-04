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

class Saturation
{
public:
    Saturation();
    ~Saturation();

    // Parse raw data to vectors
    int addData(int start, int end, string barcode, int chr_id);

    // Calculate sequencing saturation
    int calculateSaturation(string out_file);

private:
    vector< string >                 _keys;  // Unique name, composed of barcode, gene and umi
    int                              _base2i[128];
    size_t _nreads;
    char _sep;
    vector< float > _samples;
};