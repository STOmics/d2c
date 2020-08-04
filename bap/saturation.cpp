/*
 * File: saturation.cpp
 * Created Data: 2020-6-16
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "saturation.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <fstream>
#include <random>
#include <set>
#include <sstream>


Saturation::Saturation()
{
    _keys.clear();
    
    _sep         = '|';
    _nreads = 0;
    _samples = { 0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };

    _base2i['A'] = 0;
    _base2i['C'] = 1;
    _base2i['G'] = 2;
    _base2i['T'] = 3;
}

Saturation::~Saturation() {}

int Saturation::addData(int start, int end, string barcode, int chr_id)
{
    string temp = barcode+_sep+to_string(start)+_sep+to_string(end)+_sep+to_string(chr_id);
    _keys.push_back(temp);

    return 0;
}

int Saturation::calculateSaturation(string out_file)
{
    _nreads = _keys.size();
    spdlog::info("Calculate sequencing saturation");
    spdlog::debug("_nreads: {}",  _nreads);

    std::stringstream ss;
    ss << "#sample\tbar_x\tbar_y1\tbar_y2\n";

    std::random_device rd;
    std::mt19937       gen(rd());
    std::shuffle(_keys.begin(), _keys.end(), gen);

    unordered_map< string, unordered_map< string, int > > data;
    size_t p = 0;
    for (size_t i = 1; i < _samples.size(); ++i)
    {
        spdlog::info("Saturation sample:{}", _samples[i]);
        ss << _samples[i] << "\t";

        size_t size = size_t(_samples[i] * _nreads);
        for (; p < size; ++p)
        {
            string key = _keys[p];
            size_t pos = key.find(_sep);
            string barcode = key.substr(0, pos);
            string value = key.substr(pos+1);
            if (data.count(barcode) == 0)
                data[barcode] = {};
            ++data[barcode][value];
        }

        // Barcode
        size_t                   n_reads = 0;
        size_t                   n_uniq  = 0;
        vector< int >            n_frags;
        std::set< string > frags;
        for (auto& b : data)
        {
            frags.clear();
            for (auto& p : b.second)
            {
                n_reads += p.second;
                frags.insert(p.first);
                ++n_uniq;
            }
            n_frags.push_back(frags.size());
        }
        if (n_reads == 0)
        {
            spdlog::warn("Invalid data: n_reads == 0");
            continue;
        }
        
        std::nth_element(n_frags.begin(), n_frags.begin() + n_frags.size() / 2, n_frags.end());

        ss << n_reads / data.size() << "\t" << 1 - (n_uniq * 1.0 / n_reads) << "\t"
           << *(n_frags.begin() + n_frags.size() / 2) << std::endl;
    }
    
    std::ofstream ofs(out_file, std::ofstream::out);
    if (ofs.is_open())
    {
        ofs << ss.str();
        ofs.close();
        spdlog::info("Success dump saturation file:{}", out_file);
    }
    else
    {
        spdlog::error("Error opening file:{}", out_file);
    }

    _keys.clear();
    _keys.shrink_to_fit();

    return 0;
}