/*
 * File: saturation.cpp
 * Created Data: 2020-6-16
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "saturation.h"
#include "utility.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <sstream>

Saturation::Saturation()
{
    _keys.clear();

    _sep     = '|';
    _nreads  = 0;
    _samples = { 0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };

    _barcode_index = 0;
    _frags_index   = 0;
}

Saturation::~Saturation() {}

int Saturation::addData(unordered_map< string, spp::sparse_hash_map< AnnotateFragment, int > >& dups_per_cell)
{
    for (auto& [barcode, frags] : dups_per_cell)
    {
        if (_uniq_barcodes.count(barcode) == 0)
        {
            _uniq_barcodes.insert({ barcode, _barcode_index });
            ++_barcode_index;
        }
        size_t key = (size_t(_uniq_barcodes[barcode]) << 32);
        for (auto& [_, count] : frags)
        {
            size_t new_key = (key | _frags_index);
            ++_frags_index;
            for (int i = 0; i < count; ++i)
                _keys.push_back(new_key);
            // spdlog::info("{} {} {}", barcode, key, new_key);
        }
    }
    return 0;
}

int Saturation::calculateSaturation(string out_file)
{
    _nreads = _keys.size();
    spdlog::info("Calculate sequencing saturation");
    spdlog::info("Saturation total reads: {}", _nreads);

    std::stringstream ss;
    ss << "#sample_ratio\tmean_frags_per_cell\tsaturation\tmedian_uniq_frags_per_cell";

    std::random_device rd;
    std::mt19937       gen(rd());
    std::shuffle(_keys.begin(), _keys.end(), gen);

    unordered_map< int, unordered_map< int, int > >           data;
    unordered_map< int, unordered_map< int, int > >::iterator it;
    size_t                                                    p = 0;
    for (size_t i = 1; i < _samples.size(); ++i)
    {
        spdlog::info("Saturation sample:{}", _samples[i]);
        ss << "\n" << _samples[i] << "\t";

        size_t size = size_t(_samples[i] * _nreads);
        for (; p < size; ++p)
        {
            auto& frag     = _keys[p];
            int   barcode  = (frag >> 32);
            int   uniq_key = (frag & 0xFFFFFFFF);
            // it             = data.find(barcode);
            // if (it != data.end())
            // {
            //     ++(it->second[uniq_key]);
            // }
            // else
            // {
            //     ++data[barcode][uniq_key];
            // }
            ++data[barcode][uniq_key];
        }

        // Barcode
        size_t          n_reads = 0;
        size_t          n_uniq  = 0;
        map< int, int > uniq_cnt;  // Automatic sorting
        for (auto& [_, frags] : data)
        {
            int barcode_uniq = 0;
            for (auto& [_, count] : frags)
            {
                n_reads += count;
                ++barcode_uniq;
            }

            ++uniq_cnt[barcode_uniq];
        }

        if (n_reads == 0)
        {
            spdlog::warn("Invalid data: n_reads == 0");
            continue;
        }

        size_t median_num = 0;
        for (auto& [uniq_barcode_num, count] : uniq_cnt)
        {
            n_uniq += uniq_barcode_num * count;
            median_num += count;
        }

        median_num /= 2;
        size_t current_num = 0;
        int    median      = 0;
        for (auto& [uniq_barcode_num, count] : uniq_cnt)
        {
            current_num += count;
            if (current_num >= median_num)
            {
                median = uniq_barcode_num;
                break;
            }
        }
        // spdlog::info("{} {} {}", data.size(), n_reads, n_uniq);
        ss << n_reads / data.size() << "\t" << 1 - (n_uniq * 1.0 / n_reads) << "\t" << median;
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