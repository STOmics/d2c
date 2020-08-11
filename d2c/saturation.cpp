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

    _base2i['A'] = 0;
    _base2i['C'] = 1;
    _base2i['G'] = 2;
    _base2i['T'] = 3;

    _chr_num = 0;
}

Saturation::~Saturation() {}

int Saturation::addData(vector< Bedpe >& frags)
{
    FragData fd;
    fd.chr = _chr_num;
    for (auto& frag : frags)
    {
        fd.start   = frag.start;
        fd.end     = frag.end;
        fd.barcode = frag.barcode;
        _keys.push_back(fd);
    }

    ++_chr_num;

    return 0;
}

int Saturation::calculateSaturation(string out_file)
{
    _nreads = _keys.size();
    spdlog::info("Calculate sequencing saturation");
    spdlog::debug("_nreads: {}", _nreads);

    std::stringstream ss;
    ss << "#sample\tbar_x\tbar_y1\tbar_y2\n";

    std::random_device rd;
    std::mt19937       gen(rd());
    std::shuffle(_keys.begin(), _keys.end(), gen);

    unordered_map< int, vector< unordered_map< size_t, int > > >           data;
    unordered_map< int, vector< unordered_map< size_t, int > > >::iterator it;
    size_t                                                                 p = 0;
    for (size_t i = 1; i < _samples.size(); ++i)
    {
        spdlog::info("Saturation sample:{}", _samples[i]);
        ss << _samples[i] << "\t";

        size_t size = size_t(_samples[i] * _nreads);
        for (; p < size; ++p)
        {
            auto&  frag    = _keys[p];
            int    barcode = frag.barcode;
            size_t value   = (( size_t )frag.start << 32) + frag.end;
            it             = data.find(barcode);
            if (it != data.end())
            {
                ++(it->second[frag.chr][value]);
            }
            else
            {
                data[barcode].resize(_chr_num, {});
                ++data[barcode][frag.chr][value];
            }
        }

        // Barcode
        size_t          n_reads = 0;
        size_t          n_uniq  = 0;
        map< int, int > uniq_cnt;
        for (auto& b : data)
        {
            auto& frags        = b.second;
            int   barcode_uniq = 0;
            for (int i = 0; i < _chr_num; ++i)
            {
                for (auto& p : frags[i])
                {
                    n_reads += p.second;
                    ++barcode_uniq;
                }
            }
            ++uniq_cnt[barcode_uniq];
        }

        if (n_reads == 0)
        {
            spdlog::warn("Invalid data: n_reads == 0");
            continue;
        }

        size_t median_num = 0;
        for (auto& p : uniq_cnt)
        {
            n_uniq += p.first * p.second;
            median_num += p.second;
        }

        median_num /= 2;
        size_t current_num = 0;
        int    median      = 0;
        for (auto& p : uniq_cnt)
        {
            current_num += p.second;
            if (current_num >= median_num)
            {
                median = p.first;
                break;
            }
        }

        ss << n_reads / data.size() << "\t" << 1 - (n_uniq * 1.0 / n_reads) << "\t" << median << std::endl;
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