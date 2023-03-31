/*
 * File: samReader.cpp
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include <cmath>
#include <stdio.h>
#include <stdlib.h>

#include <exception>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>
using namespace std;

#include <filesystem>
namespace fs = std::filesystem;

#include <spdlog/spdlog.h>

#include "samReader.h"
#include "timer.h"
#include "utility.h"

// Class Implementations.
const int _DEFAULT_HTS_BLOCK_SIZE = 128 * (1024 * 1024);

// Class SamReader's functions.
SamReader::~SamReader()
{
    if (fp_)
    {
        Close();
    }
}

std::unique_ptr< SamReader > SamReader::FromFile(const std::string& reads_path)
{
    htsFile* fp = hts_open(reads_path.c_str(), "r");
    if (!fp)
    {
        char msg[128];
        sprintf(msg, "Could not open %s\n", reads_path.c_str());
        // printf(msg);
        std::runtime_error error(msg);
        throw std::exception(error);
    }

    // if (hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, _DEFAULT_HTS_BLOCK_SIZE) != 0)
    // {
    //     char msg[128];
    //     sprintf(msg, "Failed to set HTS_OPT_BLOCK_SIZE\n");
    //     // printf(msg);
    //     std::runtime_error error(msg);
    //     throw std::exception(error);
    // }

    bam_hdr_t* header = sam_hdr_read(fp);
    if (header == nullptr)
    {
        char msg[128];
        sprintf(msg, "Couldn't parse header for %s\n", fp->fn);
        // printf(msg);
        std::runtime_error error(msg);
    }

    // The htslib index data structure for our indexed BAM file. May be NULL if no
    // index was loaded.
    const int   threads_num = 4;
    std::string idx_path    = reads_path + ".bai";
    if (!fs::exists(idx_path) || check_file_older(idx_path, reads_path))
    {
        Timer                timer;
        [[maybe_unused]] int ret = sam_index_build3(reads_path.c_str(), idx_path.c_str(), 0, threads_num);
        spdlog::debug("Create bam index:{} rtn:{} time(s):{:.2f}", idx_path, ret, timer.toc(1000));
    }
    hts_idx_t* idx = sam_index_load(fp, fp->fn);
    return std::unique_ptr< SamReader >(new SamReader(fp, header, idx));
}

std::vector< std::pair< std::string, uint32_t > > SamReader::getContigs()
{
    return ref_;
}

bool SamReader::QueryByContig(int tid)
{
    // if (EXCLUDE_REFS.count(ref_[tid].first) != 0)
    //     return false;
    iter_ = sam_itr_queryi(idx_, tid, 0, ref_[tid].second);
    if (iter_ == nullptr || iter_->finished)
    {
        // spdlog::debug("No reads for ref:{}", ref_[tid].first);
        return false;
    }
    return true;
}

bool SamReader::next(BamRecord b)
{
    if (sam_itr_next(fp_, iter_, b) > 0)
    {
        return true;
    }
    return false;
}

bool SamReader::next(BamRecord b, hts_itr_t*& iter)
{
    int ret;
    if ((ret = sam_itr_next(fp_, iter, b)) > 0)
    {
        return true;
    }
    return false;
}

bool SamReader::next(BamRecord b, int flag)
{
    while (sam_itr_next(fp_, iter_, b) > 0)
    {
        // spdlog::debug("flag:{} qual:{}", b->core.flag, b->core.qual);
        if ((b->core.flag & flag) == flag)
            return true;
    }
    return false;
}

int SamReader::QueryOne(BamRecord b)
{
    for (size_t i = 0; i < ref_.size(); ++i)
    {
        hts_itr_t* iter = sam_itr_queryi(idx_, i, 0, ref_[i].second);
        if (iter == nullptr)
        {
            spdlog::warn("query unknown reference:{}", ref_[i].first);
            continue;
        }
        if (sam_itr_next(fp_, iter, b) > 0)
        {
            break;
        }
    }
    return 0;
}

SamReader::SamReader(htsFile* fp, bam_hdr_t* header, hts_idx_t* idx) : fp_(fp), header_(header), idx_(idx)
{
    iter_ = nullptr;
    for (int i = 0; i < header->n_targets; ++i)
        ref_.push_back({ header_->target_name[i], header_->target_len[i] });
}

int SamReader::Close()
{
    hts_close(fp_);
    fp_ = nullptr;

    if (idx_ != nullptr)
    {
        hts_idx_destroy(idx_);
        idx_ = nullptr;
    }
    if (iter_ != nullptr)
    {
        hts_itr_destroy(iter_);
        iter_ = nullptr;
    }

    bam_hdr_destroy(header_);
    header_ = nullptr;

    return 0;
}

BamHeader SamReader::getHeader()
{
    return header_;
}

void SamReader::setThreadPool(htsThreadPool* p)
{
    hts_set_opt(fp_, HTS_OPT_THREAD_POOL, p);
}