/*
 * File: types.h
 * Created Data: 2020-8-7
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <sparsepp/spp_utils.h>

struct Bedpe
{
    int start;  // start and end positions for this fragment
    int end;

    int qname1;  // record reads position for keep reads in bam
    int qname2;

    int barcode;  // encode runname,bc1,bc2 to int32: 8+12+12
};

struct UniqBarcode
{
    bool operator==(const UniqBarcode& o) const
    {
        return start == o.start && end == o.end && barcode == o.barcode;
    }
    int start;
    int end;
    int barcode;
};

namespace std
{
    // Inject specialization of std::hash for UniqBarcode into namespace std
    template<>
    struct hash<UniqBarcode>
    {
        std::size_t operator()(UniqBarcode const &p) const
        {
            std::size_t seed = 0;
            spp::hash_combine(seed, p.start);
            spp::hash_combine(seed, p.end);
            spp::hash_combine(seed, p.barcode);
            return seed;
        }
    };
}