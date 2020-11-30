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

struct AnnotateFragment
{
    AnnotateFragment(int _chr_id, int _start, int _end, int _barcode_id):
        chr_id(_chr_id), start(_start), end(_end), barcode_id(_barcode_id) {}
    bool operator==(const AnnotateFragment& o) const
    {
        return start == o.start && end == o.end && chr_id == o.chr_id && barcode_id == o.barcode_id;
    }

    int chr_id;
    int start;
    int end;
    int barcode_id;
};

namespace std
{
// Inject specialization of std::hash for UniqBarcode into namespace std
template <> struct hash< UniqBarcode >
{
    std::size_t operator()(UniqBarcode const& p) const
    {
        std::size_t seed = 0;
        spp::hash_combine(seed, p.start);
        spp::hash_combine(seed, p.end);
        spp::hash_combine(seed, p.barcode);
        return seed;
    }
};

template <> struct hash< AnnotateFragment >
{
    std::size_t operator()(AnnotateFragment const& frag) const
    {
        std::size_t seed = 0;
        spp::hash_combine(seed, frag.start);
        spp::hash_combine(seed, frag.end);
        spp::hash_combine(seed, frag.chr_id);
        spp::hash_combine(seed, frag.barcode_id);
        return seed;
    }
};
}  // namespace std