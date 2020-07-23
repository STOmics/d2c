/*
 * File: assemble_fragments.cpp
 * Created Data: 2020-7-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "assemble_fragments.h"

#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

#include <htslib/hts.h>
#include <htslib/sam.h>

typedef bam1_t*    BamRecord;

const int FLAG = 2;
const int MAX_INSERT = 2000;
const int MAPQ = 20;

inline bool get_next_alignment(htsFile* fp, bam_hdr_t* header, BamRecord b)
{
    //cout<<"in get_next_alignment"<<endl;
    while (sam_read1(fp, header, b) > 0)
    {
        //cout<<b->core.flag<<endl;
        if ((b->core.flag & FLAG) == FLAG) 
            return true;
    }
    return false;
}

inline string get_qame(BamRecord b) 
{
    return bam_get_qname(b);
}

inline bool is_paired(BamRecord b)
{
    return b->core.flag & BAM_FPAIRED;
}

inline bool is_mapped(BamRecord b)
{
    return !(b->core.flag & BAM_FUNMAP);
}

inline bool is_reverse(BamRecord b)
{
    return b->core.flag & BAM_FREVERSE;
}

void print_bed_pe(const BamRecord b1, const BamRecord b2, const vector<string>& refs, ofstream& ofs)
{
    // Initialize BEDPE variables
    string chrom1, chrom2, strand1, strand2;
    int start1, start2, end1, end2;
    start1 = start2 = end1 = end2 = -1;
    chrom1 = chrom2 = strand1 = strand2 = ".";
    int min_map_quality = 0;

    // Extract relevant info for end 1
    if (is_mapped(b1))
    {
        chrom1 = refs[b1->core.tid];
        start1 = b1->core.pos;
        end1 = bam_endpos(b1);
        if (is_reverse(b1))
            strand1 = "-";
        else
            strand1 = "+";
    }

    // Extract relevant info for end 2
    if (is_mapped(b2))
    {
        chrom2 = refs[b2->core.tid];
        start2 = b2->core.pos;
        end2 = bam_endpos(b2);
        if (is_reverse(b2))
            strand2 = "-";
        else
            strand2 = "+";
    }

    // Swap the ends if necessary
    if (chrom1 > chrom2 || ((chrom1 == chrom2) && (start1 > start2)))
    {
        swap(chrom1, chrom2);
        swap(start1, start2);
        swap(end1, end2);
        swap(strand1, strand2);
    }

    // Report BEDPE using min mapping quality
    if (is_mapped(b1) && is_mapped(b2))
        min_map_quality = min(b1->core.qual, b2->core.qual);
    
    // Filter by max insert size and mapping quality
    if ((end2 - start1) < MAX_INSERT && min_map_quality >= MAPQ)
    {
        if (strand1 == "+")
        {
            start1 += 4;
            end2 += 4;
        }
        else
        {
            start1 -= 5;
            end2 -= 5;
        }
        ofs << chrom1 <<"\t"<<start1<<"\t"<<end2<<"\t"<<get_qame(b1)<<endl;
    }
}

// Return code:
// -1 failed open bam file
// -2 failed parse bam header
int assemble_fragments(const string bam_file, const string frag_file)
{
    // Open bam file
    htsFile* fp = hts_open(bam_file.c_str(), "r");
    if (!fp)
    {
        cout<<"Could not open "<<bam_file<<endl;
        return -1;
    }
    bam_hdr_t* header = sam_hdr_read(fp);
    if (header == nullptr)
    {
        cout<<"Could not parse header for "<<bam_file<<endl;
        return -2;
    }
    vector<string> refs;
    for (int i = 0; i < header->n_targets; ++i)
        refs.push_back(header->target_name[i]);

    // Iterate all reads
    // hts_itr_t* iter;
    // iter = sam_iter_queryi()
    ofstream ofs(frag_file, std::ofstream::out);

    BamRecord b1 = bam_init1();
    BamRecord b2 = bam_init1();
    while (get_next_alignment(fp, header, b1))
    {
        get_next_alignment(fp, header, b2);
        if (get_qame(b1) != get_qame(b2))
        {
            while (get_qame(b1) != get_qame(b2))
            {
                bam_copy1(b1, b2);
                get_next_alignment(fp, header, b2);
            }
            print_bed_pe(b1, b2, refs, ofs);
        }
        else if (is_paired(b1) && is_paired(b2))
        {
            print_bed_pe(b1, b2, refs, ofs);
        }
    }
    ofs.close();

    bam_destroy1(b1);
    bam_destroy1(b2);
    hts_close(fp);
    //hts_itr_destroy(iter);
    bam_hdr_destroy(header);

    return 0;
}