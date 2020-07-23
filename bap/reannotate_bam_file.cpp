/*
 * File: reannotate_bam_file.cpp
 * Created Data: 2020-7-15
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "reannotate_bam_file.h"

#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <iostream>

#include <htslib/hts.h>
#include <htslib/sam.h>

typedef bam1_t*    BamRecord;

#include "utility.h"

bool getTag(BamRecord b, const char tag[2], std::string& value)
{
    uint8_t* data = bam_aux_get(b, tag);
    if (data == NULL)
        return false;
    value = bam_aux2Z(data);
    return true;
}


int reannotate_bam_file(string input_bam_file, string output_bam_file,
    string bead_barcode, string drop_barcode,
    string dict_file, string hq_frags_file)
{
    // Load HQ beads data
    // Handle a dictionary of bead-barcode : drop-barcode pairs
    ifstream bc_translate(dict_file, std::ifstream::in);
    string line;
    map<string, string> barcode_translate;
    while (std::getline(bc_translate, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() != 2) continue;
        barcode_translate[vec_s[0]] = vec_s[1];
    }
    bc_translate.close();

    // Load fragment file
    // Import a set of read names that are in the HQ deduplicated set
    ifstream hq_frags(hq_frags_file, std::ifstream::in);
    set<string> keep_reads;
    while (std::getline(hq_frags, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() != 5) continue;
        keep_reads.insert(vec_s[4]);
    }
    hq_frags.close();

    // Iterate through bam

    // Open bam file
    htsFile* fp = hts_open(input_bam_file.c_str(), "r");
    if (!fp)
    {
        cout<<"Could not open "<<input_bam_file<<endl;
        return -1;
    }
    bam_hdr_t* header = sam_hdr_read(fp);
    if (header == nullptr)
    {
        cout<<"Could not parse header for "<<input_bam_file<<endl;
        return -2;
    }

    htsFile* out = hts_open(output_bam_file.c_str(), "wb");
    [[maybe_unused]] int ret = sam_hdr_write(out, header);

    BamRecord b = bam_init1();
    string bead_bc, drop_bc;
    while (sam_read1(fp, header, b) > 0)
    {
        if (!getTag(b, bead_barcode.c_str(), bead_bc)) continue;
        if (barcode_translate.count(bead_bc) == 0) continue;
        drop_bc = barcode_translate[bead_bc];

        // Handle droplet barcodes that we want to consider writing out
        string qname = bam_get_qname(b);
        if (keep_reads.count(qname) == 0) continue;
        bam_aux_append(b, drop_barcode.c_str(), 'Z', drop_bc.size()+1, ( uint8_t* )drop_bc.c_str());
        [[maybe_unused]] int ret = sam_write1(out, header, b);
    }
   
    bam_destroy1(b);
    hts_close(fp);
    bam_hdr_destroy(header);
    hts_close(out);

    return 0;
}