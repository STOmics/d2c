/*
 * File: reannotate_fragments_get_summary_stats_chr.cpp
 * Created Data: 2020-7-15
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "reannotate_fragments_get_summary_stats_chr.h"

#include "gz_io.h"

#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

#include "utility.h"

inline int start(const string& s)
{
    size_t p1 = s.find("\t");
    size_t p2 = s.find("\t", p1+1);
    if (p2 == std::string::npos) return -1;
    return stoi(s.substr(p1+1, p2-p1));
}

struct cmp 
{
    bool operator() (const string& a, const string& b) const
    {
        return start(a) < start(b);
    }
};

int reannotate_fragments_get_summary_stats_chr(
    string anno_bedpe_file, string bc_translate_file,
    string clean_bedpe_file, string frag_ss_file
)
{
    // Load HQ beads data
    ifstream bc_translate(bc_translate_file, std::ifstream::in);
    string line;
    map<string, string> barcode_translate;
    while (std::getline(bc_translate, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() != 2) continue;
        barcode_translate[vec_s[0]] = vec_s[1];
    }
    bc_translate.close();

    // Load fragments
    gzFile anno_bedpe = gzopen(anno_bedpe_file.c_str(), "rb");
    map<string, string> pcr_dup;
    
    // Store n_total and n_unique as pair
    map<string, pair<int, int>> merge_ss;
    string chr;
    while (readline(anno_bedpe, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() != 5) continue;
        // Filter for eligible barcodes
        if (barcode_translate.count(vec_s[4]) == 0) continue;

        string cell_barcode = barcode_translate[vec_s[4]];
        string s = vec_s[1]+"\t"+vec_s[2]+"\t"+vec_s[3]+"\t"+cell_barcode;
        if (chr.empty()) chr = vec_s[1];
        //cout<<s<<" "<<start(s)<<endl;
        if (pcr_dup.count(s) == 0)
        {
            pcr_dup[s] = vec_s[0];
            ++merge_ss[cell_barcode].second;
        }

        ++merge_ss[cell_barcode].first;
    }
    gzclose(anno_bedpe);

    // Sort by start
    vector<string> dup_keys;
    for (auto& p : pcr_dup)
    {
        dup_keys.push_back(p.first);
    }
    std::sort(dup_keys.begin(), dup_keys.end(), cmp());

    // Export duplicated fragments
    FILE *frag_anno;
    frag_anno = fopen(clean_bedpe_file.c_str(), "w");
    for (auto& k : dup_keys)
    {
        string s = k + "\t" + pcr_dup[k] + "\n";
        fwrite(s.c_str(), 1, s.size(), frag_anno);
    }
    fclose(frag_anno);

    // Export summary statistics
    FILE *frag_ss;
    frag_ss = fopen(frag_ss_file.c_str(), "w");
    for (auto& p : merge_ss)
    {
        auto& pp = p.second;
        if (pp.first == 0 || pp.second == 0) continue;
        string s = p.first + "\t" + to_string(pp.first) + "\t" + to_string(pp.second)
            + "\t" + chr + "\n";
        fwrite(s.c_str(), 1, s.size(), frag_ss);
    }
    fclose(frag_ss);

    return 0;
}