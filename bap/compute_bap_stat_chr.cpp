/*
 * File: compute_bap_stat_chr.cpp
 * Created Data: 2020-7-13
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "compute_bap_stat_chr.h"

#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <iostream>

#include "gz_io.h"
#include "utility.h"

int compute_bap_stat_chr(int nc_thre, string HQbeads_file, string anno_bedpe_file,
    string out_frag_file, string out_nc_count_file, int regularize_thre)
{
    // Load HQ beads data
    ifstream hq_beads(HQbeads_file, std::ifstream::in);
    string line;
    set<string> beads;
    while (std::getline(hq_beads, line))
    {
        beads.insert(line);
    }
    hq_beads.close();

    // Load fragments
    gzFile anno_bedpe = gzopen(anno_bedpe_file.c_str(), "rb");
    map<unsigned long long, int> dict;
    vector<vector<string>> bedpe_data;
    while (readline(anno_bedpe, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() != 5) continue;
        // Filter for eligible barcodes
        if (beads.count(vec_s[4]) == 0) continue;

        int start = stoi(vec_s[2]);
        int end = stoi(vec_s[3]);
        ++dict[((unsigned long long)start << 32) + end];
        vec_s[0] = "";
        bedpe_data.push_back(vec_s);
    }
    gzclose(anno_bedpe);

    // Quantify NC + export
    map<int, int> cnts;
    for (auto& p : dict)
    {
        cnts[p.second] += p.second;
    }
    ofstream out_nc_count(out_nc_count_file, std::ofstream::out);
    for (auto& p : cnts)
    {
        out_nc_count<<p.first<<"\t"<<p.second<<endl;
    }
    out_nc_count.close();

    // Pull out barcode for retained fragments
    set<vector<string>> bedpe_data_filter;
    for (auto& v : bedpe_data)
    {
        int start = stoi(v[2]);
        int end = stoi(v[3]);
        if (dict[((unsigned long long)start << 32) + end] <= nc_thre)
            bedpe_data_filter.insert(v);
    }
    //cout<<"bedpe_data_filter size:"<<bedpe_data_filter.size()<<endl;
    bedpe_data.clear();
    bedpe_data.shrink_to_fit();

    map<int, vector<string>> overlap_start, overlap_end;
    for (auto& v : bedpe_data_filter)
    {
        int start = stoi(v[2]);
        int end = stoi(v[3]);
        overlap_start[start].push_back(v[4]);
        overlap_end[end].push_back(v[4]);
    }

    // devel
    size_t n = 0, no_dup = 0;
    map<string, int> bead_cnts;
    for (auto& overlap : {overlap_start, overlap_end})
    {
        for (auto& p : overlap)
        {
            auto& v = p.second;
            if (v.size() == 1) continue;
            n += v.size()*(v.size()-1);
            for (size_t i = 0; i < v.size(); ++i)
            {
                for (size_t j = i+1; j < v.size(); ++j)
                {
                    if (v[i] == v[j]) continue;
                    ++no_dup;
                    if (v[i] > v[j])
                        ++bead_cnts[v[i]+"\t"+v[j]];
                    else
                        ++bead_cnts[v[j]+"\t"+v[i]];
                }
            }
        }
    }
    //cout<<n<<" "<<no_dup<<endl;
    
    FILE * out_frag;
    out_frag = fopen(out_frag_file.c_str(), "w");
    int line_num = 0, total_num = 0;
    for (auto& p : bead_cnts)
    {
        int count = p.second;
        if (count >= regularize_thre)
        {
            ++line_num;
            total_num += count;
            string s = p.first + "\t" + to_string(count) + "\n";
            fwrite(s.c_str(), 1, s.size(), out_frag);
        }
    }
    //cout<<line_num<<" "<<total_num<<endl;
    fclose(out_frag);

    return 0;
}