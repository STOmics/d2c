/*
 * File: final_qc_se.cpp
 * Created Data: 2020-7-16
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "final_qc_se.h"

#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "gz_io.h"
#include "utility.h"

#include "ygg.hpp"

#include <cmath>

using MyInterval = std::pair< int, int >;
template < class Node > class NodeTraits : public ygg::ITreeNodeTraits< Node >
{
public:
    using key_type = int;
    static int get_lower(const Node& node)
    {
        return node.lower;
    }
    static int get_upper(const Node& node)
    {
        return node.upper;
    }
    static int get_lower(const MyInterval& i)
    {
        return std::get< 0 >(i);
    }
    static int get_upper(const MyInterval& i)
    {
        return std::get< 1 >(i);
    }
};

class Node : public ygg::ITreeNodeBase< Node, NodeTraits< Node > >
{
public:
    int         upper;
    int         lower;
    string value;
};

using MyTree = ygg::IntervalTree< Node, NodeTraits< Node > >;

struct SummaryData
{
    int overlaps; // The number of has overlap with TSS
    float tss_proportion; // overlaps / total
    float mean_insert_size;
    int median_insert_size;
    float frip;   // The mean value of peak
    vector<int> insert_size;
};

float get_dup_proportion(int t, int u)
{
    return (t-u)*1.0 / t;
}

inline float aux(float x, int c, int n)
{
    return (c / x - 1 + exp(-n / x));
}
int get_library_size(int t, int u)
{
    float m = 1;
    float M = 100;
    int n_dup = t - u + 1;
    if (u > t || aux(m*u, u, t) < 0 || u < 0 || t < 0 || n_dup < 0)
        return 0;

    while (aux(M*u, u, t) > 0)
        M *= 10.0;
    
    for (int i = 0; i <= 40; ++i)
    {
        float mid = (m + M) / 2;
        float v = aux(mid * u, u, t);
        if ( v == 0 )
            break;
        else if (v > 0)
            m = mid;
        else
            M = mid;
    }
    return round(u*(m+M)/2.0);
}

int final_qc_se(string frag_file, string tss_file, string peak_file,
    string qc_basic_file, string bc_barcode_file,
    bool species_mix, bool one_to_one)
{
    // Load tss file for finding overlaps
    ifstream tss_ifs(tss_file, std::ifstream::in);
    string line;
    vector<Node> nodes;
    while (std::getline(tss_ifs, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() < 3) continue;

        Node node;
        node.lower = stoi(vec_s[1]) - 1000;
        node.upper = stoi(vec_s[2]) + 1000;
        node.value = vec_s[0];
        //cout<<node.value<<" "<<node.lower<<" "<<node.upper<<endl;
        nodes.push_back(std::move(node));
    }
    cout<<"step1"<<endl;

    map<string, MyTree> mytrees;
    for (auto& node : nodes)
    {
        //cout<<node.value<<" "<<node.lower<<" "<<node.upper<<endl;
        mytrees[node.value].insert(node);
    }

    // Load fragments
    gzFile frags = gzopen(frag_file.c_str(), "rb");
    
    // Store has_overlap and insert size as pair
    map<string, SummaryData> summary;
    
    while (readline(frags, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() < 4) continue;

        string chr = vec_s[0];
        int start = stoi(vec_s[1]);
        int end = stoi(vec_s[2]);
        MyInterval query_range {start, end};
        auto& sd = summary[vec_s[3]];
        if (mytrees.count(chr) != 0)
        {
            const auto& tree = mytrees.at(chr);
            const auto& res = tree.query(query_range);
            
            if (res.begin() != res.end())
                ++sd.overlaps;
        }
        int insert_size = end - start;
        sd.insert_size.push_back(insert_size);
    }
    gzclose(frags);
    cout<<"step2"<<endl;

    // Deal with FRIP if we have a peak file
    if (peak_file != "")
    {
        // TODO(fxzhao): implement this option
    }

    // Summarize frag attributes
    for (auto& p : summary)
    {
        auto& sd = p.second;
        sd.mean_insert_size = std::accumulate(sd.insert_size.begin(), sd.insert_size.end(), 0.0) / 
            sd.insert_size.size();
        std::nth_element(sd.insert_size.begin(), 
            sd.insert_size.begin()+int(0.5*sd.insert_size.size()),
            sd.insert_size.end());
        sd.median_insert_size = *(sd.insert_size.begin()+int(0.5*sd.insert_size.size()));
        sd.frip = 0;
        sd.tss_proportion = sd.overlaps * 1.0 / sd.insert_size.size();
    }

    // Import existing QC stats
    ifstream qc_ifs(qc_basic_file, std::ifstream::in);
    vector<vector<string>> qc;
    // Header line
    std::getline(qc_ifs, line);
    while (std::getline(qc_ifs, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() < 5) continue;
        int total_nuclear = stoi(vec_s[1]);
        int unique_nuclear = stoi(vec_s[2]);
        vec_s.push_back(to_string(get_dup_proportion(total_nuclear, unique_nuclear)));
        vec_s.push_back(to_string(get_library_size(total_nuclear, unique_nuclear)));
        qc.push_back(vec_s);
    }
    qc_ifs.close();

    // Add barcodes back if we need to (specified with the one to one option)
    if (one_to_one)
    {
        // TODO(fxzhao)
    }

    // Create a summarized experiment if we have a valid peak file
    if (peak_file != "")
    {
        // TODO(fxzhao)
    }

    // Export QC stats
    FILE * qc_out;
    string replace = ".basicQC.tsv";
    qc_basic_file.replace(qc_basic_file.find(replace), replace.size(), ".QCstats.csv");
    qc_out = fopen(qc_basic_file.c_str(), "w");
    for (auto& q : qc)
    {
        if (summary.count(q[0]) == 0) continue;

        string s = q[0]+","+q[1]+","+q[2]+","+q[3]+","+q[4]+","+q[5]+","+q[6]+",";
        auto& sd = summary[q[0]];
        s += to_string(sd.mean_insert_size)+","+to_string(sd.median_insert_size)+","+
            to_string(sd.tss_proportion)+","+to_string(sd.frip)+"\n";
        fwrite(s.c_str(), 1, s.size(), qc_out);
    }
    fclose(qc_out);

    return 0;
}