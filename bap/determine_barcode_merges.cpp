/*
 * File: determine_barcode_merges.cpp
 * Created Data: 2020-7-14
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "determine_barcode_merges.h"

#include <filesystem>
namespace fs = std::filesystem;

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "utility.h"

#include <cmath>

int determine_barcode_merges(string dir, string for_knee_file,
    string HQbead_file, string out_csv_file, float minimum_jaccard_index,
    string name, bool one_to_one, bool barcoded_tn5)
{
    // Load all rds files
    map<string, int> sum_dt;
    for (auto& p : fs::directory_iterator(dir))
    {
        if (p.path().extension() == ".rds")
        {
            cout<<p.path()<<endl;
            ifstream ifs(p.path(), std::ifstream::in);
            string line;
            while (std::getline(ifs, line))
            {
                vector<string> vec_s = split_str(line, '\t');
                if (vec_s.size() != 3) continue;
                sum_dt[vec_s[0]+','+vec_s[1]] += stoi(vec_s[2]);
            }
            ifs.close();
        }
    }

    // Only consider merging when Tn5 is the same
    if (barcoded_tn5)
    {
        // TODO(fxzhao): implement tn5 option
    }

    // Import number of barcodes
    ifstream HQbead(HQbead_file, std::ifstream::in);
    string line;
    set<string> valid_barcodes;
    while (std::getline(HQbead, line))
    {
        valid_barcodes.insert(line);
    }
    HQbead.close();
    ifstream for_knee(for_knee_file, std::ifstream::in);
    vector<pair<string, int>> nBC;
    map<string, int> count_dict;
    while (std::getline(for_knee, line))
    {
        vector<string> vec_s = split_str(line, ',');
        if (vec_s.size() != 2) continue;
        if (valid_barcodes.count(vec_s[0]) == 0) continue;
        nBC.push_back({vec_s[0], stoi(vec_s[1])});
        count_dict[vec_s[0]] = stoi(vec_s[1]) * 2;
    }
    for_knee.close();
    // Sort by count
    std::sort(nBC.begin(), nBC.end(), [](pair<string, int>& a, pair<string, int>& b) {
        return a.second > b.second;
    });
    cout<<"nBC size:"<<nBC.size()<<" valid_barcodes size:"<<valid_barcodes.size()<<endl;
    vector<pair<string, float>> ovdf;
    for (auto& p : sum_dt)
    {
        vector<string> vec_s = split_str(p.first, ',');
        if (vec_s.size() != 2) continue;

        int N_both = p.second;
        int N_barc1 = count_dict[vec_s[0]];
        int N_barc2 = count_dict[vec_s[1]];
        float jaccard_frag = round((N_both)/(N_barc1+N_barc2-N_both+0.05), 5);
        if (jaccard_frag <= 0.0) continue;
        ovdf.push_back({p.first, jaccard_frag});
    }
    // Sort by jaccard_frag
    std::sort(ovdf.begin(), ovdf.end(), [](pair<string, float>& a, pair<string, float>& b) {
        return a.second > b.second;
    });
    
    // Call knee if we need to
    if (minimum_jaccard_index == 0.0)
    {
        // TODO(fxzhao): implement jaccard option
    }

    // Export the implicated barcodes
    FILE * tbl_out;
    tbl_out = fopen(out_csv_file.c_str(), "w");
    string header = "barc1,barc2,N_both,N_barc1,N_barc2,jaccard_frag,merged\n";
    fwrite(header.c_str(), 1, header.size(), tbl_out);
    for (size_t i = 0; i < ovdf.size(); ++i)
    {
        auto& p = ovdf[i];
        string s = p.first + ",";
        s += to_string(sum_dt[p.first]) + ",";
        vector<string> vec_s = split_str(p.first, ',');
        if (vec_s.size() != 2) continue;
        int N_barc1 = count_dict[vec_s[0]];
        int N_barc2 = count_dict[vec_s[1]];
        s += to_string(N_barc1) + "," + to_string(N_barc2) + ",";
        s += f2str(p.second, 5) + ",";
        s += p.second > minimum_jaccard_index ? "TRUE" : "FALSE";
        s += "\n";

        fwrite(s.c_str(), 1, s.size(), tbl_out);
    }
    fclose(tbl_out);

    // Filter based on the minimum_jaccard_index 
    // and prepare dict data
    map<string, set<string>> bc1, bc2;
    for (size_t i = 0; i < ovdf.size(); ++i)
    {
        auto& p = ovdf[i];
        if (p.second <= minimum_jaccard_index) continue;
        vector<string> vec_s = split_str(p.first, ',');
        if (vec_s.size() != 2) continue;
        bc1[vec_s[0]].insert(vec_s[1]);
        bc2[vec_s[1]].insert(vec_s[0]);
    }

    
    // Guess at how wide we need to make the barcodes to handle leading zeros
    int guess = ceil(log10(nBC.size()));
    // Map from barcode to position in nBC
    map<string, int> bar2pos;
    for (size_t i = 0; i < nBC.size(); ++i)
    {
        bar2pos[nBC[i].first] = i;
    }

    map<string, string> drop_barcodes;

    // Loop through and eat up barcodes
    int idx = 1;
    for (size_t i = 0; i < nBC.size(); ++i)
    {
        string barcode = nBC[i].first;
        if (barcode.empty()) continue;

        set<string> barcode_combine;
        barcode_combine.insert(barcode);

        // Find friends that are similarly implicated and append from Barcode 1
        if (bc1.count(barcode) != 0)
        {
            barcode_combine.insert(bc1[barcode].begin(), bc1[barcode].end());
        }
        // Find friends that are similarly implicated and append from Barcode 2
        if (bc2.count(barcode) != 0)
        {
            barcode_combine.insert(bc2[barcode].begin(), bc2[barcode].end());
        }

        // If user species one to one, then only remove that one barcode
        if (one_to_one) barcode_combine = {barcode};

        // Make a drop barcode and save our progress
        stringstream ss;
        if (!barcoded_tn5)
        {
            ss << name <<"_BC"<<std::setfill ('0') << std::setw (guess) << idx<<
                "_N"<<std::setw (2)<<barcode_combine.size();
        }   
        else
        {
            // TODO(fxzhao): implement this
        }
        for (auto& b : barcode_combine)
        {
            if (bar2pos.count(b) != 0)
            {
                drop_barcodes[b] = ss.str();
                nBC[bar2pos[b]].first = "";
            }
        }
        ++idx;
    }
    
    FILE * bt;
    string replace = "implicatedBarcodes.csv";
    out_csv_file.replace(out_csv_file.find(replace), replace.size(), "barcodeTranslate.tsv");
    bt = fopen(out_csv_file.c_str(), "w");
    for (auto& p : drop_barcodes)
    {
        string s = p.first + "\t" + p.second+"\n";
        fwrite(s.c_str(), 1, s.size(), bt);
    }
    fclose(bt);

    // Finally, collate the NC values
    map<int, int> nc_data;
    for (auto& p : fs::directory_iterator(dir))
    {
        if (p.path().extension() == ".tsv" && p.path().string().find("_ncCount") != std::string::npos)
        {
            cout<<p.path()<<endl;
            ifstream ifs(p.path(), std::ifstream::in);
            string line;
            while (std::getline(ifs, line))
            {
                vector<string> vec_s = split_str(line, '\t');
                if (vec_s.size() != 2) continue;
                nc_data[stoi(vec_s[0])] += stoi(vec_s[1]);
            }
            ifs.close();
        }
    }

    FILE * nc_file_out;
    replace = "barcodeTranslate.tsv";
    out_csv_file.replace(out_csv_file.find(replace), replace.size(), "NCsumstats.tsv");
    nc_file_out = fopen(out_csv_file.c_str(), "w");
    header = "NC_value\tNumberOfFragments\n";
    fwrite(header.c_str(), 1, header.size(), nc_file_out);
    for (auto& p : nc_data)
    {
        string s = to_string(p.first) + "\t" + to_string(p.second) + "\n";
        fwrite(s.c_str(), 1, s.size(), nc_file_out);
    }
    fclose(nc_file_out);


    return 0;
}
