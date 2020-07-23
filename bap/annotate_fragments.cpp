/*
 * File: annotate_fragments.cpp
 * Created Data: 2020-7-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "annotate_fragments.h"
#include "gz_io.h"

#include <fstream>
#include <map>
#include <iostream>
#include <vector>
#include <sstream>

#include "ygg.hpp"

#include "utility.h"

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

int annotate_fragments(string blacklist_file, string frag_file, string read_bead_file,
    string annotated_file, string stat_file)
{
    // Load read bead file and make dict
    gzFile rbf = gzopen(read_bead_file.c_str(), "rb");
    string line;
    map<string, string> dict;
    while (readline(rbf, line))
    {
        size_t pos = line.find("\t");
        if (pos == std::string::npos) continue;
        dict[line.substr(0, pos)] = line.substr(pos+1);
        //cout<<line.substr(0, pos)<<" "<< line.substr(pos+1)<<endl;
    }
    gzclose(rbf);

    // Load blacklist file and construct interval tree
    ifstream blf(blacklist_file, std::ifstream::in);
    vector<Node> nodes;
    while (std::getline(blf, line))
    {
        vector<string> vec_str = split_str(line, '\t');
        if (vec_str.size() != 3) continue;
        string chr = vec_str[0];
        int start = stoi(vec_str[1]);
        int end = stoi(vec_str[2]);
        Node node;
        node.lower = start;
        node.upper = end;
        node.value = chr;
        //cout<<node.value<<" "<<node.lower<<" "<<node.upper<<endl;
        nodes.push_back(std::move(node));
    }
    blf.close();

    map<string, MyTree> mytrees;
    for (auto& node : nodes)
    {
        //cout<<node.value<<" "<<node.lower<<" "<<node.upper<<endl;
        mytrees[node.value].insert(node);
    }

    // Load fragment file
    //ifstream ifs(frag_file, std::ifstream::in);
    gzFile ifs = gzopen(frag_file.c_str(), "rb");
    //ofstream ofs(annotated_file, std::ofstream::out);
    FILE * ofs;
    ofs = fopen(annotated_file.c_str(), "w");
    set<string> pcr_dup;
    while (readline(ifs, line))
    {
        vector<string> vec_str = split_str(line, '\t');
        if (vec_str.size() != 4) continue;
        string name = vec_str[3];
        //cout<<name<<endl;
        if (dict.count(name) != 0)
        {
            // Filter for fragments overlapping the blacklist
            string chr = vec_str[0];
            int start = stoi(vec_str[1]);
            int end = stoi(vec_str[2]);
            MyInterval query_range {start, end};
            if (mytrees.count(chr) != 0)
            {
                const auto& tree = mytrees.at(chr);
                const auto& res = tree.query(query_range);
                //cout<<start<<" "<<end<<endl;
                // Discard the fragments overlapping the blacklist
                if (res.begin() != res.end()) continue;
            }

            line += "\t" + dict[name] + "\n";
            //ofs << line <<endl;
            fwrite(line.c_str(), 1, line.size(), ofs);

            string other = vec_str[0]+"\t"+vec_str[1]+"\t"+vec_str[2]+"\t"+dict[name];
            pcr_dup.insert(other);
        }
    }
    //ifs.close();
    gzclose(ifs);
    //ofs.close();
    fclose(ofs);

    // Quantify the number of unique fragments per barcode
    map<string, int> bead_quant;
    for (auto& l : pcr_dup)
    {
        size_t pos = l.find_last_of("\t");
        if (pos == std::string::npos) continue;
        string bead_id = l.substr(pos+1);
        ++bead_quant[bead_id];
    }
    FILE * out_bead_quant;
    out_bead_quant = fopen(stat_file.c_str(), "w");
    for (auto& b : bead_quant)
    {
        string s = b.first + "\t" + to_string(b.second) + "\n";
        fwrite(s.c_str(), 1, s.size(), out_bead_quant);
    }
    fclose(out_bead_quant);

    return 0;
}