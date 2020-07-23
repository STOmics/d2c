/*
 * File: utility.cpp
 * Created Data: 2020-7-13
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "utility.h"
#include <sstream>
#include <iomanip>

vector<string> split_str(const std::string& str, char delim, bool skip_empty)
{
    std::istringstream iss(str);
    vector<string> res;
    for (std::string item; getline(iss, item, delim);)
        if (skip_empty && item.empty())
            continue;
        else
            res.push_back(item);
    return res;
}

float round(float f, int bits) 
{ 
    stringstream ss; 
    ss << fixed << setprecision(bits) << f; 
    ss >> f; 
    return f;
}

string f2str(float f, int bits)
{
    stringstream ss; 
    ss << fixed << setprecision(bits) << f; 
    return ss.str();
}