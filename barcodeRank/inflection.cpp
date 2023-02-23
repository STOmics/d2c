/*
 * File: inflection.cpp
 * Created Data: 2020-9-30
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "inflection.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

#include <cmath>

double inflection(vector< double >& vec, double lower, double exclude_from, bool cut)
{
    exclude_from = log10(exclude_from);

    // Sort in descending order
    vector< double > order(vec.begin(), vec.end());
    std::sort(order.begin(), order.end(), std::greater< double >());

    // Run-length encoding
    // <lengths, values>
    vector< pair< double, double > > stuff;
    rle(order, lower, stuff);

    // Get mid-rank of each run
    size_t cumsum = 0;
    std::for_each(stuff.begin(), stuff.end(), [&](pair< double, double >& p) {
        auto& [l, v] = p;
        cumsum += int(l);
        l = cumsum - (l - 1) / 2;
        l = log10(l);
        v = log10(v);
    });

    // Get the diff of the next element and current element
    vector< pair< double, double > > diff(stuff.size());
    std::adjacent_difference(stuff.begin(), stuff.end(), diff.begin(),
                             [](const pair< double, double >& a, const pair< double, double >& b) {
                                 return std::make_pair(a.first - b.first, a.second - b.second);
                             });

    // Get the diff ratio
    vector< double > d1n(diff.size());
    std::transform(diff.begin(), diff.end(), d1n.begin(),
                   [](const pair< double, double >& p) { return p.second / p.first; });

    // Some exclusion of the LHS points avoid problems with discreteness
    int count = std::count_if(stuff.begin(), stuff.end(),
                              [&](const pair< double, double >& p) { return (p.first <= exclude_from); });
    int skip  = min(int(d1n.size() - 1), count);

    
    int min_pos = -1;
    if (cut)
    {
        double min_val = 65535.0;
        for (int i = skip; i < d1n.size(); ++i)
        {
            if (min_val > fabs(d1n[i]+1))
            {
                min_val = fabs(d1n[i]+1);
                min_pos = i+1;
            }
        }
    }
    else
    {
        auto min_iter = std::min_element(d1n.begin() + skip, d1n.end());
        // The first element in diff/d1n is empty
        min_pos = std::distance(d1n.begin(), min_iter) - 1;
    }

    double inflection = 0;
    if (min_pos >= 0 && min_pos < static_cast< int >(stuff.size()))
        inflection = pow(10, stuff[min_pos].second);

    return inflection;
}

void rle(vector< double >& vec, double thre, vector< pair< double, double > >& encodes)
{
    for (size_t i = 0; i < vec.size(); ++i)
    {
        int count = 1;
        while (i < vec.size() && vec[i] == vec[i + 1])
        {
            count++;
            i++;
        }
        if (vec[i] <= thre)
            continue;
        encodes.push_back({ count, vec[i] });
    }
}
