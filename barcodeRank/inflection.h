/*
 * File: inflection.h
 * Created Data: 2020-9-30
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
#include <vector>

using namespace std;

// Run-length encode
void rle(vector< double >& vec, double thre, vector< pair< double, double > >& encodes);

double inflection(vector< double >& vec, double lower = 500, double exclude_from = 50);
