/*
 * File: utility.h
 * Created Data: 2020-7-13
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <vector>
#include <string>
using namespace std;

vector<string> split_str(const std::string& str, char delim=' ', bool skip_empty=true);

float round(float src, int bits);

string f2str(float f, int bits);