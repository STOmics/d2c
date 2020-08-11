/*
 * File: utility.h
 * Created Data: 2020-7-13
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
#include <vector>
using namespace std;

vector< string > split_str(const std::string& str, char delim = ' ', bool skip_empty = true);

float round(float src, int bits);

string f2str(float f, int bits);

// For executing shell script
int exec_shell(const char* cmd, std::vector< std::string >& resvec);

// Get physical memory used by process in real time
// Only support linux, unit is MB
size_t physical_memory_used_by_process();