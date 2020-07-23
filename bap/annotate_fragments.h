/*
 * File: annotate_fragments.h
 * Created Data: 2020-7-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
using namespace std;

// Merge read bead file to frag file and 
// find overlap reads with blacklist file, then remove
// these reads
int annotate_fragments(string blacklist_file, string frag_file, string read_bead_file,
    string annotated_file, string stat_file);
