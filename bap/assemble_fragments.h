/*
 * File: assemble_fragments.h
 * Created Data: 2020-7-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
using namespace std;

// Extract fragments from bam file.
// Using methods from samtools/bedtools .etc
int assemble_fragments(const string bam_file, const string frag_file);