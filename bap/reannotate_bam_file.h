/*
 * File: reannotate_bam_file.h
 * Created Data: 2020-7-15
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>
using namespace std;

int reannotate_bam_file(string input_bam_file, string output_bam_file,
    string bead_barcode, string drop_barcode,
    string dict_file, string hq_frags_file);