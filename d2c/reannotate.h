/*
 * File: reannotate.h
 * Created Data: 2020-10-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <string>

using namespace std;

bool reannotate(string input_bam, string barcode_translate_file, string output_bam, string barcode_tag);