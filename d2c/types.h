/*
 * File: types.h
 * Created Data: 2020-8-7
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

struct Bedpe
{
    int start;  // start and end positions for this fragment
    int end;
    
    int qname1; // record reads position for keep reads in bam
    int qname2;
    
    int barcode; // encode runname,bc1,bc2 to int32: 8+12+12
};