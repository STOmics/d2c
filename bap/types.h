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
    int start;
    int end;
    // string qname;
    int qname1;
    int qname2;
    // string barcode;
    int barcode;
};