/*
 * File: barcodeRank.cpp
 * Created Data: 2020-9-30
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "barcodeRank.h"
#include "kde.hpp"
#include "inflection.h"

#include <exception>

double barcode_rank(vector< double >& input, INFLECTION_KERNEL_TYPE inflection_kernel_type, CURVE_DATA_TYPE curve_data_type)
{
    if (inflection_kernel_type == INFLECTION_KERNEL_TYPE::BAP2)
    {
        KDE  kde;
        pair<double, double> paras;
        if (curve_data_type == CURVE_DATA_TYPE::BEAD)
            paras            = kde.run(input, "bead");
        else if (curve_data_type == CURVE_DATA_TYPE::JACCARD)
            paras            = kde.run(input, "jaccard");
        else
            throw std::runtime_error("Error curve data type: " + to_string(curve_data_type));
        return paras.first;
    }
    else if (inflection_kernel_type == INFLECTION_KERNEL_TYPE::DROPLETUTILS)
    {
        if (curve_data_type == CURVE_DATA_TYPE::BEAD)
            return inflection(input, 500);
        else if (curve_data_type == CURVE_DATA_TYPE::JACCARD)
            return inflection(input, 0.005);
        else 
            throw std::runtime_error("Error curve data type: " + to_string(curve_data_type));
    }
    else
    {
        throw std::runtime_error("Error inflection kernel type: " + to_string(inflection_kernel_type));
    }
}