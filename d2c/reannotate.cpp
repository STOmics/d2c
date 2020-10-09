/*
 * File: reannotate.cpp
 * Created Data: 2020-10-9
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "reannotate.h"
#include "utility.h"
#include "samReader.h"

#include <unordered_map>
#include <fstream>
#include <iostream>
#include <thread>

#include <htslib/thread_pool.h>

extern bool getTag(BamRecord b, const char tag[2], std::string& value);

typedef unordered_map<string, string> str_map;
constexpr auto drop_tag = "DB";

// Trim barcode sequence by seperator
inline string trim_barcode(string& barcode_seq, char sep='-')
{
    auto pos = barcode_seq.find(sep);
    return pos == string::npos ? barcode_seq : barcode_seq.substr(0, pos);
}

str_map get_translate_map(string& barcode_translate_file)
{
    str_map m;
    ifstream ifs(barcode_translate_file);
    string line;
    while (std::getline(ifs, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() != 2) continue;
        m.insert({trim_barcode(vec_s[0]), vec_s[1]});
    }
    return m;
}

bool reannotate(string input_bam, string barcode_translate_file, string output_bam, string barcode_tag)
{
    // Prepare the translate map from cell barcode to drop barcode
    str_map barcode_map = get_translate_map(barcode_translate_file);
    if (barcode_map.empty())
    {
        cerr<<"Please check the format of barcode translate file: "<<barcode_translate_file<<endl;
        return false;
    }

    // Travel over the input bam file, and search the barcode in translate map    
    htsFile* in_fp = hts_open(input_bam.c_str(), "r");
    bam_hdr_t* header = sam_hdr_read(in_fp);

    htsFile* out_fp = hts_open(output_bam.c_str(), "wb");
    [[maybe_unused]] int ret = sam_hdr_write(out_fp, header);


    // Create and share the thread pool
	htsThreadPool p = {NULL, 0};
	const int THREADS_NUM = std::thread::hardware_concurrency();
    //const int THREADS_NUM = 40;
	if (THREADS_NUM > 0)
	{
		p.pool = hts_tpool_init(THREADS_NUM);
		if (p.pool)
		{
			//samReader->setThreadPool(&p);
			//samWriter.setThreadPool(&p);
            hts_set_opt(in_fp, HTS_OPT_THREAD_POOL, &p);
            hts_set_opt(out_fp, HTS_OPT_THREAD_POOL, &p);
		}
	}

    bam1_t* b = bam_init1();
    string bead_bc, drop_bc;
    str_map::iterator iter;
    while (sam_read1(in_fp, header, b) >= 0)
    {
        if (getTag(b, barcode_tag.c_str(), bead_bc))
        {
            bead_bc = trim_barcode(bead_bc);
            if ((iter = barcode_map.find(bead_bc)) != barcode_map.end())
            {
                drop_bc = iter->second;
                bam_aux_append(b, drop_tag, 'Z', drop_bc.size() + 1, ( uint8_t* )drop_bc.c_str());
            }
        }
        
        [[maybe_unused]] int ret = sam_write1(out_fp, header, b);
    }

    hts_close(in_fp);
    in_fp = nullptr;
    bam_hdr_destroy(header);
    header = nullptr;
    hts_close(out_fp);
    out_fp = nullptr;
    bam_destroy1(b);
    if (p.pool)
	{
		hts_tpool_destroy(p.pool);
	}

    return true;
}