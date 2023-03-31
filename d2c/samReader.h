/*
 * File: samReader.h
 * Created Data: 2020-5-12
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>

typedef bam1_t*    BamRecord;
typedef bam_hdr_t* BamHeader;

class SamReader
{
public:
    // Creates a new SamReader reading from the BAM file reads_path
    static std::unique_ptr< SamReader > FromFile(const std::string& reads_path);

    ~SamReader();

    // DIsable assignment/copy operations.
    SamReader(const SamReader& other) = delete;
    SamReader& operator=(const SamReader&) = delete;

    // Close the underlying resource descriptors.
    int Close();

    // Return <chromosome,index> pairs.
    std::vector< std::pair< std::string, uint32_t > > getContigs();

    // Return header for construct other instance.
    BamHeader getHeader();

    // Query records by a given contig.
    bool QueryByContig(int tid);

    // Iterate all records from query range.
    // NOTICE: must call next() after call QueryByContig()
    bool next(BamRecord b);
    bool next(BamRecord b, hts_itr_t*& iter);
    bool next(BamRecord b, int flag);
    int  QueryOne(BamRecord b);

    void setThreadPool(htsThreadPool* p);

private:
    // Private constructor; use FromFile to safely create a SamReader from a
    // file.
    SamReader(htsFile* fp, bam_hdr_t* header, hts_idx_t* idx);

    // A pointer to the htslib file used to access the SAM/BAM data.
    htsFile* fp_;

    // A htslib header data structure obtained by parsing the header of this BAM.
    bam_hdr_t* header_;

    // The htslib index data structure for our indexed BAM file.
    // May be NULL if no index was loaded.
    hts_idx_t* idx_;

    // Store reference name and length from header.
    std::vector< std::pair< std::string, uint32_t > > ref_;

    // A pointer to the query result.
    hts_itr_t* iter_;
};
