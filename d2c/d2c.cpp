/*
 * File: d2c.cpp
 * Created Data: 2020-7-23
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "d2c.h"
#include "bamCat.h"
#include "gzIO.h"
#include "timer.h"
#include "utility.h"

#include <sparsepp/spp.h>
#include <spdlog/spdlog.h>
#include <taskflow/taskflow.hpp>

#include <algorithm>
#include <exception>
#include <fstream>
#include <thread>
#include <tuple>

#include <htslib/bgzf.h>
#include <htslib/tbx.h>

// Control the standart of uniq fragments, BOTH mean start and end are both equal
#define UNIQ_FRAG_BOTH

// Samtools view's parameter
const int FLAG       = 2;
const int MAX_INSERT = 2000;

// Bead barcode format: bc1bc2-runname
// the length of bc1 or bc2 is 10bp, and their species are 1536,
// so encode bc1,bc2 and runname to int32: 8+12+12
constexpr int BBIT  = 12;        // barcode length in bit
const char*   BSEP  = "-";       // barcode separator
constexpr int BLEN  = 10;        // single barcode length
constexpr int BBLEN = BLEN * 2;  // total barcode length
constexpr int BMASK = 0xFFF;     // mask for get single bacode value
constexpr int RMASK = 0xFF;      // mask for get runname value

// Filenames
constexpr auto PLOT_SCRIPT = "plot.pyc";
constexpr auto R_SCRIPT    = "rank.R";
// constexpr auto BEAD_THRE_SCRIPT         = "10b_knee_execute.R";
// constexpr auto JACCARD_THRE_SCRIPT      = "11b_knee_execute.R";
constexpr auto PARAM_FILE    = ".d2cCutoff.tsv";
constexpr auto SAT_FILE      = ".sequenceSaturation.tsv";
constexpr auto FRAGMENT_FILE = ".fragments.tsv.gz";
// constexpr auto BASIC_QC_FILE            = ".basicQC.tsv";
constexpr auto BARCODE_QUANT_FILE = ".barcodeCount.tsv";
// constexpr auto HQ_BEADS_FILE            = ".HQbeads.tsv";
// constexpr auto JACCARD_TMP_FILE         = ".jaccard.csv";
constexpr auto IMPLICATED_BARCODES_FILE = ".CorrelationBarcodes.tsv.gz";
constexpr auto TEMP_JACCARD_FILE        = ".tempJaccard.tsv";
constexpr auto BARCODE_TRANSLATE_FILE   = ".barcodeMerge.tsv";
// constexpr auto NC_STATS_FILE            = ".NCsumstats.tsv";
constexpr auto QC_STATS_FILE = ".Metadata.tsv";

constexpr auto FSEP = '\t';  // seperator of lines in output file

#include "ygg.hpp"

using MyInterval = std::pair< int, int >;
template < class Node > class NodeTraits : public ygg::ITreeNodeTraits< Node >
{
public:
    using key_type = int;
    static int get_lower(const Node& node)
    {
        return node.lower;
    }
    static int get_upper(const Node& node)
    {
        return node.upper;
    }
    static int get_lower(const MyInterval& i)
    {
        return std::get< 0 >(i);
    }
    static int get_upper(const MyInterval& i)
    {
        return std::get< 1 >(i);
    }
};

class Node : public ygg::ITreeNodeBase< Node, NodeTraits< Node > >
{
public:
    int    upper;
    int    lower;
    string value;
};

using MyTree = ygg::IntervalTree< Node, NodeTraits< Node > >;

bool getTag(BamRecord b, const char tag[2], std::string& value)
{
    uint8_t* data = bam_aux_get(b, tag);
    if (data == NULL)
        return false;
    value = bam_aux2Z(data);
    return true;
}

inline string getQname(BamRecord b)
{
    return bam_get_qname(b);
}

inline bool isPaired(BamRecord b)
{
    return b->core.flag & BAM_FPAIRED;
}

inline bool isMapped(BamRecord b)
{
    return !(b->core.flag & BAM_FUNMAP);
}

inline bool isReverse(BamRecord b)
{
    return b->core.flag & BAM_FREVERSE;
}
inline bool isRead1(BamRecord b)
{
    return b->core.flag & BAM_FREAD1;
}

struct SummaryData
{
    int           overlaps;        // The number of has overlap with TSS
    float         tss_proportion;  // overlaps / total
    float         mean_insert_size;
    int           median_insert_size;
    float         frip;  // The mean value of peak
    vector< int > insert_size;
};

void decodeSpeciesMix(unsigned long l, int& totalCnt, int& humanCnt, int& mouseCnt)
{
    humanCnt = l & 0xFFFFFFFF;
    mouseCnt = (l >> 32) & 0xFFFFFFFF;
    totalCnt = humanCnt + mouseCnt;
}

float get_dup_proportion(int t, int u)
{
    return (t - u) * 1.0 / t;
}

inline float aux(float x, int c, int n)
{
    return (c / x - 1 + exp(-n / x));
}
int get_library_size(int t, int u)
{
    float m     = 1;
    float M     = 100;
    int   n_dup = t - u + 1;
    if (u > t || aux(m * u, u, t) < 0 || u < 0 || t < 0 || n_dup < 0)
        return 0;

    while (aux(M * u, u, t) > 0)
        M *= 10.0;

    for (int i = 0; i <= 40; ++i)
    {
        float mid = (m + M) / 2;
        float v   = aux(mid * u, u, t);
        if (v == 0)
            break;
        else if (v > 0)
            m = mid;
        else
            M = mid;
    }
    return round(u * (m + M) / 2.0);
}

void D2C::extractBedPE(const BamRecord b1, const BamRecord b2, vector< Bedpe >& bedpes, int l1, int l2)
{
    // Initialize BEDPE variables
    string chrom1, chrom2, strand1, strand2;
    int    start1, start2, end1, end2;
    start1 = start2 = end1 = end2 = -1;
    chrom1 = chrom2 = strand1 = strand2 = ".";
    int min_map_quality = 0, max_map_quality = 0;

    // Extract relevant info for end 1
    if (isMapped(b1))
    {
        // chrom1 = refs[b1->core.tid];
        start1 = b1->core.pos;
        end1   = bam_endpos(b1);
        if (isReverse(b1))
            strand1 = "-";
        else
            strand1 = "+";
    }

    // Extract relevant info for end 2
    if (isMapped(b2))
    {
        // chrom2 = refs[b2->core.tid];
        start2 = b2->core.pos;
        end2   = bam_endpos(b2);
        if (isReverse(b2))
            strand2 = "-";
        else
            strand2 = "+";
    }

    // Swap the ends if necessary
    // if (chrom1 > chrom2 || ((chrom1 == chrom2) && (start1 > start2)))
    // If the starts of two record are euqal, use read1's strand
    if (start1 > start2 || (start1 == start2 && isRead1(b2)))
    {
        swap(chrom1, chrom2);
        swap(start1, start2);
        swap(end1, end2);
        swap(strand1, strand2);
    }

    // Report BEDPE using min mapping quality
    if (isMapped(b1) && isMapped(b2))
    {
        min_map_quality = min(b1->core.qual, b2->core.qual);
        max_map_quality = max(b1->core.qual, b2->core.qual);
    }

    // Filter by max insert size and mapping quality
    if ((end2 - start1) < MAX_INSERT && min_map_quality >= mapq && max_map_quality > mapq)
    {
        if (strand1 == "+")
        {
            start1 += 4;
            end2 += 4;
        }
        else if (strand1 == "-")
        {
            start1 -= 5;
            end2 -= 5;
        }
        // ofs << chrom1 <<"\t"<<start1<<"\t"<<end2<<"\t"<<get_qame(b1)<<endl;
        string barcode("NA");
        if (getTag(b1, barcode_tag.c_str(), barcode))
        {
            Bedpe bedpe;
            bedpe.start = start1;
            bedpe.end   = end2;
            // bedpe.qname = getQname(b1);
            bedpe.qname1 = l1;
            bedpe.qname2 = l2;

            string b1 = barcode.substr(0, BLEN);
            string b2 = barcode.substr(BLEN, BLEN);
            if (_barcode2int.count(b1) == 0 || _barcode2int.count(b2) == 0)
            {
                spdlog::warn("Invalid barcode {} not exists in barcode list", barcode);
                return;
            }
            // There is not exists tn5 sequence in bam format
            bedpe.barcode = (_barcode2int[b1] << BBIT) + _barcode2int[b2];

            // cout<<b1<<" "<<b2<<" "<<bedpe.barcode<<endl;
            bedpes.push_back(bedpe);
        }
    }
}

D2C::D2C(string input_bam, string output_path, string barcode_in_tag, string barcode_out_tag, int mapq, int cores,
         string run_name, int tn5, double min_barcode_frags, double min_jaccard_index, string ref, string mito_chr,
         string bed_genome_file, string blacklist_file, string trans_file, bool species_mix, string bin_path,
         int barcode_threshold, int jaccard_threshold, bool saturation_on, string barcode_list,
         string barcode_runname_list, int beads_force, string tn5_list, string rank)
    : input_bam(input_bam), output_path(output_path), barcode_tag(barcode_in_tag), drop_tag(barcode_out_tag),
      mapq(mapq), cores(cores), run_name(run_name), tn5(tn5), min_barcode_frags(min_barcode_frags),
      min_jaccard_index(min_jaccard_index), ref(ref), bed_genome_file(bed_genome_file), blacklist_file(blacklist_file),
      trans_file(trans_file), species_mix(species_mix), bin_path(bin_path), barcode_threshold(barcode_threshold),
      jaccard_threshold(jaccard_threshold), saturation_on(saturation_on), barcode_list(barcode_list),
      barcode_runname_list(barcode_runname_list), beads_force(beads_force), tn5_list(tn5_list), rank(rank)
{
    nc_threshold         = 6;
    regularize_threshold = 4;
    one_to_one           = false;

    if (fs::path(input_bam).extension() == ".bed")
        _is_bed = true;
    else
        _is_bed = false;

    auto pos = mito_chr.find(',');
    if (pos == std::string::npos)
    {
        single_mc = mito_chr;
        mito_chrs.insert(single_mc);
    }
    else
    {
        string chr1 = mito_chr.substr(0, pos);
        string chr2 = mito_chr.substr(pos + 1);
        if (chr1.at(0) == 'h' && chr2.at(0) == 'm')
        {
            human_mc = chr1;
            mouse_mc = chr2;
        }
        else if (chr1.at(0) == 'm' && chr2.at(0) == 'h')
        {
            human_mc = chr2;
            mouse_mc = chr1;
        }
        mito_chrs.insert(human_mc);
        mito_chrs.insert(mouse_mc);
    }
}

int D2C::run()
{
    // Prepare
    // Make sure there is index of bam file
    Timer timer;
    temp_bam_path = output_path / ".temp_bam";
    if (!fs::exists(temp_bam_path))
    {
        if (!fs::create_directories(temp_bam_path))
        {
            spdlog::warn("Failed to create directory: {}", temp_bam_path.string());
            return -1;
        }
    }
    // Parse barcode list
    if (!parseBarcodeList())
    {
        spdlog::warn("Failed to parse barcode list: {}", barcode_list);
        return -2;
    }
    // Detect barcode runname list
    if (!parseRunnameList())
    {
        spdlog::warn("Failed to parse runname list: {}", barcode_runname_list);
        return -2;
    }
    // Parse tn5 list
    parseTn5List();

    // Check is make sense for mix_species and mito chrs
    if ((species_mix && (human_mc.empty() || mouse_mc.empty()))
        || (!species_mix && !human_mc.empty() && !mouse_mc.empty()))
    {
        spdlog::warn("species_mix is {} but human_mc is {} and mouse_mc is {}", species_mix, human_mc, mouse_mc);
        return -2;
    }

    spdlog::info("Prepare time(s): {:.2f}", timer.toc(1000));

    // Execute taskflow
    taskflow();
    spdlog::info("Execute taskflow time(s): {:.2f}", timer.toc(1000));

    // Clean up
    fs::remove_all(temp_bam_path);
    spdlog::info("Clean up time(s): {:.2f}", timer.toc(1000));

    return 0;
}

int D2C::taskflow()
{
    tf::Executor executor(cores);
    tf::Taskflow taskflow;

    taskflow.name("D2C");

    // Step 1: split bam by chr
    std::vector< std::pair< std::string, unsigned int > > contigs;
    if (_is_bed)
    {
        auto temp_map = parseChrsFromBedFile(input_bam);
        for (auto& [k, _] : temp_map)
            contigs.push_back({ k, 0 });
    }
    else
    {
        auto samReader = SamReader::FromFile(input_bam);
        contigs        = samReader->getContigs();
        BamRecord b    = bam_init1();
        samReader->QueryOne(b);
        string qname = getQname(b);
        spdlog::debug("qname: {}", qname);

        string suffix = qname.substr(qname.size() - 2);
        if (suffix == "/1" || suffix == "/2")
        {
            spdlog::debug("Bam has suffix: {}", suffix);
            bam_has_suffix = true;
        }
        bam_destroy1(b);
    }
    spdlog::debug("Bam contigs num: {}", contigs.size());
    for (auto& p : contigs)
    {
        _contig_names.push_back(p.first);
        spdlog::debug("bam contig: {}", p.first);
    }

    // Verify that the supplied reference genome and the bam have overlapping chromosomes
    auto bed_chrs = parseChrsFromBedFile(bed_genome_file);
    if (bed_chrs.empty())
    {
        for (auto& ctg : _contig_names)
            bed_chrs.emplace(ctg, "");
    }
    spdlog::debug("bed_chrs num: {}", bed_chrs.size());
    for (auto& [contig, _] : bed_chrs)
        spdlog::debug("bed contig: {}", contig);

    vector< int > used_chrs;
    for (size_t i = 0; i < contigs.size(); ++i)
    {
        if (bed_chrs.count(_contig_names[i]) != 0)
            used_chrs.push_back(i);
    }
    if (used_chrs.empty())
    {
        spdlog::error("Found no overlapping chromosomes between bam and reference. Check reference genome "
                      "specification with the -r flag");
        return -1;
    }
    else
    {
        spdlog::info("Found {} chromosomes for analysis (including mitochondria)", used_chrs.size());
    }

    spdlog::debug("Initial memory(MB): {}", physical_memory_used_by_process());
    _bedpes_by_chr.resize(contigs.size());
    if (_is_bed)
    {
        map< string, int > chr_to_int;
        for (auto& id : used_chrs)
            chr_to_int[_contig_names[id]] = id;

        // Parse data for bed file
        ifstream ifs(input_bam, std::ifstream::in);
        string   line;
        while (std::getline(ifs, line))
        {
            vector< string > vec_s = split_str(line, '\t');
            if (vec_s.size() < 5)
                continue;
            if (chr_to_int.count(vec_s[0]) == 0)
                continue;
            string& barcode = vec_s[3];
            string  b1      = barcode.substr(0, BLEN);
            string  b2      = barcode.substr(BLEN, BLEN);
            if (_barcode2int.count(b1) == 0 || _barcode2int.count(b2) == 0)
            {
                spdlog::warn("Invalid barcode {} not exists in barcode list", barcode);
                continue;
            }

            Bedpe bedpe;
            bedpe.start  = stoi(vec_s[1]);
            bedpe.end    = stoi(vec_s[2]);
            bedpe.qname1 = 0;
            bedpe.qname2 = 0;

            int tail = 0;
            if (barcode.size() > BBLEN)
            {
                string tmp = barcode.substr(BBLEN);
                if (_tail2int.count(tmp) == 0)
                {
                    spdlog::warn("Invalid tn5 {} not exists in tn5.list", barcode);
                    continue;
                }
                tail = _tail2int[tmp];
            }
            bedpe.barcode = (tail << BBIT * 2) + (_barcode2int[b1] << BBIT) + _barcode2int[b2];

            auto& bedpes = _bedpes_by_chr[chr_to_int[vec_s[0]]];
            bedpes.insert(bedpes.end(), stoi(vec_s[4]), bedpe);
        }
        ifs.close();
    }
    auto [S, T] = taskflow.parallel_for(used_chrs.begin(), used_chrs.end(),
                                        [&](int chr_id)
                                        {
                                            Timer t;
                                            D2C::splitBamByChr(chr_id);
                                            spdlog::info("Split bam by chr: {} time(s): {:.2f}", contigs[chr_id].first,
                                                         t.toc(1000));
                                        });

    S.name("Start SplitBam/Assemble/Annotate");
    T.name("End SplitBam/Assemble/Annotate");
    S.for_each_successor([contigs, s = 0](tf::Task successor) mutable
                         { successor.name((contigs.begin() + (s++))->first); });

    // Step 2: determine hight quality beads
    auto determine_hq_beads = taskflow
                                  .emplace(
                                      [&]()
                                      {
                                          spdlog::info("SplitBam memory(MB): {}", physical_memory_used_by_process());
                                          Timer t;
                                          D2C::determineHQBeads();
                                          spdlog::info("Determine high-quality beads time(s): {:.2f}", t.toc(1000));
                                      })
                                  .name("Determine HQ Beads");
    determine_hq_beads.succeed(T);

    // Step 3: calculate stat by chr
    _total_bead_cnts.resize(contigs.size());
    _total_nc_cnts.resize(contigs.size());
    // cout<<"contig size: "<<_bedpes_by_chr.size()<<endl;
    // for (int i = 0; i < _bedpes_by_chr.size(); ++i)
    //     if (_bedpes_by_chr[i].size() > 0)
    //         cout<<"has data chr: "<<i<<endl;
    auto [start_cal, end_cal] = taskflow.parallel_for(used_chrs.begin(), used_chrs.end(),
                                                      [&](int chr_id)
                                                      {
                                                          // Skip the mito chrom
                                                          if (mito_chrs.count(_contig_names[chr_id]))
                                                              return;
                                                          Timer t;
                                                          D2C::computeStatByChr(chr_id);
                                                          spdlog::info("Compute stat by chr: {} time(s): {:.2f}",
                                                                       _contig_names[chr_id], t.toc(1000));
                                                      });
    start_cal.name("Start Cal");
    end_cal.name("End Cal");
    start_cal.for_each_successor([contigs, s = 0](tf::Task successor) mutable
                                 { successor.name((contigs.begin() + (s++))->first); });
    start_cal.succeed(determine_hq_beads);
    start_cal.succeed(T);

    // Step 4: barcode merge
    auto barcode_merge = taskflow
                             .emplace(
                                 [&]()
                                 {
                                     spdlog::info("ComputeByChr memory(MB): {}", physical_memory_used_by_process());
                                     Timer t;
                                     D2C::determineBarcodeMerge();
                                     spdlog::info("Determine barcode merge time(s): {:.2f}", t.toc(1000));
                                 })
                             .name("Determine Barcode Merge");
    barcode_merge.succeed(end_cal);
    barcode_merge.succeed(determine_hq_beads);

    // Step 5: reannotate frag data and get summary stats by chr
    _keep_qnames.resize(contigs.size());
    _dup_frags.resize(contigs.size());
    _frag_stats.resize(contigs.size());
    auto [start_reanno, end_reanno] = taskflow.parallel_for(
        used_chrs.begin(), used_chrs.end(),
        [&](int chr_id)
        {
            Timer t;
            D2C::reannotateFragByChr(chr_id);
            spdlog::info("Reannotate frags by chr: {} time(s): {:.2f}", _contig_names[chr_id], t.toc(1000));
        });
    start_reanno.name("Start Reannotate");
    end_reanno.name("End Reannotate");
    start_reanno.for_each_successor([contigs, s = 0](tf::Task successor) mutable
                                    { successor.name((contigs.begin() + (s++))->first); });
    start_reanno.succeed(barcode_merge);
    start_reanno.succeed(T);

    // Step 5.1: sequencing saturation
    auto sequence_saturation = taskflow
                                   .emplace(
                                       [&]()
                                       {
                                           if (saturation_on)
                                           {
                                               Timer    t;
                                               fs::path sat_out_file = output_path / (run_name + SAT_FILE);
                                               saturation.calculateSaturation(sat_out_file.string());
                                               spdlog::info("Sequencing saturation time(s): {:.2f}", t.toc(1000));
                                           }
                                       })
                                   .name("Sequence saturation");
    sequence_saturation.succeed(end_reanno);

    // Step 6: annotate bam file by chr
    auto [start_annobam, end_annobam] = taskflow.parallel_for(
        used_chrs.begin(), used_chrs.end(),
        [&](int chr_id)
        {
            // Skip the mito chrom
            // if (_contig_names[chr_id] == mito_chr)
            //    return;
            if (_is_bed)
                return;

            Timer t;
            D2C::annotateBamByChr(chr_id);
            spdlog::info("Annotate bam file by chr: {} time(s): {:.2f}", _contig_names[chr_id], t.toc(1000));
        });
    start_annobam.name("Start Annotate bam");
    end_annobam.name("End Annotate bam");
    start_annobam.for_each_successor([contigs, s = 0](tf::Task successor) mutable
                                     { successor.name((contigs.begin() + (s++))->first); });
    start_annobam.succeed(end_reanno);
    start_annobam.succeed(T);

    // Step 7: merge bam files
    auto merge_bam = taskflow
                         .emplace(
                             [&]()
                             {
                                 if (_is_bed)
                                     return;

                                 Timer t;
                                 spdlog::debug("Merge bam");
                                 std::vector< std::string > bam_files;
                                 for (auto& chr_id : used_chrs)
                                 {
                                     // Skip the mito chrom
                                     // if (_contig_names[chr_id] == mito_chr)
                                     //    continue;
                                     fs::path tmp_bam_file = temp_bam_path / (contigs[chr_id].first + ".bam");
                                     if (fs::exists(tmp_bam_file))
                                         bam_files.push_back(tmp_bam_file.string());
                                 }
                                 fs::path output_bam_file = output_path / (run_name + ".bam");
                                 spdlog::info("Be merged file size: {}", bam_files.size());
                                 // for (auto& name : bam_files)
                                 //     spdlog::debug("merge bam: {}", name);
                                 try
                                 {
                                     int cmd_rtn = bam_cat(bam_files, nullptr, output_bam_file.c_str(), nullptr, 0);
                                     if (cmd_rtn == 0)
                                     {
                                         spdlog::info("Merge bam file success");
                                         // Build index file
                                         Timer index_timer;
                                         auto  bai_file(output_bam_file);
                                         bai_file += ".bai";
                                         if (fs::exists(bai_file))
                                             fs::remove(bai_file);
                                         auto samReader = SamReader::FromFile(output_bam_file);
                                         spdlog::info("Build bam index time(s): {:.2f}", index_timer.toc(1000));
                                     }
                                     else
                                         spdlog::error("Merge bam file fail, rtn:{}", cmd_rtn);
                                 }
                                 catch (std::exception& e)
                                 {
                                     spdlog::error("Error in merge bam: {}", e.what());
                                 }
                                 catch (...)
                                 {
                                     spdlog::error("Error in merge bam");
                                 }
                                 spdlog::info("Merge bam time(s): {:.2f}", t.toc(1000));
                             })
                         .name("Merge bam");
    merge_bam.succeed(end_annobam);

    // Step 8: merge fragment files
    auto merge_frags = taskflow
                           .emplace(
                               [&]()
                               {
                                   spdlog::info("Reannotate frags memory(MB): {}", physical_memory_used_by_process());
                                   spdlog::debug("Merge frags");
                                   Timer t;
                                   for (auto& chr_id : used_chrs)
                                   {
                                       // Skip the mito chrom
                                       // if (_contig_names[chr_id] == mito_chr)
                                       //    continue;
                                       // spdlog::debug("final frags size: {} dup frags size: {} {}",
                                       // _final_frags.size(), chr_id, _dup_frags[chr_id].size());
                                       // _final_frags.insert(_final_frags.end(), _dup_frags[chr_id].begin(),
                                       // _dup_frags[chr_id].end());
                                       for (auto& frag : _dup_frags[chr_id])
                                       {
                                           // spdlog::debug("frag: {} {} {} {}", frag.start, frag.end, frag.chr_id,
                                           // frag.barcode_id);
                                           _final_frags.push_back(std::move(frag));
                                       }
                                       _dup_frags[chr_id].clear();
                                       _dup_frags[chr_id].shrink_to_fit();
                                   }
                                   spdlog::debug("Final frags size: {}", _final_frags.size());

                                   fs::path out_frag_file = output_path / (run_name + FRAGMENT_FILE);
                                   spdlog::debug("Dump frags to: {}", out_frag_file.string());
                                   BGZF* out_frag;
                                   out_frag = bgzf_open(out_frag_file.c_str(), "w");
                                   for (auto& l : _final_frags)
                                   {
                                       // Standardized format output, insert one column data
                                       string s = _contig_names[l.chr_id] + "\t" + to_string(l.start) + "\t"
                                                  + to_string(l.end) + "\t" + _idx2drop[l.barcode_id] + "\t1\n";
                                       // spdlog::debug(s);
                                       [[maybe_unused]] auto ret = bgzf_write(out_frag, s.c_str(), s.size());
                                   }
                                   bgzf_close(out_frag);
                                   spdlog::info("Merge frags time(s): {:.2f}", t.toc(1000));

                                   if (tbx_index_build(out_frag_file.c_str(), 0, &tbx_conf_bed))
                                       spdlog::warn("Failed build frags index");
                                   else
                                       spdlog::info("Build frags index time(s): {:.2f}", t.toc(1000));
                               })
                           .name("Merge fragments");
    merge_frags.succeed(end_reanno);

    // Step 9: simple qc
    auto simple_qc = taskflow
                         .emplace(
                             [&]()
                             {
                                 Timer t;
                                 D2C::simpleQC(used_chrs);
                                 spdlog::info("Simple qc time(s): {:.2f}", t.toc(1000));
                             })
                         .name("Simple qc");
    simple_qc.succeed(end_reanno);

    // Step 10: final qc
    auto final_qc = taskflow
                        .emplace(
                            [&]()
                            {
                                Timer t;
                                D2C::finalQC();
                                spdlog::info("Final qc time(s): {:.2f}", t.toc(1000));
                            })
                        .name("Final qc");
    final_qc.succeed(barcode_merge);
    final_qc.succeed(simple_qc);
    final_qc.succeed(merge_frags);

    // Step 11: plot
    auto plot = taskflow
                    .emplace(
                        [&]()
                        {
                            Timer t;
                            D2C::plot();
                            spdlog::info("Plot time(s): {:.2f}", t.toc(1000));
                        })
                    .name("Plot");
    plot.succeed(barcode_merge);
    plot.succeed(determine_hq_beads);
    plot.succeed(sequence_saturation);

    executor.run(taskflow).wait();
    // taskflow.dump(std::cout);

    return 0;
}

bool operator<(const UniqBarcode& lhs, const UniqBarcode& rhs)
{
    if (lhs.start < rhs.start)
        return true;
    else if (lhs.start > rhs.start)
        return false;
    else
    {
        if (lhs.end < rhs.end)
            return true;
        else if (lhs.end > rhs.end)
            return false;
        else
        {
            if (lhs.barcode < rhs.barcode)
                return true;
            else
                return false;
        }
    }
}

int D2C::splitBamByChr(int chr_id)
{
    string chr_str = _contig_names[chr_id];
    spdlog::debug("Call splitBamByChr: {}", chr_str);

    auto& bedpes = _bedpes_by_chr[chr_id];
    if (!_is_bed)
    {
        std::unique_ptr< SamReader > sr = SamReader::FromFile(input_bam);
        if (!sr->QueryByContig(chr_id))
            return 0;

        map< string, pair< BamRecord, int > >           pe_dict;
        map< string, pair< BamRecord, int > >::iterator it;
        BamRecord                                       bam_record = bam_init1();
        int                                             pos        = -1;
        while (sr->next(bam_record))
        {
            ++pos;
            if ((bam_record->core.flag & FLAG) != FLAG)
                continue;

            BamRecord b = bam_init1();
            bam_copy1(b, bam_record);
            string qname = getQname(b);
            if (bam_has_suffix)
            {
                qname = qname.substr(0, qname.size() - 2);
            }
            it = pe_dict.find(qname);
            if (it == pe_dict.end())
            {
                pe_dict[qname] = { b, pos };
            }
            else
            {
                extractBedPE(it->second.first, b, bedpes, it->second.second, pos);
                bam_destroy1(it->second.first);
                pe_dict.erase(it);
                bam_destroy1(b);
            }
        }
        bam_destroy1(bam_record);
        for (auto& p : pe_dict)
            bam_destroy1(p.second.first);
        pe_dict.clear();
    }

    spdlog::debug("chr: {} frags size: {}", chr_str, bedpes.size());
    spdlog::debug("chr: {} memory(MB): {}", chr_str, physical_memory_used_by_process());

    // Load blacklist file and construct interval tree
    vector< Node > nodes;
    string         line;
    if (fs::exists(blacklist_file))
    {
        ifstream blf(blacklist_file, std::ifstream::in);
        while (std::getline(blf, line))
        {
            vector< string > vec_str = split_str(line, '\t');
            if (vec_str.size() != 3)
                continue;
            string chr = vec_str[0];
            if (chr != chr_str)
                continue;
            int  start = stoi(vec_str[1]);
            int  end   = stoi(vec_str[2]);
            Node node;
            node.lower = start;
            node.upper = end;
            node.value = chr;
            nodes.push_back(std::move(node));
        }
        blf.close();
    }
    MyTree mytree;
    for (auto& node : nodes)
    {
        mytree.insert(node);
    }
    // set< UniqBarcode > pcr_dup;
    spp::sparse_hash_set< UniqBarcode > pcr_dup;
    // unordered_set<string> pcr_dup;
    // Quantify the number of unique fragments per barcode
    unordered_map< int, int > bead_quant;
    // vector<string> bead_order;
    for (auto& b : bedpes)
    {
        // Filter for fragments overlapping the blacklist
        int        start = b.start;
        int        end   = b.end;
        MyInterval query_range{ start, end };

        const auto& res = mytree.query(query_range);
        // Discard the fragments overlapping the blacklist
        if (res.begin() != res.end())
            continue;

        int barcode = b.barcode;
        // string key     = to_string(start) + '\t' + to_string(end) + '\t' + to_string(barcode);
        UniqBarcode ub;
        ub.start   = start;
        ub.end     = end;
        ub.barcode = barcode;
        if (pcr_dup.count(ub) == 0)
        {
            pcr_dup.insert(ub);
            ++bead_quant[barcode];
        }
    }
    spdlog::debug("prc dup size: {}", pcr_dup.size());
    pcr_dup.clear();

    // unordered_set<string>().swap(pcr_dup);

    // Devel
    // if (chr_str == "chrMT")
    // {
    //     ofstream ofs(output_path/"chrMT.bead_counts.tsv", std::ofstream::out);
    //     for (auto& p : bead_quant)
    //         ofs<<p.first<<'\t'<<p.second<<endl;
    //     ofs.close();
    // }

    // Devel
    // int total_cnt = 0;
    // for (auto& b : bead_quant)
    //     total_cnt += b.second;
    // spdlog::debug("chr: {} bead_quant size: {} total count: {}", _contig_names[chr_id], bead_quant.size(),
    // total_cnt);

    // Merge bead quant of all chrs
    if (!mito_chrs.count(_contig_names[chr_id]))
    {
        std::lock_guard< std::mutex > guard(_merge_chr_mutex);
        for (auto& b : bead_quant)
            _total_bead_quant[b.first] += b.second;
        spdlog::debug("_total_bead_quant size: {}", _total_bead_quant.size());
    }

    return 0;
}

map< string, string > D2C::parseChrsFromBedFile(string filename)
{
    map< string, string > chrs;
    if (!fs::exists(filename))
        return chrs;

    ifstream ifs(filename, std::ifstream::in);
    string   line;
    while (std::getline(ifs, line))
    {
        vector< string > vec_s = split_str(line, '\t');
        if (vec_s.size() < 2)
            continue;
        chrs[vec_s[0]] = vec_s[1];
    }
    ifs.close();
    return chrs;
}

pair< double, double > D2C::parseBeadThreshold(string filename)
{
    pair< double, double > res;
    ifstream               ifs(filename, std::ifstream::in);
    string                 line;
    std::getline(ifs, line);
    if (!line.empty())
        res.first = stod(line);
    std::getline(ifs, line);
    if (!line.empty())
        res.second = stod(line);
    ifs.close();
    return res;
}

int D2C::determineHQBeads()
{
    // Dump total bead quant to csv for call R script
    FILE*    out_bead_quant;
    fs::path filename = output_path;
    filename /= (run_name + BARCODE_QUANT_FILE);
    spdlog::debug("Dump bead quant to:{}", filename.string());
    out_bead_quant = fopen(filename.c_str(), "w");
    for (auto& b : _total_bead_quant)
    {
        string s = int2Barcode(b.first) + _runname + FSEP + to_string(b.second) + "\n";
        fwrite(s.c_str(), 1, s.size(), out_bead_quant);
    }
    fclose(out_bead_quant);

    // Check beads size in case program segfault
    if (_total_bead_quant.empty())
    {
        return -1;
    }

    // Calculate bead threshold
    fs::path paras_file = output_path / (run_name + PARAM_FILE);

    if (barcode_threshold > 0)
    {
        // Use the top N parameter first
        vector< double > cnts;
        for (auto& b : _total_bead_quant)
        {
            cnts.push_back(b.second);
        }
        std::sort(cnts.begin(), cnts.end(), std::greater< double >());
        min_barcode_frags =
            barcode_threshold <= static_cast< int >(cnts.size()) ? cnts[barcode_threshold - 1] : cnts.back();

        ofstream ofs(paras_file.string(), std::ofstream::out);
        ofs << "bead_cutoff" << FSEP << min_barcode_frags << endl;
        ofs.close();
    }
    else if (min_barcode_frags != 0.0)
    {
        // Use the inflection parameter second
        ofstream ofs(paras_file.string(), std::ofstream::out);
        ofs << "bead_cutoff" << FSEP << min_barcode_frags << endl;
        ofs.close();
    }
    else
    {
        // Calculate the inflection point as default
        vector< double > cnts;
        for (auto& b : _total_bead_quant)
        {
            cnts.push_back(b.second);
        }

        // Dump cutoff through Rscript
        min_barcode_frags = barcodeRank(filename.string(), rank, paras_file.string());

        if (beads_force != 0)
        {
            std::sort(cnts.begin(), cnts.end(), std::greater< double >());
            int tmp = beads_force <= static_cast< int >(cnts.size()) ? cnts[beads_force - 1] : cnts.back();
            spdlog::debug("bead_force_frags: {} min_barcode_frags: {}", tmp, min_barcode_frags);
            if (min_barcode_frags < tmp)
            {
                min_barcode_frags = tmp;

                ofstream ofs(paras_file.string(), std::ofstream::out);
                // ofs.precision(4);
                ofs << "bead_cutoff" << FSEP << min_barcode_frags << endl;
                ofs.close();
            }
        }
    }

    // Do the filter
    for (auto& p : _total_bead_quant)
    {
        if (p.second >= min_barcode_frags)
            _hq_beads.insert(p.first);
    }

    spdlog::debug("bead cutoff: {}", min_barcode_frags);
    spdlog::debug("total beads num: {} filter by min frags cutoff: {}", _total_bead_quant.size(), _hq_beads.size());
    return 0;
}

int D2C::computeStatByChr(int chr_id)
{
    spdlog::debug("in computeStatByChr: {}", _contig_names[chr_id]);
    spdlog::debug("chr: {} memory(MB): {}", _contig_names[chr_id], physical_memory_used_by_process());
    // Filter fragments
    unordered_map< unsigned long long, int > dict;
    auto&                                    frags_data = _bedpes_by_chr[chr_id];
    vector< bool >                           frags_pos(frags_data.size());
    if (frags_data.empty())
        return 0;

    for (size_t i = 0; i < frags_data.size(); ++i)
    {
        auto& bedpe = frags_data[i];
        // Filter for eligible barcodes
        if (_hq_beads.count(bedpe.barcode) == 0)
            continue;

        int start = bedpe.start;
        int end   = bedpe.end;
        ++dict[(( unsigned long long )start << 32) + end];

        frags_pos[i] = true;
    }
    spdlog::debug("dict size:{} pos size: {}", dict.size(), frags_pos.size());
    spdlog::debug("chr: {} memory(MB): {}", _contig_names[chr_id], physical_memory_used_by_process());

    // Quantify NC + export
    auto& cnts = _total_nc_cnts[chr_id];
    for (auto& p : dict)
    {
        cnts[p.second] += p.second;
    }
    spdlog::debug("chr: {} nc size: {}", _contig_names[chr_id], cnts.size());

    // Pull out barcode for retained fragments
    for (size_t i = 0; i < frags_pos.size(); ++i)
    {
        if (!frags_pos[i])
            continue;

        auto& bedpe = frags_data[i];
        int   start = bedpe.start;
        int   end   = bedpe.end;
        if (dict[(( unsigned long long )start << 32) + end] > nc_threshold)
            frags_pos[i] = false;
    }

    spp::sparse_hash_set< UniqBarcode >           uniq_frags;
    spp::sparse_hash_map< int, vector< int > >    overlap_start, overlap_end;
    spp::sparse_hash_map< size_t, vector< int > > overlap_both;
    // size_t                                        count = 0;
    for (size_t i = 0; i < frags_pos.size(); ++i)
    {
        if (!frags_pos[i])
            continue;

        auto& bedpe   = frags_data[i];
        int   start   = bedpe.start;
        int   end     = bedpe.end;
        int   barcode = bedpe.barcode;

        UniqBarcode ub;
        ub.start   = start;
        ub.end     = end;
        ub.barcode = barcode;
        // string key     = to_string(start) + '\t' + to_string(end) + '\t' + to_string(barcode);
        if (uniq_frags.count(ub) != 0)
            continue;

        uniq_frags.insert(ub);

#ifdef UNIQ_FRAG_BOTH
        overlap_both[(( size_t )start << 32) + end].push_back(barcode);
#else
        overlap_start[start].push_back(barcode);
        overlap_end[end].push_back(barcode);
        // count += 2;
#endif
    }
    // spdlog::debug("uniq_frags.size: {} count: {}", uniq_frags.size(), count);
    // spdlog::debug("chr: {} before clear uniq_frags memory(MB): {}", _contig_names[chr_id],
    // physical_memory_used_by_process());
    uniq_frags.clear();
    // spdlog::debug("uniq_frags.size: {}", uniq_frags.size());
    // spdlog::debug("overlap_start size:{} overlap_end size:{}", overlap_start.size(), overlap_end.size());
    // spdlog::debug("chr: {} memory(MB): {}", _contig_names[chr_id], physical_memory_used_by_process());

    // Double to consider left and right inserts
    unordered_map< size_t, int > bead_cnts;
#ifdef UNIQ_FRAG_BOTH
    for (auto&& overlap : { overlap_both })
#else
    for (auto&& overlap : { overlap_start, overlap_end })
#endif
    {
        for (auto&& p : overlap)
        {
            auto&& v = p.second;
            if (v.size() == 1)
                continue;
            for (size_t i = 0; i < v.size(); ++i)
            {
                for (size_t j = i + 1; j < v.size(); ++j)
                {
                    if (v[i] == v[j])
                        continue;
                    if (int2Barcode(v[i]) > int2Barcode(v[j]))
                        ++bead_cnts[(( size_t )v[i] << 32) + v[j]];
                    else
                        ++bead_cnts[(( size_t )v[j] << 32) + v[i]];
                }
            }
        }
    }
    spdlog::debug("bead_cnts size:{}", bead_cnts.size());
    bead_cnts.swap(_total_bead_cnts[chr_id]);

    overlap_start.clear();
    overlap_end.clear();
    overlap_both.clear();
    // Devel
    // int line_num = 0, total_num = 0;
    // for (auto& p : bead_cnts)
    // {
    //     int count = p.second;
    //     if (count >= regularize_threshold)
    //     {
    //         ++line_num;
    //         total_num += count;
    //     }
    // }
    spdlog::debug("chr: {} memory(MB): {}", _contig_names[chr_id], physical_memory_used_by_process());

    return 0;
}

// change default n from 6 to user input
inline string substrRight(string s, int n)
{
    return s.substr(s.size() - n);
}

// Example: ATCG,TCGA
bool D2C::checkTn5(string s)
{
    // if (!tn5)
    //     return true;

    size_t pos = s.find(',');
    if (pos == std::string::npos)
        return false;

    string s1 = s.substr(0, pos);
    string s2 = s.substr(pos + 1);
    return (substrRight(s1, tn5) == substrRight(s2, tn5));
}
// Example: int32int32
bool D2C::checkTn5(size_t l)
{
    // if (!tn5)
    //     return true;

    string s1 = int2Barcode(int(l >> 32));
    string s2 = int2Barcode(int(l & 0xFFFFFFFF));
    // std::cout<<s1<<" "<<s2<<" "<<substrRight(s1, tn5)<<" "<<substrRight(s2, tn5)<<endl;
    return (substrRight(s1, tn5) == substrRight(s2, tn5));
}

int D2C::determineBarcodeMerge()
{
    spdlog::debug("tn5: {} regularize_threshold: {}", tn5, regularize_threshold);
    // Devel
    // fs::path temp_out = output_path / "temp.txt";
    // ofstream temp_ofs(temp_out.string(), std::ofstream::out);
    // for (auto& m : _total_bead_cnts)
    // {
    //     for (auto& p : m)
    //     {
    //         temp_ofs << p.second<< "\n";
    //     }
    // }
    // temp_ofs.close();
    spdlog::debug("_total_bead_cnts size: {}", _total_bead_cnts.size());
    // Merge all barcode count
    spp::sparse_hash_map< size_t, int > sum_dt;
    if (tn5 != 0)
    {
        for (auto& m : _total_bead_cnts)
        {
            for (auto& p : m)
            {
                // Only consider merging when Tn5 is the same
                // cout<<p.first<<" "<<p.second<<endl;
                if (p.second >= regularize_threshold && checkTn5(p.first))
                    sum_dt[p.first] += p.second;
            }
            // spdlog::debug("m size: {} sum_dt size: {}", m.size(), sum_dt.size());
            m.clear();
        }
    }
    else
    {
        for (auto& m : _total_bead_cnts)
        {
            for (auto& p : m)
            {
                // Only consider merging when Tn5 is the same
                if (p.second >= regularize_threshold)
                    sum_dt[p.first] += p.second;
            }
            // spdlog::debug("m size: {} sum_dt size: {}", m.size(), sum_dt.size());
            m.clear();
        }
    }
    _total_bead_cnts.clear();
    spdlog::debug("sum_dt size: {}", sum_dt.size());

    // Check beads size in case program segfault
    if (_total_bead_quant.empty())
    {
        return -1;
    }

    // Filter and calculate nBC
    vector< pair< int, int > >       nBC;
    spp::sparse_hash_map< int, int > count_dict;
    for (auto& p : _total_bead_quant)
    {
        if (_hq_beads.count(p.first) == 0)
            continue;
        nBC.push_back({ p.first, p.second });
        count_dict[p.first] = p.second * 2;
    }

    // Release memory of hq beads
    _hq_beads.clear();
    spdlog::debug("count_dict size: {}", count_dict.size());

    // Sort by count
    std::stable_sort(nBC.begin(), nBC.end(),
                     [](const pair< int, int >& a, const pair< int, int >& b) { return a.second > b.second; });

    // Calculate jaccard frag
    vector< pair< size_t, float > > ovdf;
    for (auto& p : sum_dt)
    {
        int   N_both       = p.second;
        int   N_barc1      = count_dict[p.first >> 32];
        int   N_barc2      = count_dict[p.first & 0xFFFFFFFF];
        float jaccard_frag = round((N_both) / (N_barc1 + N_barc2 - N_both + 0.05), 5);
        if (jaccard_frag <= 0.0)
            continue;
        ovdf.push_back({ p.first, jaccard_frag });
    }

    // Sort by jaccard_frag
    std::sort(ovdf.begin(), ovdf.end(),
              [](const pair< size_t, float >& a, const pair< size_t, float >& b) { return a.second > b.second; });

    fs::path paras_file = output_path / (run_name + PARAM_FILE);
    if (jaccard_threshold > 0)
    {
        // Use the top N parameter first
        vector< double > cnts;
        int              size = min(1000000, int(ovdf.size()));
        for (int i = 0; i < size; ++i)
        {
            cnts.push_back(ovdf[i].second);
        }

        std::sort(cnts.begin(), cnts.end(), std::greater< double >());
        min_jaccard_index =
            jaccard_threshold <= static_cast< int >(cnts.size()) ? cnts[jaccard_threshold - 1] : cnts.back();

        ofstream ofs(paras_file.string(), std::ofstream::out | std::ofstream::app);
        ofs << "cor_cutoff" << FSEP << min_jaccard_index << endl;
        ofs.close();
    }
    else if (min_jaccard_index != 0.0)
    {
        // Use the inflection parameter second
        ofstream ofs(paras_file.string(), std::ofstream::out | std::ofstream::app);
        ofs << "cor_cutoff" << FSEP << min_jaccard_index << endl;
        ofs.close();
    }
    else
    {
        // Calculate the inflection point as default
        fs::path temp_path = output_path / (run_name + TEMP_JACCARD_FILE);

        ofstream ofs(temp_path.string(), std::ofstream::out);
        ofs << "jaccard_distance" << endl;
        int size = min(1000000, int(ovdf.size()));
        for (int i = 0; i < size; ++i)
        {
            ofs << ovdf[i].second << endl;
        }
        ofs.close();

        min_jaccard_index = barcodeRank(temp_path.string(), rank, paras_file.string());
        fs::remove(temp_path);
    }

    spdlog::debug("cor cutoff: {}", min_jaccard_index);

    spdlog::debug("determineBarcodeMerge memory(MB): {}", physical_memory_used_by_process());

    // Export the implicated barcodes
    cmpFile  tbl_out;
    fs::path implicated_barcode_file = output_path / (run_name + IMPLICATED_BARCODES_FILE);
    tbl_out                          = cmpOpen(implicated_barcode_file.c_str());
    string header =
        "barcode1\tbarcode2\tFragnum_overlap\tFragnum_barcode1\tFragnum_barcode2\tjaccard_distance\tmerged\n";
    cmpFunc(tbl_out, header.c_str());
    for (size_t i = 0; i < ovdf.size(); ++i)
    {
        auto&  p  = ovdf[i];
        int    b1 = p.first >> 32;
        int    b2 = p.first & 0xFFFFFFFF;
        string s  = int2Barcode(b1) + _runname + FSEP + int2Barcode(b2) + _runname + FSEP;
        s += to_string(sum_dt[p.first]) + FSEP;

        int N_barc1 = count_dict[b1];
        int N_barc2 = count_dict[b2];
        s += to_string(N_barc1) + FSEP + to_string(N_barc2) + FSEP;
        s += f2str(p.second, 5) + FSEP;
        s += p.second >= min_jaccard_index ? "TRUE" : "FALSE";
        s += "\n";

        cmpFunc(tbl_out, s.c_str());
    }
    cmpClose(tbl_out);
    spdlog::debug("dump file: {}", implicated_barcode_file.string());
    sum_dt.clear();
    count_dict.clear();

    // Filter based on the min_jaccard_index
    // and prepare dict data
    map< int, vector< int > > bc1, bc2;
    for (size_t i = 0; i < ovdf.size(); ++i)
    {
        auto& p = ovdf[i];
        if (p.second <= min_jaccard_index)
            continue;

        int b1 = p.first >> 32;
        int b2 = p.first & 0xFFFFFFFF;
        bc1[b1].push_back(b2);
        bc2[b2].push_back(b1);
    }

    // Guess at how wide we need to make the barcodes to handle leading zeros
    int guess = ceil(log10(nBC.size()));
    // Map from barcode to position in nBC
    map< int, int > bar2pos;
    for (size_t i = 0; i < nBC.size(); ++i)
    {
        bar2pos[nBC[i].first] = i;
    }

    spdlog::debug("before barcode combine memory(MB): {}", physical_memory_used_by_process());
    // Devel
    // for (size_t i = 0; i < nBC.size(); ++i)
    // {
    //     string barcode = nBC[i].first;
    //     spdlog::debug("{} {}", i+1, barcode);
    // }
    // Loop through and eat up barcodes
    int idx = 1;
    for (size_t i = 0; i < nBC.size(); ++i)
    {

        int barcode = nBC[i].first;
        if (barcode == -1)
            continue;

        vector< int > barcode_combine;
        barcode_combine.push_back(barcode);

        // Find friends that are similarly implicated and append from Barcode 1
        if (bc1.count(barcode) != 0)
        {
            barcode_combine.insert(barcode_combine.end(), bc1[barcode].begin(), bc1[barcode].end());
        }
        // Find friends that are similarly implicated and append from Barcode 2
        if (bc2.count(barcode) != 0)
        {
            barcode_combine.insert(barcode_combine.end(), bc2[barcode].begin(), bc2[barcode].end());
        }

        // If user species one to one, then only remove that one barcode
        if (one_to_one)
            barcode_combine = { barcode };

        // Remove barcode which already mapped to drop_barcode
        int size = 0;
        for (unsigned int j = 0; j < barcode_combine.size(); ++j)
        {
            auto b = barcode_combine[j];
            if (_drop_barcodes.count(b) != 0)
            {
                // It will be discarded because not exists in map bar2pos
                barcode_combine[j] = -1;
            }
            else
                ++size;
        }

        // Make a drop barcode and save our progress
        stringstream ss;
        ss << "_BC" << std::setfill('0') << std::setw(guess) << idx << "_N" << std::setw(2) << size;
        for (auto& b : barcode_combine)
        {
            if (bar2pos.count(b) != 0)
            {
                string drop_barcode;
                if (tn5 == 0)
                    drop_barcode = run_name + ss.str();
                else
                    drop_barcode = run_name + "_Tn5-" + substrRight(int2Barcode(b), tn5) + ss.str();

                if (_idx2drop.empty() || drop_barcode != _idx2drop.back())
                    _idx2drop.push_back(std::move(drop_barcode));
                _drop_barcodes[b]     = _idx2drop.size() - 1;
                nBC[bar2pos[b]].first = -1;
            }
        }
        ++idx;
    }

    spdlog::debug("after barcode combine memory(MB): {}", physical_memory_used_by_process());

    FILE*    bt;
    fs::path barcode_translate_file = output_path / (run_name + BARCODE_TRANSLATE_FILE);
    bt                              = fopen(barcode_translate_file.c_str(), "w");
    for (auto& p : _drop_barcodes)
    {
        string s = int2Barcode(p.first) + _runname + "\t" + _idx2drop[p.second] + "\n";
        fwrite(s.c_str(), 1, s.size(), bt);
    }
    fclose(bt);

    // Finally, collate the NC values
    map< int, int > nc_data;
    for (auto& d : _total_nc_cnts)
    {
        for (auto& p : d)
            nc_data[p.first] += p.second;
        d.clear();
    }
    _total_nc_cnts.clear();

    // FILE*    nc_file_out;
    // fs::path nc_sum_file = output_path / (run_name + NC_STATS_FILE);
    // nc_file_out          = fopen(nc_sum_file.c_str(), "w");
    // header               = "NC_value\tNumberOfFragments\n";
    // fwrite(header.c_str(), 1, header.size(), nc_file_out);
    // for (auto& p : nc_data)
    // {
    //     string s = to_string(p.first) + "\t" + to_string(p.second) + "\n";
    //     fwrite(s.c_str(), 1, s.size(), nc_file_out);
    // }
    // fclose(nc_file_out);
    nc_data.clear();

    return 0;
}

inline int start(const string& s)
{
    size_t p1 = s.find("\t");
    size_t p2 = s.find("\t", p1 + 1);
    if (p2 == std::string::npos)
        return -1;
    return stoi(s.substr(p1 + 1, p2 - p1));
}

struct cmp
{
    bool operator()(const AnnotateFragment& a, const AnnotateFragment& b) const
    {
        return a.start < b.start;
    }
};

int D2C::reannotateFragByChr(int chr_id)
{
    spp::sparse_hash_set< AnnotateFragment > pcr_dup;
    unordered_set< int >                     qname_dup;
    // Store n_total and n_unique as pair
    map< string, pair< unsigned long, unsigned long > >& merge_ss = _frag_stats[chr_id];
    string                                               chr      = _contig_names[chr_id];
    // for mix species
    unsigned long increment = 0;
    if (!species_mix)
        increment = 1;
    else
    {
        bool human = chr.at(0) == 'h';
        bool mouse = chr.at(0) == 'm';
        if (human)
            increment = 1;
        else if (mouse)
            increment = 0x100000000;
    }

    auto&                                                                  frags_data = _bedpes_by_chr[chr_id];
    unordered_map< string, spp::sparse_hash_map< AnnotateFragment, int > > dups_per_cell;  // only used for saturation
    unordered_map< int, int >::iterator                                    it;
    for (auto& bedpe : frags_data)
    {
        // Filter for eligible barcodes
        int barcode = bedpe.barcode;
        it          = _drop_barcodes.find(barcode);
        if (it == _drop_barcodes.end())
            continue;

        string cell_barcode = _idx2drop[it->second];
        // string s            = chr + "\t" + to_string(bedpe.start) + "\t" + to_string(bedpe.end) + "\t" +
        // cell_barcode;
        AnnotateFragment frag(chr_id, bedpe.start, bedpe.end, it->second);
        auto&            p = merge_ss[cell_barcode];
        if (pcr_dup.count(frag) == 0)
        {
            // spdlog::debug("frag1: {} {} {} {}", frag.start, frag.end, frag.chr_id, frag.barcode_id);
            pcr_dup.insert(frag);
            // qname_dup.insert(bedpe.qname);
            qname_dup.insert(bedpe.qname1);
            qname_dup.insert(bedpe.qname2);
            p.second += increment;
        }
        p.first += increment;

        if (saturation_on)
        {
            dups_per_cell[cell_barcode][frag]++;
        }
    }
    qname_dup.swap(_keep_qnames[chr_id]);

    // Sort by start
    vector< AnnotateFragment > dup_keys(pcr_dup.begin(), pcr_dup.end());
    // spdlog::debug("dup keys size: {} pcr_dup size: {}", dup_keys.size(), pcr_dup.size());
    std::sort(dup_keys.begin(), dup_keys.end(), cmp());
    // for (auto& frag : dup_keys)
    // {
    //     spdlog::debug("frag1: {} {} {} {}", frag.start, frag.end, frag.chr_id, frag.barcode_id);
    // }
    dup_keys.swap(_dup_frags[chr_id]);

    // Release the memory
    pcr_dup.clear();
    _bedpes_by_chr[chr_id].clear();
    _bedpes_by_chr[chr_id].shrink_to_fit();

    spdlog::debug("reannotateFragByChr memory(MB): {}", physical_memory_used_by_process());

    std::lock_guard< std::mutex > guard(_merge_chr_mutex);
    // Merge frags data for sequencing saturation
    if (saturation_on)
    {
        saturation.addData(dups_per_cell);
    }

    return 0;
}

int D2C::annotateBamByChr(int chr_id)
{
    auto& keep_reads = _keep_qnames[chr_id];

    std::unique_ptr< SamReader > sr = SamReader::FromFile(input_bam);
    if (!sr->QueryByContig(chr_id))
        return 0;
    spdlog::debug("Call annotateBamByChr: {}", _contig_names[chr_id]);

    fs::path             output_bam_file = temp_bam_path / (_contig_names[chr_id] + ".bam");
    htsFile*             out             = hts_open(output_bam_file.c_str(), "wb");
    auto                 header          = sr->getHeader();
    [[maybe_unused]] int ret             = sam_hdr_write(out, header);

    BamRecord                           b = bam_init1();
    string                              bead_bc, drop_bc;
    unordered_map< int, int >::iterator it;
    // Iterate through bam
    int line = -1;
    while (sr->next(b))
    {
        ++line;
        if (keep_reads.count(line) == 0)
            continue;
        if (!getTag(b, barcode_tag.c_str(), bead_bc))
            continue;

        string b1 = bead_bc.substr(0, BLEN);
        string b2 = bead_bc.substr(BLEN, BLEN);
        if (_barcode2int.count(b1) == 0 || _barcode2int.count(b2) == 0)
        {
            spdlog::debug("Invalid barcode {} not exists in barcode list", bead_bc);
            continue;
        }
        // int tail = 0;
        // if (bead_bc.size() > BBLEN)
        // {
        //     string tmp = bead_bc.substr(BBLEN);
        //     if (_tail2int.count(tmp) == 0)
        //     {
        //         spdlog::warn("Invalid tn5 {} not exists in tn5.list", barcode);
        //         continue;
        //     }
        //     tail = _tail2int[tmp];
        // }
        // int barcode = (tail << BBIT * 2) + (_barcode2int[b1] << BBIT) + _barcode2int[b2];
        int barcode = (_barcode2int[b1] << BBIT) + _barcode2int[b2];
        it          = _drop_barcodes.find(barcode);
        if (it == _drop_barcodes.end())
            continue;
        drop_bc = _idx2drop[it->second];

        // Handle droplet barcodes that we want to consider writing out
        bam_aux_append(b, drop_tag.c_str(), 'Z', drop_bc.size() + 1, ( uint8_t* )drop_bc.c_str());
        [[maybe_unused]] int ret = sam_write1(out, header, b);
    }

    bam_destroy1(b);
    hts_close(out);

    return 0;
}

int D2C::simpleQC(vector< int >& used_chrs)
{
    map< string, pair< unsigned long, unsigned long > > nuclear;
    int                                                 single_mito_pos = -1, human_mc_pos = -1, mouse_mc_pos = -1;
    for (auto& chr_id : used_chrs)
    {
        if (_contig_names[chr_id] == single_mc)
        {
            single_mito_pos = chr_id;
            continue;
        }
        else if (_contig_names[chr_id] == human_mc)
        {
            human_mc_pos = chr_id;
            continue;
        }
        else if (_contig_names[chr_id] == mouse_mc)
        {
            mouse_mc_pos = chr_id;
            continue;
        }

        for (auto& p : _frag_stats[chr_id])
        {
            if (p.second.first == 0 || p.second.second == 0)
                continue;
            nuclear[p.first].first += p.second.first;
            nuclear[p.first].second += p.second.second;
        }
        _frag_stats[chr_id].clear();
    }
    if (species_mix)
    {
        // for mixed species, there is more values to process
        if (human_mc_pos == -1 && mouse_mc_pos == -1)
        {
            for (auto& [cell_barcode, count_pair] : nuclear)
            {
                SumStat ss;
                ss.drop_barcode     = cell_barcode;
                ss.mito_total       = 0;
                ss.mito_uniq        = 0;
                ss.human_mito_total = 0;
                ss.human_mito_uniq  = 0;
                ss.mouse_mito_total = 0;
                ss.mouse_mito_uniq  = 0;

                int totalCnt, humanCnt, mouseCnt;
                decodeSpeciesMix(count_pair.first, totalCnt, humanCnt, mouseCnt);
                ss.nuclear_total = totalCnt;
                ss.human_total   = humanCnt;
                ss.mouse_total   = mouseCnt;
                decodeSpeciesMix(count_pair.second, totalCnt, humanCnt, mouseCnt);
                ss.nuclear_uniq = totalCnt;
                ss.human_uniq   = humanCnt;
                ss.mouse_uniq   = mouseCnt;

                _sum_stats.push_back(ss);
            }
        }
        else
        {
            map< string, pair< unsigned long, unsigned long > > human_stats, mouse_stats;
            if (human_mc_pos != -1)
                human_stats.swap(_frag_stats[human_mc_pos]);
            if (mouse_mc_pos != -1)
                mouse_stats.swap(_frag_stats[mouse_mc_pos]);

            map< string, pair< unsigned long, unsigned long > >::iterator it;
            for (auto& cb : _idx2drop)
            {
                it = nuclear.find(cb);
                if (it == nuclear.end())
                    continue;

                SumStat ss;
                ss.drop_barcode = cb;
                int totalCnt, humanCnt, mouseCnt;
                decodeSpeciesMix(it->second.first, totalCnt, humanCnt, mouseCnt);
                ss.nuclear_total = totalCnt;
                ss.human_total   = humanCnt;
                ss.mouse_total   = mouseCnt;
                decodeSpeciesMix(it->second.second, totalCnt, humanCnt, mouseCnt);
                ss.nuclear_uniq = totalCnt;
                ss.human_uniq   = humanCnt;
                ss.mouse_uniq   = mouseCnt;

                it = human_stats.find(cb);
                if (it != human_stats.end())
                {
                    decodeSpeciesMix(it->second.first, totalCnt, humanCnt, mouseCnt);
                    ss.human_mito_total = totalCnt;
                    decodeSpeciesMix(it->second.second, totalCnt, humanCnt, mouseCnt);
                    ss.human_mito_uniq = totalCnt;
                }

                it = mouse_stats.find(cb);
                if (it != mouse_stats.end())
                {
                    decodeSpeciesMix(it->second.first, totalCnt, humanCnt, mouseCnt);
                    ss.mouse_mito_total = totalCnt;
                    decodeSpeciesMix(it->second.second, totalCnt, humanCnt, mouseCnt);
                    ss.mouse_mito_uniq = totalCnt;
                }

                ss.mito_total = ss.human_mito_total + ss.mouse_mito_total;
                ss.mito_uniq  = ss.human_mito_uniq + ss.mouse_mito_uniq;

                _sum_stats.push_back(ss);
            }
        }
    }
    else
    {
        // No mito found
        // Stay useful infomation except mito
        if (single_mito_pos == -1)
        {
            for (auto& [cell_barcode, count_pair] : nuclear)
            {
                SumStat ss;
                ss.drop_barcode  = cell_barcode;
                ss.mito_total    = 0;
                ss.mito_uniq     = 0;
                ss.nuclear_total = count_pair.first;
                ss.nuclear_uniq  = count_pair.second;

                _sum_stats.push_back(ss);
            }
        }
        else
        {
            auto&                                                         mito = _frag_stats[single_mito_pos];
            map< string, pair< unsigned long, unsigned long > >::iterator it;
            for (auto& cb : _idx2drop)
            {
                it = nuclear.find(cb);
                if (it == nuclear.end())
                    continue;

                SumStat ss;
                ss.drop_barcode  = cb;
                ss.nuclear_total = it->second.first;
                ss.nuclear_uniq  = it->second.second;
                if (mito.count(cb))
                {
                    ss.mito_total = mito[cb].first;
                    ss.mito_uniq  = mito[cb].second;
                }
                else
                {
                    ss.mito_total = 0;
                    ss.mito_uniq  = 0;
                }
                _sum_stats.push_back(ss);
            }
        }
    }
    return 0;
}

int D2C::finalQC()
{
    // Load tss file for finding overlaps
    map< string, MyTree > mytrees;
    vector< Node >        nodes;
    if (fs::exists(trans_file))
    {
        ifstream tss_ifs(trans_file, std::ifstream::in);
        string   line;
        while (std::getline(tss_ifs, line))
        {
            vector< string > vec_s = split_str(line, '\t');
            if (vec_s.size() < 3)
                continue;

            Node node;
            node.lower = stoi(vec_s[1]) - 1000;
            node.upper = stoi(vec_s[2]) + 1000;
            node.value = vec_s[0];
            nodes.push_back(std::move(node));
        }

        for (auto& node : nodes)
        {
            mytrees[node.value].insert(node);
        }
    }

    // Store has_overlap and insert size as pair
    map< string, SummaryData > summary;
    for (auto& l : _final_frags)
    {
        string     chr   = _contig_names[l.chr_id];
        int        start = l.start;
        int        end   = l.end;
        MyInterval query_range{ start, end };
        auto&      sd = summary[_idx2drop[l.barcode_id]];
        if (mytrees.count(chr) != 0)
        {
            const auto& tree = mytrees.at(chr);
            const auto& res  = tree.query(query_range);
            if (res.begin() != res.end())
                ++sd.overlaps;
        }
        int insert_size = end - start;
        sd.insert_size.push_back(insert_size);
    }
    mytrees.clear();
    _final_frags.clear();
    _final_frags.shrink_to_fit();

    // Deal with FRIP if we have a peak file
    if (peak_file != "")
    {
        // TODO(fxzhao): implement this option
    }

    // Summarize frag attributes
    for (auto& p : summary)
    {
        auto& sd = p.second;
        // sd.mean_insert_size =
        //     std::accumulate(sd.insert_size.begin(), sd.insert_size.end(), 0.0) / sd.insert_size.size();
        // std::nth_element(sd.insert_size.begin(), sd.insert_size.begin() + int(0.5 * sd.insert_size.size()),
        //                  sd.insert_size.end());
        // sd.median_insert_size = *(sd.insert_size.begin() + int(0.5 * sd.insert_size.size()));
        // sd.frip               = 0;
        sd.tss_proportion = sd.overlaps * 1.0 / sd.insert_size.size();
    }

    for (auto& l : _sum_stats)
    {
        int nuclear_total = l.nuclear_total;
        int nuclear_uniq  = l.nuclear_uniq;
        l.dup_proportion  = get_dup_proportion(nuclear_total, nuclear_uniq);
    }

    // Add barcodes back if we need to (specified with the one to one option)
    if (one_to_one)
    {
        // TODO(fxzhao)
    }

    // Create a summarized experiment if we have a valid peak file
    if (peak_file != "")
    {
        // TODO(fxzhao)
    }

    // Export QC stats
    FILE*    qc_out;
    fs::path out_qc_file = output_path / (run_name + QC_STATS_FILE);
    qc_out               = fopen(out_qc_file.c_str(), "w");
    string header = "CellBarcode\ttotalFrags\tuniqueFrags\ttotalMitoFrags\tuniqueMitoFrags\tduplicateProportion\ttssPro"
                    "portion\tMitoProportion\n";
    if (species_mix)
        header = "CellBarcode\ttotalFrags\tuniqueFrags\ttotalMitoFrags\tuniqueMitoFrags\ttotalHumanMitoFrags\tuniqHuman"
                 "MitoFrags\ttotalMouseMitoFrags\tuniqMouseMitoFrags\ttotalHumanFrags\tuniqueHumanFrags\ttotalMouseFrag"
                 "s\tuniqueMouseFrags\tduplicateProportion\ttssProportion\tMitoProportion\n";

    fwrite(header.c_str(), 1, header.size(), qc_out);
    map< string, SummaryData >::iterator it;
    for (auto& l : _sum_stats)
    {
        it = summary.find(l.drop_barcode);
        if (it == summary.end())
            continue;

        auto&  sd = it->second;
        string s  = l.drop_barcode + FSEP + to_string(l.nuclear_total) + FSEP + to_string(l.nuclear_uniq) + FSEP
                   + to_string(l.mito_total) + FSEP + to_string(l.mito_uniq) + FSEP;

        if (species_mix)
            s += to_string(l.human_mito_total) + FSEP + to_string(l.human_mito_uniq) + FSEP
                 + to_string(l.mouse_mito_total) + FSEP + to_string(l.mouse_mito_uniq) + FSEP + to_string(l.human_total)
                 + FSEP + to_string(l.human_uniq) + FSEP + to_string(l.mouse_total) + FSEP + to_string(l.mouse_uniq)
                 + FSEP;

        s += f2str(l.dup_proportion, 4) + FSEP + f2str(sd.tss_proportion, 4) + FSEP
             + f2str((l.mito_uniq * 1.0) / l.nuclear_uniq, 4) + "\n";
        fwrite(s.c_str(), 1, s.size(), qc_out);
    }
    fclose(qc_out);

    return 0;
}

int D2C::plot()
{
    fs::path script_path = bin_path / PLOT_SCRIPT;
    string command = "python -W ignore " + script_path.string() + " " + output_path.string() + " " + run_name + " 2>&1";
    vector< string > cmd_result;
    int              cmd_rtn = exec_shell(command.c_str(), cmd_result);
    for (const auto& line : cmd_result)
        if (!line.empty())
            spdlog::info(line);
    if (cmd_rtn == 0)
        spdlog::info("Plot success");
    else
    {
        spdlog::error("Plot fail, rtn:{}", cmd_rtn);
        return -1;
    }

    return 0;
}

bool D2C::parseBarcodeList()
{
    ifstream ifs(barcode_list, std::ifstream::in);
    string   line;
    int      pos = 0;
    while (std::getline(ifs, line))
    {
        _barcode_names.push_back(line);
        _barcode2int[line] = pos++;
    }
    if (_barcode_names.empty())
        return false;
    return true;
}

bool D2C::parseTn5List()
{
    _tail_names.push_back("");
    ifstream ifs(tn5_list, std::ifstream::in);
    string   line;
    int      pos = 1;
    while (std::getline(ifs, line))
    {
        _tail_names.push_back(line);
        _tail2int[line] = pos++;
    }
    // if (_tail_names.empty())
    //     return false;
    return true;
}

inline string D2C::int2Barcode(int i)
{
    // spdlog::debug("i: {}", i);
    string res = _barcode_names[(i >> BBIT) & BMASK] + _barcode_names[i & BMASK] + _tail_names[(i >> BBIT * 2) & RMASK];
    return res;
}

inline string D2C::int2Tn5(int i)
{
    string res = _tail_names[(i >> BBIT * 2) & RMASK];
    return res;
}

bool D2C::parseRunnameList()
{
    // int pos = 1;
    // Deprecate
    // if (!barcode_runname_list.empty())
    // {
    //     ifstream ifs(barcode_runname_list, std::ifstream::in);
    //     string   line;
    //     while (std::getline(ifs, line))
    //     {
    //         if (line.substr(0, 1) != BSEP)
    //             line.insert(0, BSEP);
    //         _runnames.push_back(line);
    //         _runname2int[line] = pos++;
    //     }
    // }
    // else
    {
        // Get run names from bam file
        if (_is_bed)
        {
            // for bed, there is not exists runname in the end of barcode,
            // but there will be 10bp tn5

            // ifstream ifs(input_bam, std::ifstream::in);
            // string   line;
            // while (std::getline(ifs, line))
            // {
            //     vector< string > vec_s = split_str(line, '\t');
            //     if (vec_s.size() < 4)
            //         continue;

            //     string& barcode = vec_s[3];
            //     if (barcode.size() <= BBLEN)
            //         break;

            //     string name = barcode.substr(BBLEN);
            //     if (!name.empty() && _runname2int.count(name))
            //         continue;
            //     _runnames.push_back(name);
            //     _runname2int[name] = pos++;
            //     break;
            // }
            // ifs.close();
        }
        else
        {
            std::unique_ptr< SamReader >                          sr      = SamReader::FromFile(input_bam);
            std::vector< std::pair< std::string, unsigned int > > contigs = sr->getContigs();

            BamRecord bam_record = bam_init1();
            string    barcode("NA");
            for (size_t chr_id = 0; chr_id < contigs.size(); ++chr_id)
            {
                if (!sr->QueryByContig(chr_id))
                    continue;
                // Users insist on not mixing multiple batches of data, so we just read one record in bam
                if (sr->next(bam_record) && getTag(bam_record, barcode_tag.c_str(), barcode))
                {
                    // Barcode format: bc1bc2-runname, and the size of bc1 or bc2 is 10
                    // There is exists new format without runname, only bc1bc2
                    if (barcode.size() <= BBLEN)  // means no runname
                        _runname = "";
                    else
                        _runname = barcode.substr(BBLEN);
                    break;
                }
            }
            bam_destroy1(bam_record);
        }
    }
    spdlog::debug("Barcode runname: {}", _runname);
    return true;
}

double D2C::barcodeRank(string filename, string rank_type, string outfile)
{
    // Run Rscript to get knee or inflection
    fs::path script_path = bin_path / R_SCRIPT;
    string   command =
        "Rscript " + script_path.string() + " -i " + filename + " -t " + rank_type + " -o " + outfile + " 2>&1";
    vector< string > cmd_result;
    int              cmd_rtn = exec_shell(command.c_str(), cmd_result);
    for (const auto& line : cmd_result)
        if (!line.empty())
            spdlog::info(line);
    if (cmd_rtn == 0)
        spdlog::info("Rscript run success");
    else
    {
        spdlog::error("Rscript run fail, rtn:{}", cmd_rtn);
        throw std::runtime_error("Rscript failed. Please check if Rscript in your $PATH.");
    }

    // Parse result of Rscript
    double   res = 0;
    ifstream ifs(outfile, std::ifstream::in);
    string   line;
    while (std::getline(ifs, line))
    {
        vector< string > vec_s = split_str(line, '\t');
        if (vec_s.size() != 2)
            continue;
        res = std::stod(vec_s[1]);
    }

    return res;
}