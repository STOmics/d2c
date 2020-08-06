/*
 * File: bap.cpp
 * Created Data: 2020-7-23
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "bap.h"


#include "timer.h"
#include "utility.h"
#include "bamCat.h"
#include "gz_io.h"


#include <spdlog/spdlog.h>

#include <taskflow/taskflow.hpp>

#include <fstream>
#include <algorithm>
#include <exception>


const int FLAG = 2;
const int MAX_INSERT = 2000;

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
    int         upper;
    int         lower;
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
    int overlaps; // The number of has overlap with TSS
    float tss_proportion; // overlaps / total
    float mean_insert_size;
    int median_insert_size;
    float frip;   // The mean value of peak
    vector<int> insert_size;
};

float get_dup_proportion(int t, int u)
{
    return (t-u)*1.0 / t;
}

inline float aux(float x, int c, int n)
{
    return (c / x - 1 + exp(-n / x));
}
int get_library_size(int t, int u)
{
    float m = 1;
    float M = 100;
    int n_dup = t - u + 1;
    if (u > t || aux(m*u, u, t) < 0 || u < 0 || t < 0 || n_dup < 0)
        return 0;

    while (aux(M*u, u, t) > 0)
        M *= 10.0;
    
    for (int i = 0; i <= 40; ++i)
    {
        float mid = (m + M) / 2;
        float v = aux(mid * u, u, t);
        if ( v == 0 )
            break;
        else if (v > 0)
            m = mid;
        else
            M = mid;
    }
    return round(u*(m+M)/2.0);
}

void Bap::extractBedPE(const BamRecord b1, const BamRecord b2, vector<Bedpe>& bedpes, int l1, int l2)
{
    // Initialize BEDPE variables
    string chrom1, chrom2, strand1, strand2;
    int start1, start2, end1, end2;
    start1 = start2 = end1 = end2 = -1;
    chrom1 = chrom2 = strand1 = strand2 = ".";
    int min_map_quality = 0, max_map_quality = 0;

    // Extract relevant info for end 1
    if (isMapped(b1))
    {
        //chrom1 = refs[b1->core.tid];
        start1 = b1->core.pos;
        end1 = bam_endpos(b1);
        if (isReverse(b1))
            strand1 = "-";
        else
            strand1 = "+";
    }

    // Extract relevant info for end 2
    if (isMapped(b2))
    {
        //chrom2 = refs[b2->core.tid];
        start2 = b2->core.pos;
        end2 = bam_endpos(b2);
        if (isReverse(b2))
            strand2 = "-";
        else
            strand2 = "+";
    }

    // Swap the ends if necessary
    //if (chrom1 > chrom2 || ((chrom1 == chrom2) && (start1 > start2)))
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
        //ofs << chrom1 <<"\t"<<start1<<"\t"<<end2<<"\t"<<get_qame(b1)<<endl;
        string barcode("NA");
        if (getTag(b1, barcode_tag.c_str(), barcode))
        {
            Bedpe bedpe;
            bedpe.start = start1;
            bedpe.end = end2;
            // bedpe.qname = getQname(b1);
            bedpe.qname1 = l1;
            bedpe.qname2 = l2;
            
            string b1 = barcode.substr(0, 10);
            string b2 = barcode.substr(10, 10);
            if (_barcode2int.count(b1) == 0 || _barcode2int.count(b2) == 0)
            {
                spdlog::warn("Invalid barcode {} not exists in barcode list", barcode);
                return;
            }
            bedpe.barcode = (_barcode2int[b1] << 16) + _barcode2int[b2];
            //cout<<b1<<" "<<b2<<" "<<bedpe.barcode<<endl;
            bedpes.push_back(bedpe);
        }
    }
}

Bap::Bap(string input_bam, string output_path, string barcode_tag, int mapq, int cores, string run_name, bool tn5,
        double min_barcode_frags, double min_jaccard_index, string ref, string mito_chr, string bed_genome_file,
        string blacklist_file, string trans_file, bool species_mix, string bin_path, double barcode_threshold,
        double jaccard_threshold, bool saturation_on, string barcode_list) :
        input_bam(input_bam),
        output_path(output_path),
        barcode_tag(barcode_tag),
        mapq(mapq),
        cores(cores),
        run_name(run_name),
        tn5(tn5),
        min_barcode_frags(min_barcode_frags),
        min_jaccard_index(min_jaccard_index),
        ref(ref),
        mito_chr(mito_chr),
        bed_genome_file(bed_genome_file),
        blacklist_file(blacklist_file),
        trans_file(trans_file),
        species_mix(species_mix),
        bin_path(bin_path),
        barcode_threshold(barcode_threshold),
        jaccard_threshold(jaccard_threshold),
        saturation_on(saturation_on),
        barcode_list(barcode_list)
{
    nc_threshold = 6;
    regularize_threshold = 4;
    one_to_one = false;
    drop_tag = "DB";
}
    
int Bap::run()
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
    spdlog::info("Prepare time(s): {:.2f}", timer.toc(1000));

    // Execute taskflow
    taskflow();
    spdlog::info("Execute taskflow time(s): {:.2f}", timer.toc(1000));

    // Clean up
    //fs::remove_all(temp_bam_path);
    spdlog::info("Clean up time(s): {:.2f}", timer.toc(1000));

    return 0;
}

int Bap::taskflow()
{
    tf::Executor executor;
    tf::Taskflow taskflow;

    taskflow.name("Bap");

    // Step 1: split bam by chr
    auto samReader = SamReader::FromFile(input_bam);
    std::vector< std::pair< std::string, unsigned int > > contigs = samReader->getContigs();
    spdlog::debug("Bam contigs num: {}", contigs.size());
    for (auto& p : contigs)
        _contig_names.push_back(p.first);

    // Verify that the supplied reference genome and the bam have overlapping chromosomes
    auto bed_chrs = parseChrsFromBedFile();
    spdlog::debug("bed_chrs num: {}", bed_chrs.size());

    vector<int> used_chrs;
    for (int i = 0; i < contigs.size(); ++i)
    {
        if (bed_chrs.count(_contig_names[i]) != 0)
            used_chrs.push_back(i);
    }
    if (used_chrs.empty())
    {
        spdlog::error("Found no overlapping chromosomes between bam and reference. Check reference genome specification with the -r flag");
        return -1;
    }
    else
    {
        spdlog::info("Found {} chromosomes for analysis (including mitochondria)", used_chrs.size());
    }

    spdlog::debug("Initial memory(MB): {}", physical_memory_used_by_process());
    _bedpes_by_chr.resize(contigs.size());
    auto [S, T] = taskflow.parallel_for(
        used_chrs.begin(),
        used_chrs.end(),
        [&] (int chr_id) {
            Timer t;
            Bap::splitBamByChr(chr_id);
            spdlog::info("Split bam by chr: {} time(s): {:.2f}", contigs[chr_id].first, t.toc(1000));
        }
    );

    S.name("Start SplitBam/Assemble/Annotate");
    T.name("End SplitBam/Assemble/Annotate");
    S.for_each_successor([contigs, s=0] (tf::Task successor) mutable
    {
        successor.name((contigs.begin()+(s++))->first);
    });

    // Step 1.1: sequencing saturation
    
    auto sequence_saturation = taskflow.emplace([&] ()
    {
        if (saturation_on)
        {
            Timer t;
            fs::path sat_out_file = output_path / (run_name+".sequenceSaturation.tsv");
            saturation.calculateSaturation(sat_out_file.string());
            spdlog::info("Sequencing saturation time(s): {:.2f}", t.toc(1000));
        }
    }).name("Sequence saturation");
    sequence_saturation.succeed(T);

    // Step 2: determine hight quality beads
    auto determine_hq_beads = taskflow.emplace([&] ()
    {
        spdlog::debug("SplitBam memory(MB): {}", physical_memory_used_by_process());
        Timer t;
        Bap::determineHQBeads();
        spdlog::info("Determine high-quality beads time(s): {:.2f}", t.toc(1000));
    }).name("Determine HQ Beads");
    determine_hq_beads.succeed(T);

    // Step 3: calculate stat by chr
    _total_bead_cnts.resize(contigs.size());
    _total_nc_cnts.resize(contigs.size());
    // cout<<"contig size: "<<_bedpes_by_chr.size()<<endl;
    // for (int i = 0; i < _bedpes_by_chr.size(); ++i)
    //     if (_bedpes_by_chr[i].size() > 0)
    //         cout<<"has data chr: "<<i<<endl;
    auto [start_cal, end_cal] = taskflow.parallel_for(
        used_chrs.begin(),
        used_chrs.end(),
        [&] (int chr_id) {
            // Skip the mito chrom
            if (_contig_names[chr_id] == mito_chr) return;
            Timer t;
            Bap::computeStatByChr(chr_id);
            spdlog::info("Compute stat by chr: {} time(s): {:.2f}", _contig_names[chr_id], t.toc(1000));
        }
    );
    start_cal.name("Start Cal");
    end_cal.name("End Cal");
    start_cal.for_each_successor([contigs, s=0] (tf::Task successor) mutable
    {
        successor.name((contigs.begin()+(s++))->first);
    });
    start_cal.succeed(determine_hq_beads);
    start_cal.succeed(T);

    // Step 4: barcode merge
    auto barcode_merge = taskflow.emplace([&] ()
    {
        Timer t;
        Bap::determineBarcodeMerge();
        spdlog::info("Determine barcode merge time(s): {:.2f}", t.toc(1000));
    }).name("Determine Barcode Merge");
    barcode_merge.succeed(end_cal);
    barcode_merge.succeed(determine_hq_beads);

    // Step 5: reannotate frag data and get summary stats by chr
    _keep_qnames.resize(contigs.size());
    _dup_frags.resize(contigs.size());
    _frag_stats.resize(contigs.size());
    auto [start_reanno, end_reanno] = taskflow.parallel_for(
        used_chrs.begin(),
        used_chrs.end(),
        [&] (int chr_id) {
            Timer t;
            Bap::reannotateFragByChr(chr_id);
            spdlog::info("Reannotate frags by chr: {} time(s): {:.2f}", _contig_names[chr_id], t.toc(1000));
        }
    );
    start_reanno.name("Start Reannotate");
    end_reanno.name("End Reannotate");
    start_reanno.for_each_successor([contigs,s=0] (tf::Task successor) mutable
    {
        successor.name((contigs.begin()+(s++))->first);
    });
    start_reanno.succeed(barcode_merge);
    start_reanno.succeed(T);

    // Step 6: annotate bam file by chr
    auto [start_annobam, end_annobam] = taskflow.parallel_for(
        used_chrs.begin(),
        used_chrs.end(),
        [&] (int chr_id) {
            // Skip the mito chrom
            if (_contig_names[chr_id] == mito_chr) return;
            Timer t;
            Bap::annotateBamByChr(chr_id);
            spdlog::info("Annotate bam file by chr: {} time(s): {:.2f}", _contig_names[chr_id], t.toc(1000));
        }
    );
    start_annobam.name("Start Annotate bam");
    end_annobam.name("End Annotate bam");
    start_annobam.for_each_successor([contigs,s=0] (tf::Task successor) mutable
    {
        successor.name((contigs.begin()+(s++))->first);
    });
    start_annobam.succeed(end_reanno);
    start_annobam.succeed(T);

    // Step 7: merge bam files
    auto merge_bam = taskflow.emplace([&] ()
    {
        Timer t;
        spdlog::debug("Merge bam");
        std::vector< std::string > bam_files;
        for (auto& chr_id : used_chrs)
        {
            // Skip the mito chrom
            if (_contig_names[chr_id] == mito_chr) continue;
            fs::path tmp_bam_file = temp_bam_path / (contigs[chr_id].first+".bam");
            if (fs::exists(tmp_bam_file))
                bam_files.push_back(tmp_bam_file.string());
        }
        fs::path output_bam_file = output_path / (run_name + ".bam");
        spdlog::debug("be merged file size: {}", bam_files.size());
        // for (auto& name : bam_files)
        //     spdlog::debug("merge bam: {}", name);
        try 
        {
            int cmd_rtn = bam_cat(bam_files, nullptr, output_bam_file.c_str(), nullptr, 0);
            if (cmd_rtn == 0)
                spdlog::info("Merge bam file success");
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
    }).name("Merge bam");
    merge_bam.succeed(end_annobam);
    //end_annobam.precede(merge_bam);

    // Step 8: merge fragment files
    auto merge_frags = taskflow.emplace([&] ()
    {
        spdlog::debug("Merge frags");
        Timer t;
        for (auto& chr_id : used_chrs)
        {
            // Skip the mito chrom
            if (_contig_names[chr_id] == mito_chr) continue;
            _final_frags.insert(_final_frags.end(), _dup_frags[chr_id].begin(), _dup_frags[chr_id].end());
            _dup_frags[chr_id].clear();
            _dup_frags[chr_id].shrink_to_fit();
        }
        fs::path out_frag_file = output_path / (run_name+".fragments.tsv.gz");
        spdlog::debug("Dump frags to: {}", out_frag_file.string());
        cmpFile out_frag;
        out_frag = cmpOpen(out_frag_file.c_str());
        for (auto& l : _final_frags)
        {
            // Standardized format output
            string s = l + "\t1\n";
            cmpFunc(out_frag, s.c_str());
        }
        cmpClose(out_frag);
        spdlog::info("Merge frags time(s): {:.2f}", t.toc(1000));
    }).name("Merge fragments");
    merge_frags.succeed(end_reanno);

    // Step 9: simple qc
    auto simple_qc = taskflow.emplace([&] ()
    {
        Timer t;
        map<string, pair<int,int>> nuclear;
        int mito_pos = -1;
        for (auto& chr_id : used_chrs)
        {
            if (_contig_names[chr_id] == mito_chr)
            {
                mito_pos = chr_id;
                continue;
            }

            for (auto& p : _frag_stats[chr_id])
            {
                if (p.second.first == 0 || p.second.second == 0) continue;
                nuclear[p.first].first += p.second.first;
                nuclear[p.first].second += p.second.second;
            }
            _frag_stats[chr_id].clear();
        }
        assert(mito_pos != -1);
        auto& mito = _frag_stats[mito_pos];
        map<string, pair<int,int>>::iterator it;
        for (auto& p : mito)
        {
            if (p.second.first == 0 || p.second.second == 0) continue;
            it = nuclear.find(p.first);
            if (it == nuclear.end()) continue;

            SumStat ss;
            ss.drop_barcode = p.first;
            ss.nuclear_total = it->second.first;
            ss.nuclear_uniq = it->second.second;
            ss.mito_total = p.second.first;
            ss.mito_uniq = p.second.second;
            _sum_stats.push_back(ss);
        }

        fs::path out_ss_file = output_path / (run_name+".basicQC.tsv");
        FILE* out_ss;
        out_ss = fopen(out_ss_file.c_str(), "w");
        string header = "cell_barcode\ttotalNuclearFrags\tuniqueNuclearFrags\ttotalMitoFrags\tuniqueMitoFrags\n";
        fwrite(header.c_str(), 1, header.size(), out_ss);
        for (auto& l : _sum_stats)
        {
            string s = l.drop_barcode+'\t'+to_string(l.nuclear_total)+'\t'+to_string(l.nuclear_uniq)+'\t'+
                        to_string(l.mito_total)+'\t'+to_string(l.mito_uniq)+'\n';
            fwrite(s.c_str(), 1, s.size(), out_ss);
        }
        fclose(out_ss);

        spdlog::info("Simple qc time(s): {:.2f}", t.toc(1000));
    }).name("Simple qc");
    simple_qc.succeed(end_reanno);

    // Step 10: final qc
    auto final_qc = taskflow.emplace([&] ()
    {
        Timer t;
        Bap::finalQC();
        spdlog::info("Final qc time(s): {:.2f}", t.toc(1000));
    }).name("Final qc");
    final_qc.succeed(barcode_merge);
    final_qc.succeed(simple_qc);
    final_qc.succeed(merge_frags);

    // Step 11: plot
    auto plot = taskflow.emplace([&] ()
    {
        Timer t;
        Bap::plot();
        spdlog::info("Plot time(s): {:.2f}", t.toc(1000));
    }).name("Plot");
    plot.succeed(barcode_merge);
    plot.succeed(determine_hq_beads);

    executor.run(taskflow).wait();
    //taskflow.dump(std::cout);

    return 0;
}

int Bap::splitBamByChr(int chr_id)
{
    std::unique_ptr<SamReader> sr = SamReader::FromFile(input_bam);
    if (!sr->QueryByContig(chr_id)) return 0;
    spdlog::debug("Call splitBamByChr: {}", _contig_names[chr_id]);

    // Map {qname:barcode} used in later
    //map<string, string> qname2barcodes;

    //auto contigs = sr->getContigs();
    string chr_str = _contig_names[chr_id];
    // BamRecord bam_record1 = bam_init1();
    // BamRecord bam_record2 = bam_init1();
    // while (true)
    // {
    //     if (!sr->next(bam_record)) break;
    //     while (bam_record)
    //     if (bam_record->core.qual <= mapq) continue;
    //     string qname = bam_get_qname(bam_record);
    //     string barcode("NA");
    //     getTag(bam_record, barcode_tag.c_str(), barcode);
    //     qname2barcodes[qname] = barcode;
    // }
    map<string, pair<BamRecord, int>> pe_dict;
    map<string, pair<BamRecord, int>>::iterator it;
    auto& bedpes = _bedpes_by_chr[chr_id];
    BamRecord bam_record = bam_init1();
    int pos = -1;
    while (sr->next(bam_record))
    {
        ++pos;
        if ((bam_record->core.flag & FLAG) != FLAG) 
            continue;

        BamRecord b = bam_init1();
        bam_copy1(b, bam_record);
        string qname = getQname(b);
        it = pe_dict.find(qname);
        if (it == pe_dict.end())
        {
            pe_dict[qname] = {b, pos};
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
    // while (sr->next(bam_record1, FLAG))
    // {
    //     //spdlog::debug("qname:{}", getQname(bam_record1));
    //     sr->next(bam_record2, FLAG);
    //     if (getQname(bam_record1) != getQname(bam_record2))
    //     {
    //         while (getQname(bam_record1) != getQname(bam_record2))
    //         {
    //             bam_copy1(bam_record1, bam_record2);
    //             sr->next(bam_record2, FLAG);
    //         }
    //         extractBedPE(bam_record1, bam_record2, bedpes);
    //     }
    //     else if (isPaired(bam_record1) && isPaired(bam_record2))
    //     {
    //         extractBedPE(bam_record1, bam_record2, bedpes);
    //     }
    // }
    spdlog::debug("chr: {} frags size: {}", _contig_names[chr_id], bedpes.size());

    // Devel
    // if (chr_str == "chrMT")
    // {
    //     ofstream ofs(output_path/"chrMT.bedpe.annotated.tsv", std::ofstream::out);
    //     for (auto& p : bedpes)
    //         ofs<<p.start<<'\t'<<p.end<<'\t'<<p.qname<<'\t'<<p.barcode<<endl;
    //     ofs.close();
    // }

    // bam_destroy1(bam_record1);
    // bam_destroy1(bam_record2);

    // for (int i = 0; i < 10; ++i)
    //     cout<<bedpes[i].start<<" "<<bedpes[i].end<<" "<<bedpes[i].qname<<" "<<bedpes[i].barcode<<endl;
    // for (auto& b : bedpes)
    //     spdlog::debug("bedpe qname:{}", b.qname);

    // Load blacklist file and construct interval tree
    ifstream blf(blacklist_file, std::ifstream::in);
    vector<Node> nodes;
    string line;
    while (std::getline(blf, line))
    {
        vector<string> vec_str = split_str(line, '\t');
        if (vec_str.size() != 3) continue;
        string chr = vec_str[0];
        if (chr != chr_str) continue;
        int start = stoi(vec_str[1]);
        int end = stoi(vec_str[2]);
        Node node;
        node.lower = start;
        node.upper = end;
        node.value = chr;
        nodes.push_back(std::move(node));
    }
    blf.close();

    MyTree mytree;
    for (auto& node : nodes)
    {
        mytree.insert(node);
    }

    set<string> pcr_dup;
    // Quantify the number of unique fragments per barcode
    map<int, int> bead_quant;
    //vector<string> bead_order;
    for (auto& b : bedpes)
    {
        // Filter for fragments overlapping the blacklist
        int start = b.start;
        int end = b.end;
        MyInterval query_range {start, end};
        
        const auto& res = mytree.query(query_range);
        // Discard the fragments overlapping the blacklist
        if (res.begin() != res.end()) continue;

        int barcode = b.barcode;
        string key = to_string(start)+'\t'+to_string(end)+'\t'+to_string(barcode);
        if (pcr_dup.count(key) == 0)
        {
            pcr_dup.insert(key);
            // if (bead_quant.count(barcode) == 0)
            //     bead_order.push_back(barcode);
            ++bead_quant[barcode];
        }
    }
    pcr_dup.clear();

    // Devel
    // if (chr_str == "chrMT")
    // {
    //     ofstream ofs(output_path/"chrMT.bead_counts.tsv", std::ofstream::out);
    //     for (auto& p : bead_quant)
    //         ofs<<p.first<<'\t'<<p.second<<endl;
    //     ofs.close();
    // }
    
    // Devel
    int total_cnt = 0;
    for (auto& b : bead_quant)
        total_cnt += b.second;
    spdlog::debug("chr: {} bead_quant size: {} total count: {}", _contig_names[chr_id], bead_quant.size(), total_cnt);

    // Merge bead quant of all chrs
    std::lock_guard<std::mutex> guard(_merge_chr_mutex);
    for (auto& b : bead_quant)
        _total_bead_quant[b.first] += b.second;
    // Merge frags data for sequencing saturation
    // if (saturation_on)
    // {
    //     for (auto& bedpe : bedpes)
    //         saturation.addData(bedpe.start, bedpe.end, bedpe.barcode, chr_id);
    // }

    // for (auto& b : bead_order)
    // {
    //     if (_total_bead_quant.count(b) == 0)
    //         _total_bead_order.push_back(b);
    //     _total_bead_quant[b] += bead_quant[b];
    // }

    return 0;
}

map<string, string> Bap::parseChrsFromBedFile()
{
    map<string, string> chrs;
    ifstream ifs(bed_genome_file, std::ifstream::in);
    string line;
    while (std::getline(ifs, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() != 2) continue;
        chrs[vec_s[0]] = vec_s[1];
    }
    ifs.close();
    return chrs;
}

pair<double, double> Bap::parseBeadThreshold(string filename)
{
    pair<double, double> res;
    ifstream ifs(filename, std::ifstream::in);
    string line;
    std::getline(ifs, line);
    if (!line.empty()) res.first = stod(line);
    std::getline(ifs, line);
    if (!line.empty()) res.second = stod(line);
    ifs.close();
    return res;
}

int Bap::determineHQBeads()
{
    // Dump total bead quant to csv for call R script
    FILE * out_bead_quant;
    fs::path filename = output_path;
    filename /= (run_name+".barcodeQuantSimple.csv"); 
    spdlog::debug("Dump bead quant to:{}", filename.string());
    out_bead_quant = fopen(filename.c_str(), "w");
    for (auto& b : _total_bead_quant)
    //for (auto& b : _total_bead_order)
    {
        string s = int2Barcode(b.first) + "," + to_string(b.second) + "\n";
        //string s = b + "," + to_string(_total_bead_quant[b]) + "\n";
        fwrite(s.c_str(), 1, s.size(), out_bead_quant);
    }
    fclose(out_bead_quant);

    // TODO(fxzhao): change it to cpp
    // Call R script to calculate bead threshold
    fs::path paras_file = output_path / (run_name+".bapParam.csv");
    ofstream ofs(paras_file.string(), std::ofstream::out);
    ofs.precision(15);
    if (min_barcode_frags == 0.0)
    {
        fs::path script_path = bin_path / "10b_knee_execute.R";
        string command = "Rscript "+script_path.string()+" "+filename.string()+" 1 V2";
        vector<string> cmd_result;
        int cmd_rtn = exec_shell(command.c_str(), cmd_result);
        for (const auto& line : cmd_result)
            if (!line.empty())
                spdlog::error(line);
        if (cmd_rtn == 0)
            spdlog::info("Execute Rscript success");
        else
        {
            spdlog::error("Execute Rscript fail, rtn:{}", cmd_rtn);
            return -1;
        }

        // Define the set of high-quality bead barcodes
        auto paras = parseBeadThreshold(filename.string()+"_kneeValue.txt");
        min_barcode_frags = paras.first;
        double call_threshold = paras.second;
        
        ofs << "bead_threshold_nosafety,"<<call_threshold<<endl;
    }

    // Add threshold for reducing data to be processed
    vector<int> bead_nums;
    for (auto& p : _total_bead_quant)
        bead_nums.push_back(p.second);
    std::nth_element(bead_nums.begin(), bead_nums.begin()+int(barcode_threshold*bead_nums.size()), bead_nums.end(), std::greater<int>());
    int calculated_barcode_frags = *(bead_nums.begin()+int(barcode_threshold*bead_nums.size()));
    spdlog::debug("barcode filter rate: {} calculated frags: {} min barcode frags: {}", barcode_threshold, calculated_barcode_frags, min_barcode_frags);
    if (calculated_barcode_frags > min_barcode_frags)
    {
        spdlog::debug("Set min barcode frags from {} to {}", min_barcode_frags, calculated_barcode_frags);
        min_barcode_frags = calculated_barcode_frags;
    }

    ofs << "bead_threshold,"<<min_barcode_frags<<endl;
    ofs.close();
				
    // Do the filter   
    for (auto& p : _total_bead_quant)
    {
        if (p.second >= min_barcode_frags)
            _hq_beads.insert(p.first);
    }
    // Export high-quality beads
    fs::path out_hq_file = output_path / (run_name+".HQbeads.tsv");
    ofstream out_hq(out_hq_file.string(), std::ofstream::out);
    for (auto& b : _hq_beads)
        out_hq << int2Barcode(b) << "\n";
    out_hq.close();

    spdlog::debug("bead threshold: {}", min_barcode_frags);
    spdlog::debug("total beads num: {} filter by min frags threshold: {}", _total_bead_quant.size(), _hq_beads.size());
    return 0;
}

int Bap::computeStatByChr(int chr_id)
{
    spdlog::debug("in computeStatByChr: {}", _contig_names[chr_id]);
    // Filter fragments
    map<unsigned long long, int> dict;
    vector<int> frags_pos;
    auto& frags_data = _bedpes_by_chr[chr_id];
    if (frags_data.empty()) return 0;

    for (size_t i = 0; i < frags_data.size(); ++i)
    {
        auto& bedpe = frags_data[i];
        // Filter for eligible barcodes
        if (_hq_beads.count(bedpe.barcode) == 0) continue;

        int start = bedpe.start;
        int end = bedpe.end;
        ++dict[((unsigned long long)start << 32) + end];
        
        frags_pos.push_back(i);
    }

    // Quantify NC + export
    auto& cnts = _total_nc_cnts[chr_id];
    int total = 0;
    for (auto& p : dict)
    {
        cnts[p.second] += p.second;
        total += p.second;
    }
    
    spdlog::debug("chr: {} nc size: {} total: {}", _contig_names[chr_id], cnts.size(), total);
    // ofstream out_nc_count(out_nc_count_file, std::ofstream::out);
    // for (auto& p : cnts)
    // {
    //     out_nc_count<<p.first<<"\t"<<p.second<<endl;
    // }
    // out_nc_count.close();

    // Pull out barcode for retained fragments
    for (size_t i = 0; i < frags_pos.size(); ++i)
    {
        auto& bedpe = frags_data[frags_pos[i]];
        int start = bedpe.start;
        int end = bedpe.end;
        if (dict[((unsigned long long)start << 32) + end] > nc_threshold)
            frags_pos[i] = -1;
    }

    set<string> uniq_frags;
    map<int, vector<int>> overlap_start, overlap_end;
    for (auto& pos : frags_pos)
    {
        if (pos == -1) continue;

        auto& bedpe = frags_data[pos];
        int start = bedpe.start;
        int end = bedpe.end;
        int barcode = bedpe.barcode;
        string key = to_string(start)+'\t'+to_string(end)+'\t'+to_string(barcode);
        if (uniq_frags.count(key) != 0) continue;

        uniq_frags.insert(key);
        overlap_start[start].push_back(barcode);
        overlap_end[end].push_back(barcode);
    }

    // Double to consider left and right inserts
    map<size_t, int> bead_cnts;
    for (auto& overlap : {overlap_start, overlap_end})
    {
        for (auto& p : overlap)
        {
            auto& v = p.second;
            if (v.size() == 1) continue;
            for (size_t i = 0; i < v.size(); ++i)
            {
                for (size_t j = i+1; j < v.size(); ++j)
                {
                    if (v[i] == v[j]) continue;
                    if (v[i] > v[j])
                        ++bead_cnts[((size_t) v[i] << 32) + v[j]];
                    else
                        ++bead_cnts[((size_t) v[j] << 32) + v[i]];
                }
            }
        }
    }
    bead_cnts.swap(_total_bead_cnts[chr_id]);
    
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
    
   
    return 0;
}

inline string substrRight(string s, int n = 6)
{
    return s.substr(s.size()-n);
}

// Example: ATCG,TCGA
bool Bap::checkTn5(string s)
{
    if (!tn5) return true;

    size_t pos = s.find(',');
    if (pos == std::string::npos) return false;

    string s1 = s.substr(0, pos);
    string s2 = s.substr(pos+1);
    return (substrRight(s1) == substrRight(s2));
}
// Example: int32int32
bool Bap::checkTn5(size_t l)
{
    if (!tn5) return true;

    string s1 = _barcode_names[l >> 48] + _barcode_names[(l >> 32) & 0xFFFF];
    string s2 = _barcode_names[(l >> 16) & 0xFFFF] + _barcode_names[l & 0xFFFF];
    return (substrRight(s1) == substrRight(s2));
}

int Bap::determineBarcodeMerge()
{
    spdlog::debug("tn5: {}", tn5);
    // Merge all barcode count
    map<size_t, int> sum_dt;
    for (auto& m : _total_bead_cnts)
    {
        for (auto& p : m)
        {
            // Only consider merging when Tn5 is the same
            if (p.second >= regularize_threshold && checkTn5(p.first))
                sum_dt[p.first] += p.second;
        }
        m.clear();
    }
    _total_bead_cnts.clear();
     
    // Filter and calculate nBC
    vector<pair<int, int>> nBC;
    map<int, int> count_dict;
    for (auto& p : _total_bead_quant)
    {
        if (_hq_beads.count(p.first) == 0) continue;
        nBC.push_back({p.first, p.second});
        count_dict[p.first] = p.second * 2;
    }

    // Devel
    // ifstream tmp("/hwfssz5/ST_BIGDATA/USER/zhaofuxiang/test_bap2/data2/bapfile/knee/Parotid_gland_20200520_AJ34.barcodeQuantSimple.csv");
    // string line;
    // while (std::getline(tmp, line))
    // {
    //     vector<string> vec_s = split_str(line, ',');
        
    //     if (_hq_beads.count(vec_s[0]) == 0) continue;
    //     nBC.push_back({vec_s[0], stoi(vec_s[1])});
    // }

    // Release memory of hq beads
    _hq_beads.clear();

    // Sort by count
    std::stable_sort(nBC.begin(), nBC.end(), [](const pair<int, int>& a, const pair<int, int>& b) {
        return a.second > b.second;
    });


   
    // Calculate jaccard frag
    vector<pair<size_t, float>> ovdf;
    for (auto& p : sum_dt)
    {
        int N_both = p.second;
        int N_barc1 = count_dict[p.first >> 32];
        int N_barc2 = count_dict[p.first & 0xFFFFFFFF];
        float jaccard_frag = round((N_both)/(N_barc1+N_barc2-N_both+0.05), 5);
        if (jaccard_frag <= 0.0) continue;
        ovdf.push_back({p.first, jaccard_frag});
    }
    
    // Sort by jaccard_frag
    std::sort(ovdf.begin(), ovdf.end(), [](const pair<size_t, float>& a, const pair<size_t, float>& b) {
        return a.second > b.second;
    });

    fs::path paras_file = output_path / (run_name+".bapParam.csv");
    ofstream ofs(paras_file.string(), std::ofstream::out | std::ofstream::app);
    ofs.precision(15);
    // Call knee if we need to
    if (min_jaccard_index == 0.0)
    {
        // Prepare jaccard frag data for calling R script
        fs::path filename = output_path / (run_name+".jaccard.csv");
        ofstream jaccard_out(filename.string(), std::ofstream::out);
        int size = min(1000000, int(ovdf.size()));
        for (int i = 0; i < size; ++i)
            jaccard_out<<ovdf[i].second<<"\n";
        jaccard_out.close();

        fs::path script_path = bin_path / "11b_knee_execute.R";
        string command = "Rscript "+script_path.string()+" "+filename.string()+" 1 V1";
        vector<string> cmd_result;
        int cmd_rtn = exec_shell(command.c_str(), cmd_result);
        for (const auto& line : cmd_result)
            if (!line.empty())
                spdlog::error(line);
        if (cmd_rtn == 0)
            spdlog::info("Execute Rscript success");
        else
        {
            spdlog::info("Execute Rscript fail, rtn:{}", cmd_rtn);
            return -1;
        }

        // Define the set of high-quality bead barcodes
        auto paras = parseBeadThreshold(filename.string()+"_kneeValue.txt");
        min_jaccard_index = paras.first;
        double call_threshold = paras.second;

        ofs << "jaccard_threshold_nosafety,"<<call_threshold<<endl;
    }

    float calculated_jaccard_index = ovdf[int(ovdf.size()*jaccard_threshold)].second;
    spdlog::debug("jaccard filter rate: {} calculated jaccard: {} min barcode frags: {}", jaccard_threshold, calculated_jaccard_index, min_jaccard_index);
    if (calculated_jaccard_index > min_jaccard_index)
    {
        spdlog::debug("Set min jaccard index from {} to {}", min_jaccard_index, calculated_jaccard_index);
        min_jaccard_index = calculated_jaccard_index;
    }

    ofs << "jaccard_threshold,"<<min_jaccard_index<<endl;
    ofs.close();
    spdlog::debug("jaccard threshold: {}", min_jaccard_index);

    // Export the implicated barcodes
    cmpFile tbl_out;
    fs::path implicated_barcode_file = output_path / (run_name+".implicatedBarcodes.csv.gz"); 
    tbl_out = cmpOpen(implicated_barcode_file.c_str());
    string header = "barc1,barc2,N_both,N_barc1,N_barc2,jaccard_frag,merged\n";
    cmpFunc(tbl_out, header.c_str());
    for (size_t i = 0; i < ovdf.size(); ++i)
    {
        auto& p = ovdf[i];
        int b1 = p.first >> 32;
        int b2 = p.first & 0xFFFFFFFF;
        string s = int2Barcode(b1)+","+int2Barcode(b2) + ",";
        s += to_string(sum_dt[p.first]) + ",";
        
        int N_barc1 = count_dict[b1];
        int N_barc2 = count_dict[b2];
        s += to_string(N_barc1) + "," + to_string(N_barc2) + ",";
        s += f2str(p.second, 5) + ",";
        s += p.second > min_jaccard_index ? "TRUE" : "FALSE";
        s += "\n";

        cmpFunc(tbl_out, s.c_str());
    }
    cmpClose(tbl_out);

    // Filter based on the min_jaccard_index 
    // and prepare dict data
    map<int, vector<int>> bc1, bc2;
    for (size_t i = 0; i < ovdf.size(); ++i)
    {
        auto& p = ovdf[i];
        if (p.second <= min_jaccard_index) continue;
        
        int b1 = p.first >> 32;
        int b2 = p.first & 0xFFFFFFFF;
        bc1[b1].push_back(b2);
        bc2[b2].push_back(b1);
    }

    // Guess at how wide we need to make the barcodes to handle leading zeros
    int guess = ceil(log10(nBC.size()));
    // Map from barcode to position in nBC
    map<int, int> bar2pos;
    for (size_t i = 0; i < nBC.size(); ++i)
    {
        bar2pos[nBC[i].first] = i;
    }

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
        if (barcode == -1) continue;
        
        vector<int> barcode_combine;
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
        if (one_to_one) barcode_combine = {barcode};

        // Make a drop barcode and save our progress
        stringstream ss;
        ss << "_BC"<<std::setfill ('0') << std::setw (guess) << idx<<
                "_N"<<std::setw (2)<<barcode_combine.size();
        for (auto& b : barcode_combine)
        {
            if (bar2pos.count(b) != 0)
            {
                string drop_barcode;
                if (!tn5)
                    drop_barcode = run_name + ss.str();
                else
                    drop_barcode = run_name + "_Tn5-" + substrRight(_barcode_names[b]) + ss.str();
                _drop_barcodes[b] = drop_barcode;
                nBC[bar2pos[b]].first = -1;
            }
        }
        ++idx;
    }
    
    FILE * bt;
    fs::path barcode_translate_file = output_path / (run_name+".barcodeTranslate.tsv");
    bt = fopen(barcode_translate_file.c_str(), "w");
    for (auto& p : _drop_barcodes)
    {
        string s = int2Barcode(p.first) + "\t" + p.second+"\n";
        fwrite(s.c_str(), 1, s.size(), bt);
    }
    fclose(bt);

    // Finally, collate the NC values
    map<int, int> nc_data;
    for (auto& d : _total_nc_cnts)
    {
        for (auto& p : d)
            nc_data[p.first] += p.second;  
        d.clear();
    }
    _total_nc_cnts.clear();

    FILE * nc_file_out;
    fs::path nc_sum_file = output_path / (run_name+".NCsumstats.tsv");
    nc_file_out = fopen(nc_sum_file.c_str(), "w");
    header = "NC_value\tNumberOfFragments\n";
    fwrite(header.c_str(), 1, header.size(), nc_file_out);
    for (auto& p : nc_data)
    {
        string s = to_string(p.first) + "\t" + to_string(p.second) + "\n";
        fwrite(s.c_str(), 1, s.size(), nc_file_out);
    }
    fclose(nc_file_out);
    nc_data.clear();

    return 0;
}

inline int start(const string& s)
{
    size_t p1 = s.find("\t");
    size_t p2 = s.find("\t", p1+1);
    if (p2 == std::string::npos) return -1;
    return stoi(s.substr(p1+1, p2-p1));
}

struct cmp 
{
    bool operator() (const string& a, const string& b) const
    {
        return start(a) < start(b);
    }
};

int Bap::reannotateFragByChr(int chr_id)
{
    set<string> pcr_dup;
    set<int> qname_dup;
    // Store n_total and n_unique as pair
    map<string, pair<int, int>>& merge_ss = _frag_stats[chr_id];
    string chr = _contig_names[chr_id];
    auto& frags_data = _bedpes_by_chr[chr_id];
    map<int, string>::iterator it;
    for (auto& bedpe : frags_data)
    {
        // Filter for eligible barcodes
        int barcode = bedpe.barcode;
        it = _drop_barcodes.find(barcode);
        if (it == _drop_barcodes.end()) continue;

        string cell_barcode = it->second;
        string s = chr+"\t"+to_string(bedpe.start)+"\t"+to_string(bedpe.end)+"\t"+cell_barcode;
        if (pcr_dup.count(s) == 0)
        {
            pcr_dup.insert(s);
            // qname_dup.insert(bedpe.qname);
            qname_dup.insert(bedpe.qname1);
            qname_dup.insert(bedpe.qname2);
            ++merge_ss[cell_barcode].second;
        }

        ++merge_ss[cell_barcode].first;
    }
    qname_dup.swap(_keep_qnames[chr_id]);

    // Sort by start
    vector<string> dup_keys(pcr_dup.begin(), pcr_dup.end());
    std::sort(dup_keys.begin(), dup_keys.end(), cmp());
    dup_keys.swap(_dup_frags[chr_id]);

    pcr_dup.clear();
    // Release the memory
    
    _bedpes_by_chr[chr_id].clear();
    _bedpes_by_chr[chr_id].shrink_to_fit();

    //cout<<"chr "<<chr_id<<" keep_qname size: "<<_keep_qnames[chr_id].size()<<" dup_frags size: "<<_dup_frags[chr_id].size()<<endl;

    // Export duplicated fragments
    // FILE *frag_anno;
    // fs::path path = output_path;
    // fs::path frag_anno_file = path / "frag.bedpe..tsv";
    // frag_anno = fopen(frag_anno_file.c_str(), "w");
    // for (auto& k : dup_keys)
    // {
    //     string s = k + "\t" + pcr_dup[k] + "\n";
    //     fwrite(s.c_str(), 1, s.size(), frag_anno);
    // }
    // fclose(frag_anno);

    // Export summary statistics
    // FILE *frag_ss;
    // fs::path frag_ss_file = path / "frag.sumstats.tsv";
    // frag_ss = fopen(frag_ss_file.c_str(), "w");
    // for (auto& p : merge_ss)
    // {
    //     auto& pp = p.second;
    //     if (pp.first == 0 || pp.second == 0) continue;
    //     string s = p.first + "\t" + to_string(pp.first) + "\t" + to_string(pp.second)
    //         + "\t" + chr + "\n";
    //     fwrite(s.c_str(), 1, s.size(), frag_ss);
    // }
    // fclose(frag_ss);


    return 0;
}

int Bap::annotateBamByChr(int chr_id)
{
    auto& keep_reads = _keep_qnames[chr_id];

    std::unique_ptr<SamReader> sr = SamReader::FromFile(input_bam);
    if (!sr->QueryByContig(chr_id)) return 0;
    spdlog::debug("Call annotateBamByChr: {}", _contig_names[chr_id]);

    fs::path output_bam_file = temp_bam_path / (_contig_names[chr_id] + ".bam");
    htsFile* out = hts_open(output_bam_file.c_str(), "wb");
    auto header = sr->getHeader();
    [[maybe_unused]] int ret = sam_hdr_write(out, header);

    BamRecord b = bam_init1();
    string bead_bc, drop_bc;
    map<int, string>::iterator it;
    // Iterate through bam
    int line = -1;
    while (sr->next(b))
    {
        ++line;
        if (!getTag(b, barcode_tag.c_str(), bead_bc)) continue;
        
        string b1 = bead_bc.substr(0,10);
        string b2 = bead_bc.substr(10, 10);
        if (_barcode2int.count(b1) == 0 || _barcode2int.count(b2) == 0)
        {
            spdlog::debug("Invalid barcode {} not exists in barcode list", bead_bc);
            continue;
        }
        int barcode = (_barcode2int[b1] << 16) + _barcode2int[b2];
        it = _drop_barcodes.find(barcode);
        if (it == _drop_barcodes.end()) continue;
        drop_bc = it->second;

        // Handle droplet barcodes that we want to consider writing out
        // string qname = bam_get_qname(b);
        // if (keep_reads.count(qname) == 0) continue;
        if (keep_reads.count(line) == 0) continue;
        bam_aux_append(b, drop_tag.c_str(), 'Z', drop_bc.size()+1, ( uint8_t* )drop_bc.c_str());
        [[maybe_unused]] int ret = sam_write1(out, header, b);
    }
   
    bam_destroy1(b);
    hts_close(out);

    return 0;
}

int Bap::finalQC()
{
    // Load tss file for finding overlaps
    ifstream tss_ifs(trans_file, std::ifstream::in);
    string line;
    vector<Node> nodes;
    while (std::getline(tss_ifs, line))
    {
        vector<string> vec_s = split_str(line, '\t');
        if (vec_s.size() < 3) continue;

        Node node;
        node.lower = stoi(vec_s[1]) - 1000;
        node.upper = stoi(vec_s[2]) + 1000;
        node.value = vec_s[0];
        nodes.push_back(std::move(node));
    }

    map<string, MyTree> mytrees;
    for (auto& node : nodes)
    {
        mytrees[node.value].insert(node);
    }

    // Store has_overlap and insert size as pair
    map<string, SummaryData> summary;
    
    for (auto& l : _final_frags)
    {
        vector<string> vec_s = split_str(l, '\t');
        if (vec_s.size() < 4) continue;

        string chr = vec_s[0];
        int start = stoi(vec_s[1]);
        int end = stoi(vec_s[2]);
        MyInterval query_range {start, end};
        auto& sd = summary[vec_s[3]];
        if (mytrees.count(chr) != 0)
        {
            const auto& tree = mytrees.at(chr);
            const auto& res = tree.query(query_range);
            
            if (res.begin() != res.end())
                ++sd.overlaps;
        }
        int insert_size = end - start;
        sd.insert_size.push_back(insert_size);
    }
    mytrees.clear();

    // Deal with FRIP if we have a peak file
    if (peak_file != "")
    {
        // TODO(fxzhao): implement this option
    }

    // Summarize frag attributes
    for (auto& p : summary)
    {
        auto& sd = p.second;
        sd.mean_insert_size = std::accumulate(sd.insert_size.begin(), sd.insert_size.end(), 0.0) / 
            sd.insert_size.size();
        std::nth_element(sd.insert_size.begin(), 
            sd.insert_size.begin()+int(0.5*sd.insert_size.size()),
            sd.insert_size.end());
        sd.median_insert_size = *(sd.insert_size.begin()+int(0.5*sd.insert_size.size()));
        sd.frip = 0;
        sd.tss_proportion = sd.overlaps * 1.0 / sd.insert_size.size();
    }

    for (auto& l : _sum_stats)
    {
        int nuclear_total = l.nuclear_total;
        int nuclear_uniq = l.nuclear_uniq;
        l.dup_proportion = get_dup_proportion(nuclear_total, nuclear_uniq);
        l.library_size = get_library_size(nuclear_total, nuclear_uniq);
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
    FILE * qc_out;
    fs::path out_qc_file = output_path / (run_name+".QCstats.csv");
    qc_out = fopen(out_qc_file.c_str(), "w");
    string header = "DropBarcode,totalNuclearFrags,uniqueNuclearFrags,totalMitoFrags,uniqueMitoFrags,duplicateProportion,librarySize,meanInsertSize,medianInsertSize,tssProportion,FRIP\n";
    fwrite(header.c_str(), 1, header.size(), qc_out);
    map<string, SummaryData>::iterator it;
    for (auto& l : _sum_stats)
    {
        it = summary.find(l.drop_barcode); 
        if (it == summary.end()) continue;

        string s = l.drop_barcode+","+to_string(l.nuclear_total)+","+to_string(l.nuclear_uniq)+","+
                    to_string(l.mito_total)+","+to_string(l.mito_uniq)+","+f2str(l.dup_proportion, 3)+","+
                    to_string(l.library_size)+",";
        auto& sd = it->second;
        s += f2str(sd.mean_insert_size, 1)+","+to_string(sd.median_insert_size)+","+
            f2str(sd.tss_proportion, 4)+","+f2str(sd.frip, 4)+"\n";
        fwrite(s.c_str(), 1, s.size(), qc_out);
    }
    fclose(qc_out);

    return 0;
}

int Bap::plot()
{
    fs::path script_path = bin_path / "19_makeKneePlots.R";
    fs::path parameteter_file = output_path / (run_name+".bapParam.csv");
    fs::path implicated_barcodes_file = output_path / (run_name+".implicatedBarcodes.csv.gz");
    fs::path barcode_quant_file = output_path / (run_name+".barcodeQuantSimple.csv");
    string command = "Rscript "+script_path.string()+" "+parameteter_file.string()+" "+
                barcode_quant_file.string()+" "+implicated_barcodes_file.string();
    vector<string> cmd_result;
    int cmd_rtn = exec_shell(command.c_str(), cmd_result);
    for (const auto& line : cmd_result)
        if (!line.empty())
            spdlog::error(line);
    if (cmd_rtn == 0)
        spdlog::info("Execute Rscript success");
    else
    {
        spdlog::error("Execute Rscript fail, rtn:{}", cmd_rtn);
        return -1;
    }

    return 0;
}

bool Bap::parseBarcodeList()
{
    ifstream ifs(barcode_list, std::ifstream::in);
    string line;
    int pos = 0;
    while (std::getline(ifs, line))
    {
        _barcode_names.push_back(line);
        _barcode2int[line] = pos++;
    }
    if (_barcode_names.empty()) return false;
    return true;
}

inline string Bap::int2Barcode(int i)
{
    //spdlog::debug("i: {}", i);
    string res = _barcode_names[i >> 16] + _barcode_names[i & 0xFFFF];
    return res;
}