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

#include <spdlog/spdlog.h>

#include <taskflow/taskflow.hpp>

#include <fstream>
#include <algorithm>


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

void Bap::extractBedPE(const BamRecord b1, const BamRecord b2, vector<Bedpe>& bedpes)
{
    // Initialize BEDPE variables
    string chrom1, chrom2, strand1, strand2;
    int start1, start2, end1, end2;
    start1 = start2 = end1 = end2 = -1;
    chrom1 = chrom2 = strand1 = strand2 = ".";
    int min_map_quality = 0;

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
    if (chrom1 > chrom2 || ((chrom1 == chrom2) && (start1 > start2)))
    {
        swap(chrom1, chrom2);
        swap(start1, start2);
        swap(end1, end2);
        swap(strand1, strand2);
    }

    // Report BEDPE using min mapping quality
    if (isMapped(b1) && isMapped(b2))
        min_map_quality = min(b1->core.qual, b2->core.qual);
    
    // Filter by max insert size and mapping quality
    if ((end2 - start1) < MAX_INSERT && min_map_quality >= mapq)
    {
        if (strand1 == "+")
        {
            start1 += 4;
            end2 += 4;
        }
        else
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
            bedpe.qname = getQname(b1);
            bedpe.barcode = barcode;
            bedpes.push_back(bedpe);
        }
    }
}

Bap::Bap(string input_bam, string output_path, string barcode_tag, int mapq, int cores, string run_name, bool tn5,
        double min_barcode_frags, double min_jaccard_index, string ref, string mito_chr, string bed_genome_file,
        string blacklist_file, string trans_file, bool species_mix, string bin_path) :
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
        bin_path(bin_path) 
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
            cout<<"Failed to create directory: "<<temp_bam_path.string()<<endl;
            return -1;
        }
    }
    spdlog::info("Prepare time(s):{:.2f}", timer.toc(1000));

    // Execute taskflow
    taskflow();
    spdlog::info("Taskflow time(s):{:.2f}", timer.toc(1000));

    // Clean up
    //fs::remove_all(temp_bam_path);
    spdlog::info("Finish taskflow.");

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
    spdlog::debug("Bam contigs num:{}", contigs.size());

    for (auto& p : contigs)
        _contig_names.push_back(p.first);

    // Verify that the supplied reference genome and the bam have overlapping chromosomes
    // TODO(fxzhao)

    _bedpes_by_chr.resize(contigs.size());

    set<int> used_chrs;
    auto bed_chrs = parseChrsFromBedFile();
    cout<<"bed_chrs size: "<<bed_chrs.size()<<endl;
    auto [S, T] = taskflow.parallel_for(
        0,
        int(contigs.size()),
        1,
        [&] (int chr_id) {
            if (bed_chrs.count(contigs[chr_id].first) == 0) return;
            if (!samReader->QueryByContig(chr_id)) return;
            used_chrs.insert(chr_id);
            Bap::splitBamByChr(chr_id);
        }
    );

    S.name("Start SplitBam/Assemble/Annotate");
    T.name("End SplitBam/Assemble/Annotate");
    S.for_each_successor([contigs, s=0] (tf::Task successor) mutable
    {
        successor.name((contigs.begin()+(s++))->first);
    });

    // Step 2: determine hight quality beads
    auto determine_hq_beads = taskflow.emplace([&] ()
    {
        Bap::determineHQBeads();
    }).name("Determine HQ Beads");
    determine_hq_beads.succeed(T);

    // Step 3: calculate stat by chr
    _total_bead_cnts.resize(contigs.size());
    _total_nc_cnts.resize(contigs.size());
    cout<<"_total_bead_cnts size: "<<_total_bead_cnts.size()<<endl;
    // cout<<"contig size: "<<_bedpes_by_chr.size()<<endl;
    // for (int i = 0; i < _bedpes_by_chr.size(); ++i)
    //     if (_bedpes_by_chr[i].size() > 0)
    //         cout<<"has data chr: "<<i<<endl;
    auto [start_cal, end_cal] = taskflow.parallel_for(
        0,
        int(contigs.size()),
        1,
        [&] (int chr_id) {
            if (used_chrs.count(chr_id) == 0) return;
            Bap::computeStatByChr(chr_id);
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
        Bap::determineBarcodeMerge();
    }).name("Determine Barcode Merge");
    barcode_merge.succeed(end_cal);
    barcode_merge.succeed(determine_hq_beads);

    // Step 5: reannotate frag data and get summary stats by chr
    _keep_qnames.resize(contigs.size());
    _dup_frags.resize(contigs.size());
    _frag_stats.resize(contigs.size());
    auto [start_reanno, end_reanno] = taskflow.parallel_for(
        0,
        int(contigs.size()),
        1,
        [&] (int chr_id) {
            if (used_chrs.count(chr_id) == 0) return;
            Bap::reannotateFragByChr(chr_id);
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
        0,
        int(contigs.size()),
        1,
        [&] (int chr_id) {
            if (used_chrs.count(chr_id) == 0) return;
            Bap::annotateBamByChr(chr_id);
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
        std::vector< std::string > bam_files;
        for (int chr_id = 0; chr_id < contigs.size(); ++chr_id)
        {
            if (used_chrs.count(chr_id) == 0) continue;
            fs::path tmp_bam_file = temp_bam_path / (contigs[chr_id].first+".bam");
            if (fs::exists(tmp_bam_file))
                bam_files.push_back(tmp_bam_file.string());
        }
        fs::path output_bam_file = output_path / (run_name + ".bam");
        int cmd_rtn = bam_cat(bam_files, nullptr, output_bam_file.c_str(), nullptr, 0);
        if (cmd_rtn == 0)
            spdlog::info("Merge bam file success");
        else
            spdlog::info("Merge bam file fail, rtn:{}", cmd_rtn);
    }).name("Merge bam");
    merge_bam.succeed(end_annobam);

    // Step 8: merge fragment files
    auto merge_frags = taskflow.emplace([&] ()
    {
        for (int chr_id = 0; chr_id < contigs.size(); ++chr_id)
        {
            if (used_chrs.count(chr_id) == 0) continue;
            _final_frags.insert(_final_frags.end(), _dup_frags[chr_id].begin(), _dup_frags[chr_id].end());
            _dup_frags[chr_id].clear();
            _dup_frags[chr_id].shrink_to_fit();
        }
        fs::path out_frag_file = output_path / (run_name+".fragment.tsv");
        FILE* out_frag;
        out_frag = fopen(out_frag_file.c_str(), "w");
        char new_line = '\n';
        for (auto& l : _final_frags)
        {
            fwrite(l.c_str(), 1, l.size(), out_frag);
            fwrite(&new_line, 1, 1, out_frag);
        }
        fclose(out_frag);
        std::cout<<"merge frags"<<std::endl;
    }).name("Merge fragments");
    merge_frags.succeed(end_reanno);

    // Step 9: simple qc
    auto simple_qc = taskflow.emplace([&] ()
    {
        map<string, pair<int,int>> nuclear;
        int mito_pos = -1;
        for (int chr_id = 0; chr_id < contigs.size(); ++chr_id)
        {
            if (used_chrs.count(chr_id) == 0) continue;
            cout<<"contig name: "<<chr_id<<" "<<_contig_names[chr_id]<<" "<<mito_chr<<endl;
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

        std::cout<<"Simple qc"<<std::endl;
    }).name("Simple qc");
    simple_qc.succeed(end_reanno);

    // Step 10: final qc
    auto final_qc = taskflow.emplace([&] ()
    {
        Bap::finalQC();
        std::cout<<"Final qc"<<std::endl;
    }).name("Final qc");
    final_qc.succeed(barcode_merge);
    final_qc.succeed(simple_qc);
    final_qc.succeed(merge_frags);

    // Step 11: plot
    auto plot = taskflow.emplace([&] ()
    {
        Bap::plot();
        std::cout<<"Plot"<<std::endl;
    }).name("Plot");
    plot.succeed(barcode_merge);
    plot.succeed(determine_hq_beads);

    executor.run(taskflow).wait();

    return 0;
}

int Bap::splitBamByChr(int chr_id)
{
    std::unique_ptr<SamReader> sr = SamReader::FromFile(input_bam);
    if (!sr->QueryByContig(chr_id)) return 0;
    spdlog::debug("Call splitBamByChr:{}", chr_id);

    // Map {qname:barcode} used in later
    map<string, string> qname2barcodes;

    auto contigs = sr->getContigs();
    string chr_str = contigs[chr_id].first;
    BamRecord bam_record1 = bam_init1();
    BamRecord bam_record2 = bam_init1();
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
    map<string, BamRecord> pe_dict;
    map<string, BamRecord>::iterator it;
    auto& bedpes = _bedpes_by_chr[chr_id];
    while (true)
    {
        BamRecord bam_record = bam_init1();
        if(sr->next(bam_record, FLAG))
        {
            string qname = getQname(bam_record);
            it = pe_dict.find(qname);
            if (it == pe_dict.end())
            {
                pe_dict[qname] = bam_record;
            }
            else
            {
                extractBedPE(bam_record, it->second, bedpes);
                bam_destroy1(it->second);
                pe_dict.erase(it);
                bam_destroy1(bam_record);
            }
        }
        else
            break;
    }
    cout<<"pe_dict size: "<<pe_dict.size()<<endl;
    for (auto& p : pe_dict)
        bam_destroy1(p.second);
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
    cout<<"frags size: "<<bedpes.size()<<endl;

    bam_destroy1(bam_record1);
    bam_destroy1(bam_record2);

    spdlog::debug("chr:{} frags size:{}", chr_id, bedpes.size());
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
    map<string, int> bead_quant;
    for (auto& b : bedpes)
    {
        // Filter for fragments overlapping the blacklist
        int start = b.start;
        int end = b.end;
        MyInterval query_range {start, end};
        
        const auto& res = mytree.query(query_range);
        // Discard the fragments overlapping the blacklist
        if (res.begin() != res.end()) continue;

        string barcode = b.barcode;
        string key = to_string(start)+'\t'+to_string(end)+barcode;
        if (pcr_dup.count(key) == 0)
        {
            pcr_dup.insert(key);
            ++bead_quant[barcode];
        }
    }
    
    // Devel
    cout<<"bead_quant size:"<<bead_quant.size()<<endl;
    int total_cnt = 0;
    for (auto& b : bead_quant)
        total_cnt += b.second;
    cout<<"bead_quant total count:"<<total_cnt<<endl;

    // Merge bead quant of all chrs
    std::lock_guard<std::mutex> guard(_merge_chr_mutex);
    for (auto& b : bead_quant)
        _total_bead_quant[b.first] += b.second;

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
    {
        string s = b.first + "," + to_string(b.second) + "\n";
        fwrite(s.c_str(), 1, s.size(), out_bead_quant);
    }
    fclose(out_bead_quant);

    // TODO(fxzhao): change it to cpp
    // TODO(fxzhao): other conditions
    // Call R script to calculate bead threshold
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
        spdlog::info("Execute Rscript fail, rtn:{}", cmd_rtn);
        return -1;
    }

    // Define the set of high-quality bead barcodes
    auto paras = parseBeadThreshold(filename.string()+"_kneeValue.txt");
    double bead_threshold = paras.first, call_threshold = paras.second;

    fs::path paras_file = output_path / (run_name+".bapParam.csv");
    ofstream ofs(paras_file.string(), std::ofstream::out);
    ofs << "bead_threshold_nosafety,"<<call_threshold<<endl;
    ofs << "bead_threshold,"<<bead_threshold<<endl;
    ofs.close();
				
    // Devel
    //bead_threshold = 0;

    for (auto& p : _total_bead_quant)
    {
        if (p.second >= bead_threshold)
            _hq_beads.insert(p.first);
    }
    cout<<"_total_bead_quant size: "<<_total_bead_quant.size()<<endl;
    cout<<"bead_threshold: "<<bead_threshold<<endl;
    return 0;
}

int Bap::computeStatByChr(int chr_id)
{
    spdlog::debug("in computeStatByChr: {}", chr_id);
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
    cout<<"_hq_beads size: "<<_hq_beads.size()<<endl;
    cout<<"frags data size: "<<frags_data.size()<<endl;
    cout<<"frags pos size: "<<frags_pos.size()<<endl;

    // Quantify NC + export
    auto& cnts = _total_nc_cnts[chr_id];
    for (auto& p : dict)
    {
        cnts[p.second] += p.second;
    }
    cout<<"cnts size: "<<cnts.size()<<endl;
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
    map<int, vector<string>> overlap_start, overlap_end;
    for (auto& pos : frags_pos)
    {
        if (pos == -1) continue;

        auto& bedpe = frags_data[pos];
        int start = bedpe.start;
        int end = bedpe.end;
        string barcode = bedpe.barcode;
        string key = to_string(start)+'\t'+to_string(end)+barcode;
        if (uniq_frags.count(key) != 0) continue;

        uniq_frags.insert(key);
        overlap_start[start].push_back(barcode);
        overlap_end[end].push_back(barcode);
    }

    // Double to consider left and right inserts
    map<string, int> bead_cnts;
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
                        ++bead_cnts[v[i]+","+v[j]];
                    else
                        ++bead_cnts[v[j]+","+v[i]];
                }
            }
        }
    }
    cout<<"test "<<bead_cnts.size()<<endl;
    bead_cnts.swap(_total_bead_cnts[chr_id]);
    cout<<"_total_bead_cnts size: "<<_total_bead_cnts.size()<<endl;
    cout<<"test "<<bead_cnts.size()<<endl;
    cout<<"test "<<_total_bead_cnts[chr_id].size()<<endl;

    // Devel
    int line_num = 0, total_num = 0;
    for (auto& p : bead_cnts)
    {
        int count = p.second;
        if (count >= regularize_threshold)
        {
            ++line_num;
            total_num += count;
        }
    }
    cout<<"compute stat by chr: "<<chr_id<<" "<<line_num<<" "<<total_num<<endl;

    return 0;
}

int Bap::determineBarcodeMerge()
{
    // Merge all barcode count
    map<string, int> sum_dt;
    for (auto& m : _total_bead_cnts)
        for (auto& p : m)
            sum_dt[p.first] += p.second;
    cout<<"sum_dt size: "<<sum_dt.size()<<endl;
    // Only consider merging when Tn5 is the same
    if (tn5)
    {
        // TODO(fxzhao): implement tn5 option
    }
     
    // Filter and calculate nBC
    vector<pair<string, int>> nBC;
    map<string, int> count_dict;
    for (auto& p : _total_bead_quant)
    {
        if (_hq_beads.count(p.first) == 0) continue;
        nBC.push_back({p.first, p.second});
        count_dict[p.first] = p.second * 2;
    }

    // Release memory of hq beads
    _hq_beads.clear();

    // Sort by count
    std::sort(nBC.begin(), nBC.end(), [](pair<string, int>& a, pair<string, int>& b) {
        return a.second > b.second;
    });
    cout<<"nBC size:"<<nBC.size()<<endl;

    // Calculate jaccard frag
    vector<pair<string, float>> ovdf;
    for (auto& p : sum_dt)
    {
        vector<string> vec_s = split_str(p.first, ',');
        if (vec_s.size() != 2) continue;

        int N_both = p.second;
        int N_barc1 = count_dict[vec_s[0]];
        int N_barc2 = count_dict[vec_s[1]];
        float jaccard_frag = round((N_both)/(N_barc1+N_barc2-N_both+0.05), 5);
        if (jaccard_frag <= 0.0) continue;
        ovdf.push_back({p.first, jaccard_frag});
    }
    cout<<"ovdf size: "<<ovdf.size()<<endl;
    // Sort by jaccard_frag
    std::sort(ovdf.begin(), ovdf.end(), [](pair<string, float>& a, pair<string, float>& b) {
        return a.second > b.second;
    });

    fs::path paras_file = output_path / (run_name+".bapParam.csv");
    ofstream ofs(paras_file.string(), std::ofstream::out | std::ofstream::app);
    // Call knee if we need to
    if (min_jaccard_index == 0.0)
    {
        // TODO(fxzhao): implement jaccard option
        // fs::path script_path = bin_path / "10b_knee_execute.R";
        // string command = "Rscript "+script_path.string()+" "+filename.string()+" 1 V2";
        // vector<string> cmd_result;
        // int cmd_rtn = exec_shell(command.c_str(), cmd_result);
        // for (const auto& line : cmd_result)
        //     if (!line.empty())
        //         spdlog::error(line);
        // if (cmd_rtn == 0)
        //     spdlog::info("Execute Rscript success");
        // else
        // {
        //     spdlog::info("Execute Rscript fail, rtn:{}", cmd_rtn);
        //     return -1;
        // }

        // // Define the set of high-quality bead barcodes
        // double bead_threshold = parseBeadThreshold(filename.string()+"_kneeValue.txt");

        ofs << "bead_threshold_nosafety,"<<min_jaccard_index<<endl;
    }
    ofs << "jaccard_threshold,"<<min_jaccard_index<<endl;
    ofs.close();

    // Export the implicated barcodes
    FILE * tbl_out;
    fs::path implicated_barcode_file = output_path / (run_name+".implicatedBarcodes.csv"); 
    tbl_out = fopen(implicated_barcode_file.c_str(), "w");
    string header = "barc1,barc2,N_both,N_barc1,N_barc2,jaccard_frag,merged\n";
    fwrite(header.c_str(), 1, header.size(), tbl_out);
    for (size_t i = 0; i < ovdf.size(); ++i)
    {
        auto& p = ovdf[i];
        string s = p.first + ",";
        s += to_string(sum_dt[p.first]) + ",";
        vector<string> vec_s = split_str(p.first, ',');
        if (vec_s.size() != 2) continue;
        int N_barc1 = count_dict[vec_s[0]];
        int N_barc2 = count_dict[vec_s[1]];
        s += to_string(N_barc1) + "," + to_string(N_barc2) + ",";
        s += f2str(p.second, 5) + ",";
        s += p.second > min_jaccard_index ? "TRUE" : "FALSE";
        s += "\n";

        fwrite(s.c_str(), 1, s.size(), tbl_out);
    }
    fclose(tbl_out);

    // Filter based on the min_jaccard_index 
    // and prepare dict data
    map<string, set<string>> bc1, bc2;
    for (size_t i = 0; i < ovdf.size(); ++i)
    {
        auto& p = ovdf[i];
        if (p.second <= min_jaccard_index) continue;
        vector<string> vec_s = split_str(p.first, ',');
        if (vec_s.size() != 2) continue;
        bc1[vec_s[0]].insert(vec_s[1]);
        bc2[vec_s[1]].insert(vec_s[0]);
    }

    // Guess at how wide we need to make the barcodes to handle leading zeros
    int guess = ceil(log10(nBC.size()));
    // Map from barcode to position in nBC
    map<string, int> bar2pos;
    for (size_t i = 0; i < nBC.size(); ++i)
    {
        bar2pos[nBC[i].first] = i;
    }

    // Loop through and eat up barcodes
    int idx = 1;
    for (size_t i = 0; i < nBC.size(); ++i)
    {
        string barcode = nBC[i].first;
        if (barcode.empty()) continue;

        set<string> barcode_combine;
        barcode_combine.insert(barcode);

        // Find friends that are similarly implicated and append from Barcode 1
        if (bc1.count(barcode) != 0)
        {
            barcode_combine.insert(bc1[barcode].begin(), bc1[barcode].end());
        }
        // Find friends that are similarly implicated and append from Barcode 2
        if (bc2.count(barcode) != 0)
        {
            barcode_combine.insert(bc2[barcode].begin(), bc2[barcode].end());
        }

        // If user species one to one, then only remove that one barcode
        if (one_to_one) barcode_combine = {barcode};

        // Make a drop barcode and save our progress
        stringstream ss;
        if (!tn5)
        {
            ss << run_name <<"_BC"<<std::setfill ('0') << std::setw (guess) << idx<<
                "_N"<<std::setw (2)<<barcode_combine.size();
        }   
        else
        {
            // TODO(fxzhao): implement this
        }
        for (auto& b : barcode_combine)
        {
            if (bar2pos.count(b) != 0)
            {
                _drop_barcodes[b] = ss.str();
                nBC[bar2pos[b]].first = "";
            }
        }
        ++idx;
    }
    
    FILE * bt;
    fs::path barcode_translate_file = output_path / (run_name+".barcodeTranslate.tsv");
    bt = fopen(barcode_translate_file.c_str(), "w");
    for (auto& p : _drop_barcodes)
    {
        string s = p.first + "\t" + p.second+"\n";
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
    set<string> qname_dup;
    // Store n_total and n_unique as pair
    map<string, pair<int, int>>& merge_ss = _frag_stats[chr_id];
    string chr = _contig_names[chr_id];
    auto& frags_data = _bedpes_by_chr[chr_id];
    map<string, string>::iterator it;
    for (auto& bedpe : frags_data)
    {
        // Filter for eligible barcodes
        string barcode = bedpe.barcode;
        it = _drop_barcodes.find(barcode);
        if (it == _drop_barcodes.end()) continue;

        string cell_barcode = it->second;
        string s = chr+"\t"+to_string(bedpe.start)+"\t"+to_string(bedpe.end)+"\t"+cell_barcode;
        if (pcr_dup.count(s) == 0)
        {
            pcr_dup.insert(s);
            qname_dup.insert(bedpe.qname);
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

    cout<<"chr "<<chr_id<<" keep_qname size: "<<_keep_qnames[chr_id].size()<<" dup_frags size: "<<_dup_frags[chr_id].size()<<endl;

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
    spdlog::debug("Call annotateBamByChr:{}", chr_id);

    fs::path output_bam_file = temp_bam_path / (_contig_names[chr_id] + ".bam");
    htsFile* out = hts_open(output_bam_file.c_str(), "wb");
    auto header = sr->getHeader();
    [[maybe_unused]] int ret = sam_hdr_write(out, header);

    BamRecord b = bam_init1();
    string bead_bc, drop_bc;
    map<string, string>::iterator it;
    // Iterate through bam
    while (sr->next(b))
    {
        if (!getTag(b, barcode_tag.c_str(), bead_bc)) continue;
        it = _drop_barcodes.find(bead_bc);
        if (it == _drop_barcodes.end()) continue;
        drop_bc = it->second;

        // Handle droplet barcodes that we want to consider writing out
        string qname = bam_get_qname(b);
        if (keep_reads.count(qname) == 0) continue;
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
                    to_string(l.mito_total)+","+to_string(l.mito_uniq)+","+to_string(l.dup_proportion)+","+
                    to_string(l.library_size)+",";
        auto& sd = it->second;
        s += to_string(sd.mean_insert_size)+","+to_string(sd.median_insert_size)+","+
            to_string(sd.tss_proportion)+","+to_string(sd.frip)+"\n";
        fwrite(s.c_str(), 1, s.size(), qc_out);
    }
    fclose(qc_out);

    return 0;
}

int Bap::plot()
{
    fs::path script_path = bin_path / "19_makeKneePlots.R";
    fs::path parameteter_file = output_path / (run_name+".bapParam.csv");
    fs::path implicated_barcodes_file = output_path / (run_name+".implicatedBarcodes.csv");
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
        spdlog::info("Execute Rscript fail, rtn:{}", cmd_rtn);
        return -1;
    }

    return 0;
}