/*
 * File: taskflow.cpp
 * Created Data: 2020-7-23
 * Author: fxzhao
 * Contact: <zhaofuxiang@genomics.cn>
 *
 * Copyright (c) 2020 BGI
 */

#include "taskflow.h"

#include <taskflow/taskflow.hpp>  // Taskflow is header-only
#include <iostream>
#include <string>
using std::string;
using namespace std;

int taskflow()
{
    tf::Executor executor;
    tf::Taskflow taskflow;

    taskflow.name("Bap");

    // Step 1: split bam by chr
    auto v = {"chr1", "chr2", "chr3", "chr4"};
    auto [S, T] = taskflow.parallel_for(
          v.begin(),
          v.end(),
          [] (const char* s) {
            std::cout<<"split bam "<<s<<std::endl;
          }
    );

    S.name("Start SplitBam/Assemble/Annotate");
    T.name("End SplitBam/Assemble/Annotate");
    S.for_each_successor([v,s=0] (tf::Task successor) mutable
    {
        successor.name(*(v.begin()+(s++)));
    });


    // Step 2: determine HQ beads
    auto determineHQBeads = taskflow.emplace([] ()
    {
        std::cout<<"Determine HQ beads"<<std::endl;
    }).name("Determine HQ Beads");
    determineHQBeads.succeed(T);

    // Step 3:
    auto [start_cal, end_cal] = taskflow.parallel_for(
        v.begin(),
        v.end(),
        [] (const char* s) {
            std::cout<<"cal chr "<<s<<std::endl;
        }
    );
    start_cal.name("Start Cal");
    end_cal.name("End Cal");
    start_cal.for_each_successor([v,s=0] (tf::Task successor) mutable
    {
        successor.name(*(v.begin()+(s++)));
    });
    start_cal.succeed(determineHQBeads);
    start_cal.succeed(T);

    // Step 4:
    auto barcodeMerge = taskflow.emplace([] ()
    {
        std::cout<<"Determine Barcode Merge"<<std::endl;
    }).name("Determine Barcode Merge");
    barcodeMerge.succeed(end_cal);
    barcodeMerge.succeed(determineHQBeads);

    // Step 5:
    auto [start_reanno, end_reanno] = taskflow.parallel_for(
        v.begin(),
        v.end(),
        [] (const char* s) {
            std::cout<<"reannotate chr "<<s<<std::endl;
        }
    );
    start_reanno.name("Start Reannotate");
    end_reanno.name("End Reannotate");
    start_reanno.for_each_successor([v,s=0] (tf::Task successor) mutable
    {
        successor.name(*(v.begin()+(s++)));
    });
    start_reanno.succeed(barcodeMerge);
    start_reanno.succeed(T);

    // Step 6:
    auto [start_annobam, end_annobam] = taskflow.parallel_for(
        v.begin(),
        v.end(),
        [] (const char* s) {
            std::cout<<"reannotate chr "<<s<<std::endl;
        }
    );
    start_annobam.name("Start Annotate bam");
    end_annobam.name("End Annotate bam");
    start_annobam.for_each_successor([v,s=0] (tf::Task successor) mutable
    {
        successor.name(*(v.begin()+(s++)));
    });
    start_annobam.succeed(end_reanno);
    start_annobam.succeed(T);

    // Step 7
    auto merge_bam = taskflow.emplace([] ()
    {
        std::cout<<"Merge Bam file"<<std::endl;
    }).name("Merge bam");
    merge_bam.succeed(end_annobam);

    auto simple_qc = taskflow.emplace([] ()
    {
        std::cout<<"Simple qc"<<std::endl;
    }).name("Simple qc");
    simple_qc.succeed(end_reanno);

    // Step 8
    auto final_qc = taskflow.emplace([] ()
    {
        std::cout<<"Final qc"<<std::endl;
    }).name("Final qc");
    final_qc.succeed(simple_qc);
    final_qc.succeed(barcodeMerge);

    auto plot = taskflow.emplace([] ()
    {
        std::cout<<"Plot"<<std::endl;
    }).name("Plot");
    plot.succeed(barcodeMerge);
    plot.succeed(determineHQBeads);




    //executor.run(taskflow).wait();

    //taskflow.dump(std::cout);

    return 0;
}