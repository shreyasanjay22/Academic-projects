#!/usr/bin/env bash
#mergeAll.sh

#redirecting the output to bamIn.txt to produce list
ls bam/Aip*.sorted.bam > bamIn.txt
#samtools to merge the files and the parameter -b is used to pass in bamIn.txt
samtools merge -b bamIn.txt bam/AipAll.bam \
    #writing standard outputs and errors in merge.log and merge.err files and running in background
    1>bam/merge.log 2>bam/merge.err &
