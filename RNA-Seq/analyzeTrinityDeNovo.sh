#!/usr/bin/env bash
# analyzeTrinity.sh

#running TrinityStats.pl with the assembled transcriptome as input
/usr/local/programs/trinityrnaseq-Trinity-v2.8.4/util/TrinityStats.pl \
    #directory and file for output
    trinity_de-novo/Trinity.fasta
