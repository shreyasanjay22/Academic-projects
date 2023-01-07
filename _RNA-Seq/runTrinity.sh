#!/usr/bin/env bash
#runTrinity.sh

#path of the trinity assembler 
nice -n19 /usr/local/programs/trinityrnaseq-Trinity-v2.8.4/Trinity \
    #indicates that this is a genome-guided assembly and the merged bam file
    --genome_guided_bam bam/AipAll.bam \
        #specifies the maximum separation distance trinity will allow for segments of transcript to span introns
    --genome_guided_max_intron 10000 \
        #maximum memory to be used by trinity for the assembly
    --max_memory 10G --CPU 4 \
        #files for standard outputs and error and running the program in background
    1>trinity.log 2>trinity.err &
