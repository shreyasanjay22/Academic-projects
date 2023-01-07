#!/usr/bin/env bash
# trinityDeNovo.sh


# Get the list of left reads and store as $leftReads
leftReads="$(ls -q Paired/Aip*.R1.fastq)"
# Store echo of $leftReads as $leftReads to get rid of line breaks
leftReads=$(echo $leftReads)
# Replace spaces in the list of reads with comma
leftReads="${leftReads// /,}"
echo $leftReads
# Get the list of right reads and store as $rightReads
rightReads="$(ls -q Paired/Aip*.R2.fastq)"
# Store echo of $rightReads as $rightReads to get rid of line breaks
rightReads=$(echo $rightReads)
# Replace spaces in the list of reads with comma
rightReads="${rightReads// /,}"
echo $rightReads

#path of trinity assembly
nice -n19 /usr/local/programs/trinityrnaseq-Trinity-v2.8.4/Trinity \
    --seqType fq \
    #specifying an output directory to avoid writing over genome-guided assembly
    --output trinity_de-novo \
   #specify the max memory to be used for the assembly
    --max_memory 10G --CPU 4 \
    # passing variable $leftReads and $rightReads for -left and -right
    --left $leftReads \
    --right $rightReads \
    #writing outputs and errors and running in background
    1>trinity_dn.log 2>trinity_dn.err &
