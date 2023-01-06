#!/usr/bin/env bash
# buildIIT.sh

#to give less priority than the critical function
nice -n19 iit_store \ 
-G /scratch/AiptasiaMiSeq/\
GCA_001417965.1_Aiptasia_genome_1.1_genomic.gff \
#set the output file for storing output
-o AiptasiaGmapIIT \

#set files to write the success output and error outputs
1>buildIIT.log 2>buildIIT.err &

