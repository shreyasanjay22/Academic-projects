#!/usr/bin/env bash
# AipBuild.sh
#-D . indicates  to use the working directory to build database
nice -n19 gmap_build -D . \
#indicates the name of the database
-d AiptasiaGmapDb \
#set location of the Aiptasia genome in fasta format
/scratch/AiptasiaMiSeq/\
GCA_001417965.1_Aiptasia_genome_1.1_genomic.fna \
#Write the success and error output and run in background
1>AipBuild.log 2>AipBuild.err &
