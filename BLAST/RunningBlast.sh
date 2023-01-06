#!/usr/bin/env bash
# RunningBlast

#specifying the query and the type of blast 
blastx -query /scratch/SampleDataFiles/Seq.fasta \
#specifying the database
-db /blastDB/swissprot \
#specifying the number of alignments 
-num_alignments 5 \
 #obtaining the output in tabular format
-outfmt 6
