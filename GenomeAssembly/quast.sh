#!/usr/bin/env bash
#Quast.sh



DIRECTORY="Quast/"
 # Create the output directories
  mkdir -p $DIRECTORY

# assign function for basic analysis of geenomic report and to obtain N50 contigs using quast
function run_quast {
  nice -n19 /usr/bin/quast.py \
  -o $DIRECTORY \
  --threads 6 \
  -s \
  Rhodo/scaffolds.fasta

 }

 run_quast
#write standard output and error
run_quast 1>quast.log 2>quast.err 
