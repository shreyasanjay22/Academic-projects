#!/usr/bin/env bash
#blastPep.sh

#specifying the type of blast program and query sequence
blastp -query Trinity.fasta.transdecoder_dir/longest_orfs.pep  \
 #specifying the database and the target alignment
-db swissprot  -max_target_seqs 1 \
 #specifying the threads and the format for output. Writing the standard output file and error files
-outfmt 6 -evalue 1e-5 -num_threads 4 1> blastp.outfmt6 \
2>blastp.err &
