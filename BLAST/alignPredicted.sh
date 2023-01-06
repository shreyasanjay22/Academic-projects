#!/usr/bin/env bash
#alignPredicted.sh

#specifying the blast program and query sequence
blastp -query ./Trinity.fasta.transdecoder.pep  \
#specifying the database for alignment
-db swissprot \
 #specifying the number of threads and setting the  E-value 
-evalue 1e-10 -num_threads 4 \
#arranging the output in tabular format with the following fields
-outfmt "6 qseqid sacc qlen slen length nident pident evalue stitle" \
#writing the standard output and error files
    1>alignPredicted.txt 2>alignPredicted.err &
