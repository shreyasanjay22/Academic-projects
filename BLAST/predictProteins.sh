#!/usr/bin/env bash
# predictProteins.sh

#specifying the path and version of the software
/usr/local/programs/TransDecoder-5.0.1/TransDecoder.Predict \
#specifying the path of the sequences
-t ../RNA-Seq/trinity_de-novo/Trinity.fasta \
#storing the hits from pfam and blastp
--retain_pfam_hits pfam.domtblout \
--retain_blastp_hits blastp.outfmt6 \
#writing the standdard output and error files
1>predict.log 2>predict.err &
