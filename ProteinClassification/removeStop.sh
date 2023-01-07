#usr/bin/env bash
#removeStop.sh


#path of the protein file
#extract the first 500 lines from protein fasta file
#remove stop codon from file
cat ~/BINF6308/BLAST/Trinity.fasta.transdecoder.pep | head -n500 | sed 's/*//g' > ./proteins.fasta
