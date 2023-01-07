#usr/bin/env bash
#runIPS.sh 

/usr/local/programs/interproscan-5.26-65.0/interproscan.sh \ #use interproscan with the inpput path
 -i proteins.fasta \ #input file
 -o proteins.tsv \ #convert to tsv output file
 -f TSV #using tsv format
