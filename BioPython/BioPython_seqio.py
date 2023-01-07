#! /usr/bin/env python3
# BioPython_seqio.py

from Bio.Seq import Seq
from Bio import SeqIO
import sys
from Bio.Alphabet import generic_dna
   
'''
checking if the arguments is equal to 2
'''
if __name__ == "__main__":
   arg_count = len(sys.argv) - 1
   if arg_count != 2 :
       raise Exception("This script requires 2 arguments: Name of the original FASTA file and the new multi-seq FASTA file ")

   sequence = sys.argv[1]
   output = sys.argv[2]

seqe = []
"""
reading the multi-sequence FASTA file with SeqIO and output the reverse complement of the FASTA file
"""
for record in SeqIO.parse(sequence, "fasta"):
       sequence_1 = record.seq
       sequence_1 = sequence_1.reverse_complement()
       record.seq = sequence_1
       seqe.append(record)
SeqIO.write(seqe, output , "fasta")
