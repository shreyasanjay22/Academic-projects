#!/usr/bin/env python3
# translate_APOE.py

# importing BIo, Bio.Seq, Bio.Alphabet and re
from Bio import SeqIO
from Bio.Seq import Seq
import re
from Bio.Alphabet import generic_dna

# initializing variable to store the sequence
seqe = []


"""
Translating the sequences from nucleotide sequences to amino acid sequences using BioPython
"""
for record in SeqIO.parse("/home/shirodkar.sh/BINF6309/MSA/APOE_refseq_transcript.fasta","fasta"):
        dna = record.seq
        # carrying out transcription of the sequences
        rna = dna.transcribe()
        # obtaining the longest ORF
        orf = re.search('AUG([AUGC]{3})+?(UAA|UAG|UGA)',str(rna)).group()
        # carrying out translation of the sequences
        protein = Seq(orf).translate()
        record.seq = protein
        # appending the sequences
        seqe.append(record)

# writing the output in a fasta file
with open("apoe_aa.fasta", "w") as output_handle:
        SeqIO.write(seqe, output_handle, "fasta")
