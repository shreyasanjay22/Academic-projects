#! /usr/bin/env python3
# BioPython_genbank.py

from Bio import Entrez
from Bio import SeqIO

# entering the email to download the sequence
Entrez.email = "shirodkar.sh@northeastern.edu"

# using Entrez.efetch to obatain the GenBank entry with IDs
with Entrez.efetch(
        db="nucleotide", rettype="gb", retmode = "text", id= "515056,J01673.1"
) as handle:
    for seq_record in  SeqIO.parse(handle,"gb"):
        print("Name of the Sequence {}".format(seq_record.name))
        print(seq_record.seq)
        for k in seq_record.features:
            print("Sequence type {} \n Sequence location {} \n Sequence strand {}".format(k.type,k.location,k.strand))

