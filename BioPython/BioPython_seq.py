#! /usr/bin/env python3
# BioPython_seq.py

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

"""
creating a SeqRecord object
"""

simple_seq= Seq("aaaatgggggggggggccccgtt", generic_dna)

simple_seq_r = SeqRecord(simple_seq, id = '#12345', description = 'example 1')

print(simple_seq_r)
print(simple_seq_r.seq)
print(simple_seq_r.id)
print(simple_seq_r.description)
"""
writing the sequence file in GenBack format
"""
SeqIO.write(simple_seq_r,"BioPython_seq.gb", "genbank")
