#! /usr/bin/env python3
# sliding_window_fasta.py

import re
import sys
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio import SeqIO

def sliding_window(k,string):
    """ 
    Returns the list of all Kmers in the string
    """
    kmers = []
    end = len(string) - k + 1
    for start in range (0, end):
        kmers.append(string[start:start + k ])
    return kmers


def gc_content(dna):
    """
    Returns the fractions of the GCs in the string
    """
    dna = dna.lower()
    gc = 0
    for nucleotide in dna:
        if nucleotide in ['g','c']:
           gc += 1
    return gc/len(dna)


if __name__ == "__main__":
    # checks if two arguments are passed
  arg_count =len(sys.argv) - 1
  if arg_count != 2:
         raise Exception ("The script requires exactly two arguments")
  elif int(sys.argv[1]) <= 0:
          sys.exit("Value of Kmer must be greater than zero")
  # stores the arguments in the variables
  k = int(sys.argv[1])
  files = sys.argv[2]

# reads a fasta file and separates the headers from the strings
  string = " "
  for record in SeqIO.parse(files,"fasta"):
      print(record.id)
      string = record.seq

# prints the k-mer followed by the GCs fraction
  for i in sliding_window(k,string):
      x = 0
      x = gc_content(i)
      print("{}\t{:.2f}".format(i ,x))
   
