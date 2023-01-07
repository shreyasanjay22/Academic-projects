#!/usr/bin/env python3
#interleave.py
#Import Seq, SeqRecord and SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
#import itertools to take a slice of list
import itertools
#read both the files where the first parameter is the file location and the second is the file type
leftReads =SeqIO.parse("/scratch/AiptasiaMiSeq/fastq/Aip02.R1.fastq", "fastq")
rightReads=SeqIO.parse("/scratch/AiptasiaMiSeq/fastq/Aip02.R2.fastq", "fastq")
#initialise an empty list  to store the sequences
sequence= [ ]
#iterate over the two files
for LeftSeq,RightSeq in zip(leftReads,rightReads):
    #append the left and right sequences to the initialised list
      sequence.append(LeftSeq)
      sequence.append(RightSeq)
      #open the file and pass the write parameter to the file
with open("Interleaved.fasta","w") as output:
    #storing the appended list from above to the file in fasta format
    SeqIO.write(sequence,output,"fasta")
    #initialising the empty list again to avoid repeating any elements of the list
    sequence=[ ]
