#!/usr/bin/env bash
# trinityDeNovo1.sh

# store left reads  as $leftReads
leftReads="$(ls -q Paired/Aip*.R1.fastq)"
# Store echo of $leftReads as $leftReads 
leftReads=$(echo $leftReads)
# Replace spaces in the list with comma
leftReads="${leftReads// /,}"
#print leftReads
echo $leftReads
# store right reads as $rightReads
rightReads="$(ls -q Paired/Aip*.R2.fastq)"
# Store echo of $rightReads as $rightReads to get rid of line breaks
rightReads=$(echo $rightReads)
# Replace spaces in the list  with comma
rightReads="${rightReads// /,}"
#print rightReads
echo $rightReads
