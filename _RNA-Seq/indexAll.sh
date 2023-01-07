#!/usr/bin/env bash
#indexAll.sh

#initialize filepath to contain the directory of sorted bam files
filepath="./bam/"
#initialize sorted to contain the sorted bam files.
sorted=".sorted.bam"

#loop through all the sorted bam files in $filepath
function indexAll {

       for leftInFile in $filepath*$sorted
       do 
           #remove the path from the filename and assign to pathRemoved
           pathRemoved="${leftInFile/$filepath/}"
           #remove the suffix sorted.bam from the pathRemoved
           sampleName="${pathRemoved/$sorted/}"
           # apply samtools to index the files in $filepath$sampleName$sorted
           samtools index $filepath$sampleName$sorted
       done

}

#set output files to store the number of successful output and error output
indexAll 1>indexAll.log 2>indexAll.err &
