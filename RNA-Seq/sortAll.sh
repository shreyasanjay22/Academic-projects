#!/usr/bin/env bash
#sortAll.sh

#initialize variable to contain the files from sam directory
filepath="./sam/"
#initialize variable to contain the .sam format
left=".sam"
#initialize variable to contain the format of the output
newFormat=".sorted.bam"
out="bam/"

#create output directory
mkdir -p $out

#loop through all the left-reads in $filepath
 function sortAll {
     
    for leftInFile in $filepath*$left
    do
        #remove the path from the filename and assign to pathRemoved
        pathRemoved="${leftInFile/$filepath/}"
        #remove the left-read suffix from $pathremoved and assign to sampleNAme
        sampleName="${pathRemoved/$left/}"
        #Use samtools to sort the sam files
        samtools sort \
         #parameter for files to be used to sort 
        $filepath$sampleName$left \
            #write the final sorted output to $output$sampleNAme$newFormat
        -o $out$sampleName$newFormat
    done
}

#set files for the success output and error outputs
sortAll 1>sortAll.log 2>sortAll.err &
