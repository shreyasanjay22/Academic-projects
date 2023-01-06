#!/usr/bin/env bash
#alignAll.sh

#initialize variable to contain the trimmed fastq files
filepath="./Paired/"
#initialize variables to contain left and right reads
left=".R1.fastq"
right=".R2.fastq"
#initialize variable to store the format of the output file
newFormat=".sam"
output="sam/"

#create the output directory
mkdir -p $output
#loop through all the left-reads files in $filepath

function alignReads {
   for leftInFile in $filepath*$left
           do
               #Remove the path from the filename and assign to pathRemoved
             pathRemoved="${leftInFile/$filepath/}"
             #Remove the left reads suffix from pathRemoved and assign to sampleName
             sampleName="${pathRemoved/$left/}"
             #print sampleName
             echo $sampleName
             #inform  server to give lower priority than the critical server function
             nice -n19 gsnap \
              #instructs gsnap to create sam alignment program
             -A sam \
             #indicates directory to use working directory for the program
             -D . \
             #indicates the name of the database used for aligning
             -d AiptasiaGmapDb \
             -s AiptasiaGmapIIT.iit \
             #files to be used to create the sam alignment outputs
             $filepath$sampleName$left \
             $filepath$sampleName$right \
             #redirect sam alignment information to $output$sampleName$newFormat
             1>$output$sampleName$newFormat
     done 
 }

 #set ouput files to store success outputs and error outputs 
 alignReads 1>alignAll.log 2>alignAll.err &
