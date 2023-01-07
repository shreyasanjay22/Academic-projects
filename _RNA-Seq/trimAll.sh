#!/usr/bin/env bash
# trimAll.sh

# Initialize variable to contain the directory of un-trimmed fastq files 
fastqPath="/scratch/AiptasiaMiSeq/fastq/"
# Initialize variable to contain the suffix for the left reads and right reads
leftSuffix=".R1.fastq"
rightSuffix=".R2.fastq"
pairedOutPath="Paired/"
unpairedOutPath="Unpaired/"
# Create the output directories
mkdir -p $pairedOutPath
mkdir -p $unpairedOutPath

# Loop through all the left-read fastq files in $fastqPath
function trimAll {
    for leftInFile in $fastqPath*$leftSuffix
            do
              # Remove the path from the filename and assign to pathRemoved
                pathRemoved="${leftInFile/$fastqPath/}"
              # Remove the left-read suffix from $pathRemoved and assign to suffixRemoved
                sampleName="${pathRemoved/$leftSuffix/}"
              # Print $sampleName to see what it contains after removing the path
                echo $sampleName
               #set the location of the java command with paratmeter PE to indicate pair-end reads. Set thread to 1 to indicate 1 server thread
                 nice -n19 java -jar /usr/local/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
                 -threads 1 -phred33 \
                   # parameter for left and right read file 
                    $fastqPath$sampleName$leftSuffix \
                    $fastqPath$sampleName$rightSuffix \
                   # parameter indicating the  output, Paired and unpaired files.
                    $pairedOutPath$sampleName$leftSuffix \
                    $unpairedOutPath$sampleName$leftSuffix \
                    $pairedOutPath$sampleName$rightSuffix \
                    $unpairedOutPath$sampleName$rightSuffix \
                   #HEADCROP to indicate the number of bases to be removed at the beginning 
                     HEADCROP:0 \
                   #ILLUMINACLIP to specifies a file of adapter sequences, and the number of mismatches allowed
                     ILLUMINACLIP:/usr/local/programs/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
                    #LEADING AND TRAILING specify the minimum quality for trimming the start and end of reads
                    # SLIDINGWINDOW to indicate the size of the sliding window and  minimum average quality for the bases in that window.
                    # MINLEN to specify the minimum length for a read to be kept
                     LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:36
    done
}
#write success output and error output and run in background
trimAll 1>trimAll.log 2>trimAll.err &
