#! /usr/bin/env bash
# new.sh

# setting the file path
filepath="/scratch/SampleDataFiles/Paired/"

# assigning the left and right reads
left=".R1.paired.fastq"
right=".R2.paired.fastq"

# setting output directory
outDir='quant/'


function align {
   for leftInFile in $filepath*$left
   do
   # obtain sample names
      pathRemoved="${leftInFile/$filepath/}"
      sampleName="${pathRemoved/$left/}"
      echo $sampleName
      # Mapping files
      salmon quant -l IU \
         -1 /scratch/SampleDataFiles/Paired/$sampleName.R1.paired.fastq \
         -2 /scratch/SampleDataFiles/Paired/$sampleName.R2.paired.fastq \
         -i AipIndex \
         --validateMappings \
         -o $outDir$sampleName
   done
}

   align 1>align.log 2>align.err &
