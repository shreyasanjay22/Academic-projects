#! /usr/bin/env bash
# new.sh


filepath="/scratch/SampleDataFiles/Paired/"

left=".R1.paired.fastq"
right=".R2.paired.fastq"

outDir='quant/'


function align {
   for leftInFile in $filepath*$left
   do
      pathRemoved="${leftInFile/$filepath/}"
      sampleName="${pathRemoved/$left/}"
      echo $sampleName
      salmon quant -l IU \
         -1 /scratch/SampleDataFiles/Paired/$sampleName.R1.paired.fastq \
         -2 /scratch/SampleDataFiles/Paired/$sampleName.R2.paired.fastq \
         -i AipIndex \
         --validateMappings \
         -o $outDir$sampleName
   done
}

   align 1>align.log 2>align.err &
