# Short Read Alignment
This directory is created to perform basic steps for running GSNAP for Aiptasia reads. This includes trimmig the Aiptasia  reads, building a genomic GMAP database, aligning the read, and sort and index the reads.

# METHODS
## Quality trimming of the reads:
The input data for trimming was collected from /scratch/AiptasiaMiSeq/fastq. The trimming of read was carried out with Trimmomatic. Trimmomatic includes a variety of processing steps for read trimming and filtering, but the main algorithmic innovations are related to identification of adapter sequences and quality filtering. In trimmomatic, we have used the Sliding window quality filtering approach in which the scanning of work is done from 5' end and removed from 3'end of the read when the average quality of the bases drop below a certain threshold.
### Description of trimming tasks of Trimmomatic:
The shell command nice - 19 tells the server to give lower priority for that function than the critical function
java jar command runs Java jar files as Trimmomatic is a java program
/usr/local/programs/Trimmomatic-0.36/trimmomatic-0.36.jar in the program is the location of the java jar
PE indicates the paired end reads and -threads is set to 1 to indicate that 1 server thread is used for the job, -phred33 indicates quality code method.
HEADCROP is set to 0 which indicates the number of bases to be removed from the beginning regardless of quality score.
ILLUMINACLIP indicates the file of adapter sequence and number of mismatches allowed
LEADING an TRAILING are set to 20 which specifies the quality of trimming at the beginning and the end.
SLIDINGWINDOW is set to 4:30 which indicates the size of the sliding window and the minimum average quality for the bases in that window.
MINLEN is set to 36 which specifies the minimum length for a read to be kept.

### Steps followed:
First the path of the file,suffix of left and right reads,output  directory, are initialized.The output directory are created. The function trimAll is set and ran in the loop through the left-read fastq files in the input file.In the function, the path of the file is removed and assigned to pathRemoved. The Further the left suffix is removed from these files and assigned to $sampleName and printed.Next the trimmomatic tasks mentioned above are run  for the left and right read files and stores the output in  PAired and Unpaired files. The last line indicates the file to write the success output(1>) and error output(2>).

## Aligning of the reads:
The trimmed samples are aligned against the GMAP database.gmap_build is used to create a GMAP database from the Aiptasia genome.GSNAP uses this databse to perform alignment of RNA-Seq reads.SAM stands for Sequence Alignment/Map format. In this program the output of alignment in SAM format.It consists of an optional header section and an alignment section.Each alignment section has 11 mandatory fields essential for alignment information.


### Steps followed:
Similar to the steps in trimming, the path of the file,suffix of reads, format of the file and the output files are initialised. The output directory is created. The function alignReads is run for a loop through all the left reads files in the input file. Here the input file is the quality trimmed reads stored in Paired directory. In the function similar to trimAll the path of the file and suffix is removed and stored in $sampleName. Next the command -A indicates gsnap to create sam alignment program. The working directory is set for the program.The name of the database used for aligning is set using -d command.The next parameters indicates the files to be used to create sam alignment outputs. The information from  sam alignment are directed to a new directory in sam format.The loop is completed and the success and errors are written in .log .err files.

## Sorting and Indexing of the Aligned reads:
The SAM files are sorted and converted to BAM format. samtools utilities  are used for this purpose. 
### Steps followed:
Similar to the trim and align program, the input path is initialised to contain the files from  sam directory. The output format is initialised to .sorted.bam. The output directory for this file is bam. The function sortAll is run for a loop through all the sam files in the sam directory. The path of the file and suffix is removed. From the samtool sort command is used to sort the files from sam directory $filepath$sampleNAme$left and the output is redirected to the bam directory. .log and .err files are written for success and error outputs.
Similar steps are followed for indexing where samtool index command is used for indexing and the results are in the form of bai.files in the sam directory.

## Citations:
Anthony M. Bolger, Marc Lohse, Bjoern Usadel, Trimmomatic: a flexible trimmer for Illumina sequence data, Bioinformatics, Volume 30, Issue 15, 1 August 2014, Pages 2114–2120.
