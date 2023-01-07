# Transcriptome Assembly
In this module, the reads are used to perform reference-guided transcriptome assembly and de-novo transcriptome assembly.

# Methods

## Trinity
 Both methods use Trinity for assembling the transcriptome. Trinity consists of 3 consecutive module: Inchworm, Chrysalis and Butterfly.
First the Inchworm module is used to construct contigs. Overlapping k-mers are extracted from RNA-Seq reads and examined. The transcript contigs are generated based on 9k-10-mer overlaps.
The Chrysalis clusters related Inchworm contigs. The structural complexity of clustered Inchworm is encoded by building a de Bruijn graph for each cluster and making partitios amongst the cluster. 
The Butterfly module process individual graphs in parallel and reconstructs transcript sequences in manner to reflect the original cDNA molecule.

## Reference-guided assembly of transcriptome
A reference-guided assembly uses alignments of reads to the genome that are produced by specialised alignment tool, In this case GSNAP to identify reads that represent potential transcripts. An assembly is built for these alignment.

The assembly is done by using Trinity. Trinity requires single sorted bam files as input for this purpose merging of the bam file is performed before the assembly.Samtools is used to merge the files.
The assembly is carried out using the merged file. First, the nice -n19 parameter sets lower priority to this program than the critical function. Next, the path of the Trinity assembler is set. -genome_guided_bam is used to specify genome-guided assembly. -genome_guided_max_intron indicates the maximum separation distance that trinity will allow for segments of a transcript.-max_memory indicates the maximum memory to be used for this program by Trinity. The success and errors are directed to trinity.log and trinity.err files and the program is run in the background.
The statistics of the assembly is checked by running TrinityStats.pl with the assembled transcriptome as input. The output here will show the size of assembled contigs. 

## De-novo assembly:
 In de-novo assembly, the assembly is created without the aid of a reference genome.

A script to save comma-separated lists in two variables is performed as Trinity requires two comma-separated lists of files for the left and right reads. This list is used to perform the assembly. An output directory icalled trinity_de-novo is specified to avoid overwriting with the gemone-guided assembly.The comma separated variables are passed. --max_memory specifies the maximum memory to used by Trinity in this case 4 CPU. The success and errors are written in trinity_dn.log and trinity_dn.err and run in the background.
 

## Citations:
Haas BJ, Papanicolaou A, Yassour M, et al. De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nat Protoc. 2013;8(8):1494–1512. doi:10.1038/nprot.2013.084
