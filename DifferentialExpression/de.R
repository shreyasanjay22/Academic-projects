#!/usr/bin/env Rscript
# de.R

library(tximport)
library(readr)
library(DESeq2)
tx2gene <- read.csv("tx2gene.csv")
head(tx2gene)

# setting the file path
samples <- read.csv("/scratch/SampleDataFiles/Samples.csv", header=TRUE)
head(samples)

# importing the counts with tximport
files <- file.path("quant", samples$Sample, "quant.sf")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ Menthol + Vibrio)

# performing DESeq
dds$Vibrio <- relevel(dds$Vibrio, ref = "Control")
dds$Menthol <- relevel(dds$Menthol, ref = "Control")
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds)


padj <- .05
minLog2FoldChange <- .5
dfAll <- data.frame()
# Get all DE results except Intercept, and "flatten" into a single file.
for (result in resultsNames(dds)){
     if(result != 'Intercept'){
        res <- results(dds, alpha=.05, name=result)
        dfRes <- as.data.frame(res)
        dfRes <- subset(subset(dfRes, select=c(log2FoldChange, padj)))
        dfRes$Factor <- result
        dfAll <- rbind(dfAll, dfRes)
                    }
}


write.csv(dfAll, file="dfAll.csv")
head(dfAll)

path <-read.table("/scratch/SampleDataFiles/Annotation/path.txt",sep="\t",header=FALSE)
colnames(path) <- c("ko","pathway")

pathway <-read.table("/home/shirodkar.sh/BINF6308/KEGG_GO/ko", sep="\t",header=FALSE)
colnames(pathway)<-c("pathway","path")

mergeko=read.csv("dfAll.csv")
colnames(mergeko)=c("ko","log","padj","factor")

finalMerge1<-merge(path,pathway)
finalMerge<-merge(finalMerge1,mergeko)

write.table(finalMerge,file="deAnnotated.csv")





