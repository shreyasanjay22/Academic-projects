#!/usr/bin/env Rscript
# mergeAll.R


# Load BLAST results as a table using tab (\t) as separator.
# There is no header with column names, so set header=FALSE
blast <- read.table("../BLAST/alignPredicted.txt", sep="\t", header=FALSE)
# Set column names to match fields selected in BLAST
colnames(blast) <- c("qseqid", "sacc", "qlen", "slen", "length", "nident", "pident", "evalue", "stitle")
# Calculate the coverage as nident/slen

blast$cov <- blast$nident/blast$slen
# Select only blast rows where coverage is greater than .9
blast <- subset(blast, cov > .9, select=-c(stitle))
# Read kegg.txt, produced by get_kegg_ids, as a table
kegg <- read.table("kegg.txt", sep="\t",fill = TRUE , header=FALSE)

# Set the column names for kegg
colnames(kegg) <- c("sacc", "kegg")
# Remove the up: prefix from the accession number so it will match the BLAST
# subject accession (sacc)

kegg$sacc <- gsub("up:", "", kegg$sacc)
blast_kegg <- merge(blast, kegg)
head(blast_kegg)

go <- read.csv("sp_go.txt", sep="\t", header=FALSE)
head(go)

#Read kegg.txt as a table
kegg <- read.table("kegg.txt", sep="\t",fill = TRUE, header=FALSE)
#set column names for kegg
colnames(kegg) <- c("sacc", "kegg")

#read ko.txt as a table
kotxt <- read.table("ko.txt", sep="\t",fill = TRUE, header=FALSE)
# Set the column names for kotxt
colnames(kotxt) <- c("kegg", "KO")

#read kegg_path.txt as table
Keggpath <- read.table("kegg_path.txt", sep="\t",fill = TRUE, header=FALSE)
# Set the column names for Keggpath
colnames(Keggpath) <- c("KO", "pathID")

#read ko as a table
Path<- read.table("ko", sep="\t",fill = TRUE, header=FALSE)
# Set the column names for Path
colnames(Path) <- c("pathID", "Pathdescription")

#merge all the tables
Merge1<- merge(kegg ,kotxt)
Merge2<- merge(Merge1, Keggpath)
Merge3<- merge(Merge2, Path)

#show the first few lines of the results
head(Merge3)



