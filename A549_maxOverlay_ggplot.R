# Gets read coverage and percent G for PROseq data

source("bigWig.R")
source("granges.R")

# Load libraries
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPpeakAnno)
library(stringr)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(reshape2)

# Add genome information
genome.hg38.ensembl = BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(genome.hg38.ensembl) = "Ensembl"
data("TSS.human.GRCh38")

# Read in Granges list of CUT&Tag peak regions
CTpeaks = GRanges(read.csv("A549_top100_GRanges.csv"))

# Annotate peak regions
annoPeaks = annotatePeakInBatch(CTpeaks, AnnotationData = TSS.human.GRCh38)

# Retrieve DNA sequence for each peak region
annotate.gene.loc = function(annoPeaks) {
  seqlevelsStyle(annoPeaks) = "Ensembl"
  annoPeaks$seq = as.character(getSeq(genome.hg38.ensembl, annoPeaks))
}
annoPeaks = sort(annoPeaks)
add_sequence = annotate.gene.loc(annoPeaks)

# Adds "chr" to chromosome name and provides location of bigWig files
seqlevelsStyle(annoPeaks) = "UCSC"
g = annoPeaks
bw.file = "C:/Users/delkerda/Desktop/G4_sequences/A549_PROseq/A549_sorted"

# Computes average read depth for each position
read.depth = as.data.frame(bw.coverage.stranded(bw.file, g, flip.minus=TRUE))
read.stats = t(rowMeans(read.depth))

# Computes G percentage of strings of DNA
gc.percentage = function(s)
  str_count(s, "[G]") / str_count(s, "[ACGT]")

# Computes G4 stat function over 10 bp windows
gc.windowed.per.gene = function(s, f, window.size=10) {
  
  # we assume these are all the same length
  n = nchar(s[1])
  
  r = NULL
  
  for(i in 1:(n + 1 - window.size)) {
    if (i %% 100 == 0)
      cat(i, "")
    y = f(substr(s, i, i + window.size - 1))
    r = rbind(r, y)
  }
  
  as.matrix(r)
}

# Calculate %G by windowing
read.seq = as.data.frame(t(rbind(read.depth,add_sequence)))
percentG_matrix <- gc.windowed.per.gene(read.seq$V501, f = gc.percentage)

# Replace NaN and NA values with 0
percentG_matrix[is.nan(percentG_matrix)]<-0
percentG_matrix[is.na(percentG_matrix)]<-0

# Calculate means and stdev for all coordinate positions
A549_pctG <- as.data.frame(rowMeans(percentG_matrix))
A549_pctG$std <- rowSds(percentG_matrix)

# Adjust read coverage to fit %G output
read.cov = as.data.frame((read.seq[,5:495]))
new.cov <- as.data.frame(lapply(read.cov, as.numeric))
read.mean <- as.data.frame(colMeans(new.cov))

# Combine read coverage and percent G stats
mydata <- cbind(A549_pctG,read.mean)

# Create columns for plot coordinates
names(mydata) = c("percentG", "std", "reads")
mydata$coordinates = seq(nrow(mydata))
write.table(mydata,"A549_pause_percentG_stats.txt", sep="\t", row.names = FALSE)

# Plot PROseq reads 
ggplot(mydata, aes(x=coordinates, y=reads, colour=percentG, fill=percentG)) + 
  geom_col(size=2) + geom_vline(xintercept=246, colour="red", linetype="longdash") +
  theme(panel.background=element_blank(), axis.line=element_line(color="black", 
  linetype=1)) + xlab("pause max") + ylab("RPM")




