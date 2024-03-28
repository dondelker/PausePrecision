# Gets maximum locations of PRO-seq 3' reads and 200 bp matrix

library(rtracklayer)
library(Rsamtools)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

source("bigWig.R")
source("coverage.R")

# Read in GRanges object and create vector of gene names
genes <- GRanges(read.csv("NELFE_Pol2_1kb.csv"))
names <- read.csv("NELFE_Pol2_1kb.csv") [ ,4]

# Read in stranded bigWig files
bw.file = "G:/NS20041-Watts/hg38_alignments/bigWigs/combinedBigWig/Bcell_libs1234"

# Get read coverage values and save data matrix with gene names 
read.depth <- bw.coverage.stranded(bw.file, genes, flip.minus=TRUE)
read.mat <- rle.to.matrix(read.depth)
read.mat <- as.data.frame(read.mat, row.names = names)

# Get max pause site and value
max.loc <- colnames(read.mat)[apply(read.mat, 1, which.max)]
read.mat$max.loc <- gsub("V", "", max.loc)
read.mat$max = apply(read.mat[1:1001], 1, max, na.rm = TRUE)

# Combine pause read stats and matrix into single dataframe
metadata <- as.data.frame(genes)
anno.matrix <- cbind(read.mat[202:203], metadata) 

# Write to file
write.csv(anno.matrix, "Bcell_pauseMax_NELFE_Pol2_regions.csv")

# Plot read stats
plot(read.mat$max.loc, read.mat$max)

