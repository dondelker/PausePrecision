library(ggplot2)
library(ggseqlogo)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)


genome = BSgenome.Hsapiens.UCSC.hg38

# Load precise pause coordinates
pauses = import("pause_quartile1.bed")

# Get 21nt sequences centered at max pause
sequences = getSeq(genome, resize(pauses, 21, fix="center"))

# Improve the names in the file
names(sequences) = paste0(as.character(pauses), " ", pauses$name)
export(sequences, "sequence_at_pause.fa", "fasta")
sequences.as.strings = as.character(sequences)

# Make ggplot of sequence logo
ggplot() + geom_logo(sequences.as.strings) + ylim(0,1.6) + theme_logo()


