# Center pause max and add +/- 500 bases to create new bed file

# Read in max.loc data
myregions <- read.delim("eRNA_geneBody_max.txt", header = TRUE, row.names = 1)

# Create new column with pause max center coordinate
myregions$center <- myregions$start + myregions$max.loc

# Add 250 bases on either side of the pause max
myregions$newstart = myregions$center - 500
myregions$newend = myregions$center + 500

# Write new coordinate data to file 
final <- myregions[,c(3,10,11,2,7)]
write.table(final, "eRNA_pauseMax_centered.txt", row.names = TRUE, sep = "\t")