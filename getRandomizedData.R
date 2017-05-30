#!/usr/bin/Rscript

#####################################################################################
#####################################################################################
####### Function to prepare randomized inout dataset for NMF (Factorization Ranks) ##
#####################################################################################
#####################################################################################

library(fpc);
library(NMF);

args <- commandArgs(trailingOnly = TRUE)
options("scipen"=1000)
input.datamatrix <- args[1];     
output.dir <- args[2];

## Import the real data matrix and split between genomic coordinates and numeric values
mx <- read.table(file=input.datamatrix, header=FALSE, stringsAsFactors=FALSE, sep="\t")
genomic.bins <- mx[,c(1:3)];
mx <- as.matrix(mx[,c(4:ncol(mx))] )
#colnames(genomic.bins) <- c("chrom", "chromStart", "chromEnd")
rmx <- randomize(mx);


# Define the randomized data matrix file name:
parts <- unlist(strsplit(input.datamatrix, split="/"))
filename <- parts[length(parts)]
outputfilename <- gsub("Filtered", "Random", filename);
if (!(dir.exists(file.path(output.dir, "random") ) ) ) { dir.create(file.path(output.dir,"random")) }
random.matrix <- data.frame(genomic.bins, rmx, stringsAsFactors=FALSE)

# Export output file
output.dir <- file.path(output.dir, "random")
write.table(random.matrix, file=file.path(output.dir, outputfilename), col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, dec=".")

## End script



