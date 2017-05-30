#!/usr/bin/Rscript

###############################################################################################################################
### Preliminary analysis for the Data-Integration (ENCODE/Roadmap Epigenomics) project: #######################################
###############################################################################################################################
### statistical survey on the overlap between different ChIP-seq data and genes from the GenCode database   @Francesco Gandolfi
###############################################################################################################################


#NB (Task); Convert processed BED files into Normalized (TPM or RPKM) Coverage signal tracks.
args <- commandArgs(trailingOnly = TRUE)
options("scipen"=10, "digits"=4)
require(GenomicRanges);
require(IRanges);
library(rtracklayer);





####################
# 0. Arguments #####
####################
sum.expected.file <- args[1];
sum.observed.file <- args[2];
und <- unlist(strsplit(sum.expected.file, split="/"))
und <- und[c(2:(length(und)-1))]
output.dir <- paste("/", paste(und, collapse="/"), sep="")







#############################
# 2. Import Input files #####
#############################
imported.sumScores = read.table(sum.expected.file, header=F, sep="\t", stringsAsFactors=F, fill=TRUE)
summarizedExScores <- rowSums(imported.sumScores[,c(4:ncol(imported.sumScores))])
imported.sumScores = read.table(sum.observed.file, header=F, sep="\t", stringsAsFactors=F, fill=TRUE)
summarizedObScores <- rowSums(imported.sumScores[,c(4:ncol(imported.sumScores))])
# 1.2 calculate the normalizedCoverageSignal
sumNormSignals <- round((summarizedObScores/summarizedExScores),3)
# Replace Infinite values by observed normalized signal
sumNormSignals[which(is.infinite(sumNormSignals))] <- summarizedObScores[which(is.infinite(sumNormSignals))] 
# Replace NaN values with 0
sumNormSignals[which(is.nan(sumNormSignals))] <- rep(0,sum(is.nan(sumNormSignals)) );



#############################################
# 2. Export SummarizedNormalizedSignals #####
#############################################
filepathexport <- file.path(output.dir, "pooledNormalizedSignal.bed")
cat("exporting pooled normalizedSignal bed file", "\n", sep=" ")
write.table( data.frame( chrom=imported.sumScores[,1],
  					     chromStart=imported.sumScores[,2],
  						 chromEnd=imported.sumScores[,3],
  						 #name=paste( imported.sumScores[,1], imported.sumScores[,2], sep="_"),
  						 score=sumNormSignals,
  						 #strand=rep("*", nrow(imported.sumScores)),
  						 stringsAsFactors=FALSE),
						 file=filepathexport, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);
  	


### End of script ####
