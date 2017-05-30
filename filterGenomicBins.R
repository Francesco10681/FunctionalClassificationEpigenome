#!/usr/bin/Rscript
####################################################################################
####################################################################################
##### Script to select and filter genomic bins according to P-value threshold  #####
####################################################################################
####################################################################################
require(MASS);
args <- commandArgs(trailingOnly = TRUE)
options("scipen"=10, "digits"=4)
inputdata.dir <- args[1];     # This is the directory where sub-directories for each epigenetic mark are placed.
dist.model <- args[2];        # distribution assumed ('Poisson'  or 'NegativeBinomial')
pval.thrs <- args[3];         # Pvalue threshold




#########################
# 0 Input arguments;#####
#########################
#0. `Define Epigenetic tracks
chipseq.marks <- c("CTCF", "Pol2", "DNase-seq" , "H2A.Z",  "H3K4me1",  "H3K4me2",  "H3K4me3",  "H3K27ac",  "H3K27me3", "H3K36me3", "H3K9ac", "H3K9me3", "H3K79me1")





#######################################################################################
# 1. Import tracks and merge into a unique matrix, then export the whole matrix: ######
#######################################################################################

if (file.exists(file.path(inputdata.dir, "normalizedCompleteMatrix.txt") ) ) {

dat <- read.table(file.path(inputdata.dir, "normalizedCompleteMatrix.txt"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(dat) <- c(c("chrom", "start", "end"), chipseq.marks)  }   else  {  # if 'normalizedCompleteMatrix.txt' have not been yet generated...


chipseq.av.marks <- NULL
for (i in c(1:length(chipseq.marks))) {
chipseq.mark <- chipseq.marks[i]
normalizedfile.dir<- file.path( inputdata.dir, chipseq.mark )
if (!(file.exists(normalizedfile.dir))) { cat (paste("directory for track:", chipseq.mark, sep=""), "does not exist!", "\n", sep=" ") 
										next;	} 
chipseq.av.marks <- c(chipseq.av.marks, chipseq.mark)
normalizedfile <- file.path(normalizedfile.dir, "pooledNormalizedSignal.bed")
normalizedsignal <- read.table(normalizedfile, header=FALSE, sep="\t", stringsAsFactors=FALSE)
if (i == 1) { dat <- normalizedsignal[,c(1:3)]
			colnames(dat) <- c("chrom", "start", "end" ) } 
dat <- cbind(dat, normalizedsignal[,4])
colnames(dat)[ncol(dat)] <- chipseq.mark
cat(chipseq.mark, "imported",  "\n", sep=" ")
rm(normalizedsignal)
}
# Export the complete datamatrix:
cat("exporting normalized full matrix..", "\n", sep=" ");
write.table(dat, file=file.path(inputdata.dir, "normalizedCompleteMatrix.txt"), col.names=F, row.names=F, sep="\t", quote=F, dec="." )  }





##################################################################
# 2. Remove completely empty bins  (keep genomic coordinates) ####
##################################################################
cat("remove null bins...", "\n", sep=" ");
rs <- rowSums(as.matrix(dat[,c(4:ncol(dat))]) )
nonempty.bins <- dat[which(rs>0),]
rm(dat);gc();





#########################################
# 3. Convert into Pvalue matrix: ########
#########################################
cat("filter genomic bins using P-value cutoff", pval.thrs, "on a", dist.model, "distribution", "\n", sep=" ");


if ( dist.model == "Poisson" ) { ####  Function to calculate p-values based on the 0.01% of a Poisson-tail ######

# Create a new table with equal dimensions of 'dataf' but reporting the poisson pvalue from the lowertail
pvalmx <- data.frame(bin=c(1:nrow(nonempty.bins)), stringsAsFactors=F)
for (chipseq.mark in chipseq.av.marks) {
cat(chipseq.mark, "\n", sep=" ")
chipseqmark.mean <- mean(nonempty.bins[,chipseq.mark])
pvals <- round(ppois( nonempty.bins[,chipseq.mark], lambda=chipseqmark.mean, lower.tail=FALSE ),8)
pvalmx <- data.frame(pvalmx, pval=pvals, stringsAsFactors=F)
colnames(pvalmx)[ncol(pvalmx)] <- chipseq.mark
}


## Convert the Poisson matrix 'geninfo' in a 0/1 matrix  (cutoff <= 10^-4)
pvalmx <- as.matrix(pvalmx[,c(2:ncol(pvalmx))] )
callmx <- pvalmx
callmx[callmx <= pval.thrs] <- 10
callmx[(callmx > pval.thrs & callmx <= 1)] <- 0
callmx[callmx == 10] <- 1
# Compute sums per row
callmxRS <- rowSums(callmx)
# Filter out bins without at least one significant bin:
fl.xdata <- nonempty.bins[which(callmxRS>0),]

}  else
 
 
 
if ( dist.model == "NegativeBinomial" )  {   ###Function to calculate pvalues based on a Negative-Binomial Distribution ######


# Estimate gamma-shape for each mark
xdata <- as.matrix(nonempty.bins[,c(4:ncol(nonempty.bins))] )
shapes <- NULL
for (j in c(1:ncol(xdata))) {
cat("fit on Gamma-distribution for", chipseq.av.marks[j], "\n", sep=" ")
x <- (xdata[,j] + 0.0001)  # add min. small constant
gammafit <- glm(x~1,family=Gamma)
cat("Calculate maximum-likelihood estimate for the shape gamma parameter") 
gammaestimate  <- gamma.shape(gammafit, verbose=TRUE);
shape <- as.numeric(gammaestimate[1])
shapes <- c(shapes, shape)
}
names(shapes) <- chipseq.av.marks;
rm(xdata);gc();



# Estimate probabilities:
pvalmx <- data.frame(bin=c(1:nrow(nonempty.bins)), stringsAsFactors=F)
for (l in c(1:length(chipseq.av.marks))) {
chipseq.mark <- chipseq.av.marks[l]
cat(chipseq.mark, "\n", sep=" ")
x <- nonempty.bins[,chipseq.mark]
chipseqmark.mean <- mean(x)
# Parametrized Negative-Binomial via "mu" and "size" (gamma shape parameter)
pvals <- pnbinom(q=x ,mu=chipseqmark.mean, size=shapes[l], lower.tail=FALSE)
pvalmx <- data.frame(pvalmx, pval=pvals, stringsAsFactors=F)
colnames(pvalmx)[ncol(pvalmx)] <- chipseq.mark
}


## Convert the Negative Binomial matrix in a 0/1 matrix  (cutoff <= 10^-4)
pvalmx <- as.matrix(pvalmx[,c(2:ncol(pvalmx))] )
callmx <- pvalmx
callmx[callmx <= pval.thrs] <- 10
callmx[(callmx > pval.thrs & callmx <= 1)] <- 0
callmx[callmx == 10] <- 1
# Compute sums per row
callmxRS <- rowSums(callmx)
# Filter out bins without at least one significant bin:
fl.xdata <- nonempty.bins[which(callmxRS>0),]

}
 
 




#########################################
# 4. Export filtered genomic bins: ######
#########################################
# export matrix with bins having one or more significant epigenetic modification
write.table(fl.xdata, file.path(inputdata.dir, paste(dist.model, "normalizedFilteredMatrix.txt",sep=".") ), col.names=F, row.names=F, sep="\t", quote=F, dec="." )

## end script;









