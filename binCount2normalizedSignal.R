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


input.datafile <- args[1];
datasheetsamples.file <- args[2];
mappabilitytrack.file <- args[3];
binnedgenome.file <- args[4];
binsize <- args[5];
totMappedAllExpFile <- args[6];
# Retrieve the number of informative reads mapped for each sample Bi in the group P  (cell-linX/Signal-trackY) 
totMappedAllExp <- scan(totMappedAllExpFile, what="character")
totMappedAllExp <- as.numeric(unlist(strsplit(totMappedAllExp, split=" ")))




####################################
# 0. Import mappability track ######
####################################
cat("import uniqueness mappability track", "\n", sep=" ") 
mappabilitytrack <- read.table(mappabilitytrack.file, header=F, sep="\t", stringsAsFactors=F)
colnames(mappabilitytrack) <- c("chrom" , "chromStart", "chromEnd", "score" )





######################################
# 1. Import datasheet samples file ###
######################################
cat("import", "datasheet samples file", "\n", sep=" ")
datasheetsamples = read.table(datasheetsamples.file, header=F, sep="\t", stringsAsFactors=F, fill=TRUE)
colnames(datasheetsamples) <- c("cell_line" , "data_type", "assay", "sampleID", "alignfile", "map_type" ,"replicate" ,"rawfile" ,"lab")






###############################
# 2. Import mapping files ####  >>> to count the number of mapped reads in each  'processed.bed' (sample) of the group (Bi)
###############################
input.data <- scan(input.datafile, what="character" , sep="\t" )

totMappedSamples = NULL;

for (filepath in input.data) {
cat("import", filepath, "to count nr. of mapped sequencing tags", "\n", sep=" ")
pathsections <- unlist(strsplit(filepath, split="[/]"))
sampleid <- unlist(strsplit(pathsections[length(pathsections)], split="[.]"))[1]
alignfilename <- datasheetsamples[match( sampleid, datasheetsamples$sampleID), "alignfile"]
processedfilename <- gsub("bed", "processed.bed", alignfilename)
processedfilename <- gsub("bam", "processed.bed", processedfilename)
# Reconstruct the whole filepath for the 'processed.bed' file
processedfilepath <- paste( c(  pathsections[2], 
								pathsections[3], 
								pathsections[4], 
								pathsections[5], 
								"data", 
								pathsections[7], 
								pathsections[8], 
								sampleid, 
								processedfilename),
								collapse="/")
processedfilepath <- paste("/",processedfilepath, sep="");
# NB: function 'import.bed' directly imports a bed file in R as a GRanges object!
processedbed <- import.bed(con=processedfilepath, asRangedData=FALSE)
# Save the total number of informative mapped reads in each sample
totMappedSamples = c(totMappedSamples, length(processedbed ))   # vector containing total number of tags in each sample in the group
rm(processedbed);   
 }


# Regenerate the common directory and export txt file with the read totals per sample:
 samplegroup.path <- paste(pathsections[c(1:(length(pathsections)-1))], collapse="/")  # 'samplegroup.path' corresponds to >>> ./genomic_survey/cell.line/signaltrack
 write(totMappedSamples, file=file.path(samplegroup.path, "sample_tag_totals.txt"), sep="\n",  append=FALSE)
 
 
 
 
 
 
 
#######################################
# 3. Import genomic bins file #########
cat("import binned genome file", "\n", sep=" ")
binned.genome = read.table(binnedgenome.file, header=F, sep="\t", stringsAsFactors=F, fill=TRUE)
colnames(binned.genome) <- c("chrom" , "chromStart" ,"chromEnd")
binned.genome <- data.frame(binned.genome, counts=rep(0, nrow(binned.genome)), stringsAsFactors=FALSE)
binned.genome <- data.frame(binned.genome, name=paste(binned.genome$chrom, binned.genome$chromStart, sep="_"), stringsAsFactors=F)
 
 
 
 
 
 
 
 
 

#########################################################################
### 5. Import raw-counts and mappability track into GRanges objects #####
#########################################################################

# 5.1 Convert Mappability track in GRanges object
grmap <- GRanges(   
					seqnames=( Rle(names(table(mappabilitytrack$chrom)[unique(mappabilitytrack$chrom)]), as.numeric(table(mappabilitytrack$chrom)[unique(mappabilitytrack$chrom)]) ) ),
					ranges=IRanges(mappabilitytrack$chromStart, mappabilitytrack$chromEnd),
					strand=Rle('*', nrow(mappabilitytrack)),
					score=mappabilitytrack$score
			 	)

# 5.2 Calculate mappable genome size:
mappabilityIntervalSizes <- width(ranges(grmap));
mappabilityGenomeSize <- sum(as.numeric(mappabilityIntervalSizes));




for (i in c(1:length(input.data))) {

	filepath <- input.data[i]
    cat("import", filepath, "for estimating normalized signal, observed score,expected score ", "\n", sep=" ")
	# 5.3 import single rawCount file (only bins with >= 1 tag)
	
	bincounts <- read.table(file=filepath, header=FALSE, sep="\t", stringsAsFactors=FALSE)
	colnames(bincounts) <- c("chrom", "chromStart", "chromEnd", "counts" )
	bincounts <- data.frame(bincounts, name=paste(bincounts$chrom, bincounts$chromStart, sep="_"), stringsAsFactors=F)
	
	# Add empty bins to the bincounts table (genomic bins with no tags)
	binned.genome2 <- binned.genome
	binned.genome2[match(bincounts$name, binned.genome2$name), "counts"] <- bincounts$counts
	bincounts <- binned.genome2
	rm(binned.genome2)
	
	# 5.4 Convert bincounts in GRanges object!
	### NB: this 'GRanges' command is good for all inputs that are already sorted by chromosome.
	grbins <- GRanges( 
					seqnames=( Rle(names(table(bincounts$chrom)[unique(bincounts$chrom)]), as.numeric(table(bincounts$chrom)[unique(bincounts$chrom)]) ) ),
					ranges=IRanges(bincounts$chromStart, bincounts$chromEnd),
					strand=Rle('*', nrow(bincounts)),
					score=bincounts$counts
					 )
		
		
					 
	# 5.5 Find overlap between bins and uniqueness mappability frames
	cat("examine overlap between genomic bins and mappability regions", "\n", sep=" ")
	bin2mapHits <- findOverlaps( query=grbins, subject=grmap, maxgap = 0L, minoverlap = 1L,ignore.strand = TRUE)
	
	# Determine sizes of overlap regions for each query with any match
	w <- width(ranges(bin2mapHits, ranges(grbins), ranges(grmap))) 

	
	# Build a three-column table; queryIndexes, subjectIndexes, subjectItemSizes
	queryToSubject <- cbind( bin2mapHits@queryHits, bin2mapHits@subjectHits, w)
	colnames(queryToSubject) <- c("queryidx", "subjectidx", "size")
	# Estiamte number of unique-mappability positions within each genomic bin
	queryToSubject <- as.data.frame(queryToSubject, stringsAsFactors=F)
	uqmapp <- as.numeric(tapply(queryToSubject$size, queryToSubject$queryidx, sum) )
	# Estimate the saling factor on sample size
  	ScalingFactorLibrarySize = round( (totMappedSamples[i]/mean(totMappedAllExp)), 3) # totMappedSamples[i] = tot.mapped tags in sample [i]
	# Estimate the normalized expected fragment count in each bin
	expCounts= round( (uqmapp * (totMappedSamples[i] /mappabilityGenomeSize ) )* ScalingFactorLibrarySize, 3)
	# Estimate the normalized observed count in each bin
  	obsCounts= round ( (score(grbins)[unique(bin2mapHits@queryHits)] ) , 3)
  	NormalizedSignals <- round((obsCounts/expCounts), 3)
  	
  	
  	
  	
  	# 5.6 Export: 
  	filepathexport <- gsub( ".rawCounts.bed", ".normalizedSignal.bed", filepath)

  	cat("exporting normalizedSignal bed file for", filepath, "\n", sep=" ")
  	write.table( data.frame( chrom=as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
  							 chromStart=start(ranges(grbins))[unique(bin2mapHits@queryHits)],
  							 chromEnd=end(ranges(grbins))[unique(bin2mapHits@queryHits)],
  							 #name=paste(as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)], 
  							 #			start(ranges(grbins))[unique(bin2mapHits@queryHits)], sep="_"),
  							 score=NormalizedSignals,
  							 #strand=rep("*", length(unique(bin2mapHits@queryHits))),
  							 stringsAsFactors=FALSE),
			file=filepathexport, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);
  	
  	
  	
  	
  	
  	# 5.7 Export only expCount values
  	filepathexportExp <- gsub( ".rawCounts.bed", ".expectedScore.bg", filepath)
  	cat("exporting expected signal bed file for", filepath, "\n", sep=" ")
  	write.table( data.frame( chrom=as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
  							 chromStart=start(ranges(grbins))[unique(bin2mapHits@queryHits)],
  							 chromEnd=end(ranges(grbins))[unique(bin2mapHits@queryHits)],
  							 score=expCounts,
  							 stringsAsFactors=FALSE),
			file=filepathexportExp, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);
  	
  	
  	
  	
  	
  	# 5.8 Export only obsCount values
  	filepathexportObs <- gsub( ".rawCounts.bed", ".observedScore.bg", filepath)
  	cat("exporting observed signal bed file for", filepath, "\n", sep=" ")
  	write.table( data.frame( chrom=as.character(seqnames(grbins))[unique(bin2mapHits@queryHits)],
  							 chromStart=start(ranges(grbins))[unique(bin2mapHits@queryHits)],
  							 chromEnd=end(ranges(grbins))[unique(bin2mapHits@queryHits)],
  							 score=obsCounts,
  							 stringsAsFactors=FALSE),
			file=filepathexportObs, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE);
  	
	};


### End of script ####
