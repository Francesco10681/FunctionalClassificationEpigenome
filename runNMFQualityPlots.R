#!/usr/bin/Rscript

#####################################################################################################
#####################################################################################################
####### Function to run NMF and generate Quality plot for increasing k (Factorization Ranks) ########
#####################################################################################################
#####################################################################################################

library(fpc);
library(RColorBrewer);
library(gplots);
library(ggplot2);
library(foreach);
library(doParallel);
library(NMF);


#################################################################
###### Function to generate multiplots in the same pdf ##########
#################################################################




### Function to perform multiplots in the same PDF ###########
##############################################################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}






### Define function to report quality plots for a series of factorization ranks ###########
###########################################################################################
options(scipen=0.1)

createQualityPlots <- function(estimname,  mx, ranks, transform.method, distribtype, plots.dir) {    #name of the resulting object;; #original matrix;; #vector of factorization ranks tested;;
# get input data:
estimates <- get(estimname)


cophcoef.v <- NULL
rss.v <- NULL
basis.sparseness.v <- NULL
coef.sparseness.v <- NULL


for (k in c(1:length(ranks))) {

rank.factor <- ranks[k];
cat("extract", "nmf results", "for factorization rank", rank.factor, "\n", sep=" ")
singlerank.results <- estimates[[k]]
# Extract the three elements from the list:
nmfobject <- unlist(singlerank.results[[1]]);
# get cophenetic coefficient
cophcoef <- as.numeric(unlist(singlerank.results[[3]]));
# get sparseness of H
coef.sparseness <- as.numeric( sparseness(nmfobject))[2];
# get sparseness of W
basis.sparseness <- as.numeric( sparseness(nmfobject))[1];
# get residual Sum of squares
rss <- rss(nmfobject, mx);

# Collect quality measurement:
cophcoef.v <- c(cophcoef.v, cophcoef)
rss.v <- c(rss.v, rss)
coef.sparseness.v <- c(coef.sparseness.v, coef.sparseness)
basis.sparseness.v <- c(basis.sparseness.v, basis.sparseness)  }


# Generate Quality Plots PDF
pdf(file=file.path(plots.dir, paste( transform.method, distribtype, "quality.test.NMF.pdf", sep=".")), width=15, height=4)

# Build the main dataframe
dat <- data.frame(ranks=ranks,
                 coph=cophcoef.v,
                 rss=rss.v,
                 stringsAsFactors=F)
                  
# Build the 'sparseness' dataframe           
sparseness <- data.frame( ranks=rep(ranks,times=2),
						spar=c(coef.sparseness.v, basis.sparseness.v),
						class=rep(c('coef', 'basis'), each=length(ranks)),
						stringsAsFactors=F)                  
						     
# Cophenetic correlation coefficient
plot1 <- ggplot(dat, aes(ranks, coph, color="green3")) + geom_line(color="yellowgreen", size=1.1) + geom_point(color="darkgreen", size=3) + ggtitle("Cophenetic coeff.");
# Sparseness
plot2 <- ggplot(sparseness, aes(ranks, spar, group=class, color=factor(class)) ) + geom_line(size=1.1) + geom_point(size=3) + ggtitle("Sparseness of W and H matrixes");
# RSS
plot3 <- ggplot(dat, aes(ranks, rss) ) + geom_line(color="red2", size=1.1) + geom_point(color="red4", size=3) +  ggtitle("Residual sum of squares(RSS)");
# Call the 'multiplot' function;
multiplot(plot1,plot2,plot3,cols=3);

dev.off();   }










##############################
###  Import arguments ########
##############################
chipseq.marks <- c("CTCF", "Pol2", "DNase-seq" , "H2A.Z",  "H3K4me1",  "H3K4me2",  "H3K4me3",  "H3K27ac",  "H3K27me3", "H3K36me3", "H3K9ac", "H3K9me3", "H3K79me1")


args <- commandArgs(trailingOnly = TRUE)
options("scipen"=1000)
input.datamatrix <- args[1];     # The path to the input data-matrix (cols 1-3: chrom/start/end;; cols 4-n: epigenetic normalized tracks)
transform.method <- args[2];     # Type of data transformation
ranks.range <- args[3];			 # Range of factorization ranks to test
nruns <- as.numeric(args[4]);			# Number of runs to perform for the NMF algorithm
output.dir <- args[5];           # Define the output directory
#conversion <- as.numeric(args[6]);			# Convert ( or not ) the input matrix into a data-matrix containing only coverage signals
fact.rank <- as.numeric(args[6]);    # Define the factorization rank for which NMF results have to be extracted 


###############################
###  Import data-matrix #######
###############################
mx <- read.table(file=input.datamatrix, header=FALSE, stringsAsFactors=FALSE, sep="\t")
genomic.bins <- mx[,c(1:3)];
colnames(genomic.bins) <- c("chrom", "chromStart", "chromEnd")
mx <- as.matrix(mx[c(4:ncol(mx))]) #} else 
#if (conversion == 0 )  { mx <- as.matrix(mx) }
parts <- unlist(strsplit(input.datamatrix, split="/"))
filename <- parts[length(parts)]
distribtype <- unlist(strsplit(filename, split="[.]"))[1]




#############################
##### Data transformation ###
#############################
cat("Data transformation required:", transform.method, "\n", sep=" ")

if (transform.method == "raw") { colnames(mx) <- chipseq.marks } else

# Apply sigmoid function transformation of the data:
if (transform.method == "sigmoid" )  { 
						for (j in c(1:ncol(mx) ) ) {
						cat("column", j, "\n", sep=" ")
						xj <- mx[,j]
						perc95th <- as.numeric(quantile(xj, probs=.95))
					    C = 2/(1+2.72^(-2*xj/perc95th)) - 1 
					    C[is.na(C)] <- 0
						mx[,j] <- C  }
						rm(C);rm(xj);gc();
						colnames(mx) <- chipseq.marks  }  else
						
if (transform.method == "invhyperbolic") {
						  mx = log2(mx + sqrt(mx^2+1)) 
						  colnames(mx) <- chipseq.marks }
											




########################
##### NMF package; #####
########################

##1. estimate best ranks
#nmfranks <- nmfEstimateRank( xdata, range=c(3:6), method='brunet' , verbose=T , .options='p4')
# get best factorization rank
# 2. run NMF
#nmfres <- nmf( xdata,  rank=3:6, method='brunet', nrun=20, .options='pv' )  # launch in parallel on  4 cores 
# 3. get cluster memberships




## Random test (repeat the randomization process until no full-zero rows are found):
#flag <- "white"
#while( flag == "white") {
#rand.matrix <- randomize(nbmatrix)
#rs <- rowSums(rand.matrix)
#if (!(any(rs==0) )  ) { flag <- "red"} }  # or poimatrix
#mx <- rand.matrix
# Then sigmoid transformation  # or NOT...


###############################################################################
## Run NMF analysis!! #########################################################
## nmf + lapply: set of n-runs (foreach k) will be launched sequentially! #####
###############################################################################
ranks.range <- unlist(strsplit(ranks.range, split="[,]"))
ranks.range <- as.numeric(ranks.range)
ranks <- c(min(ranks.range):max(ranks.range))

estim <- lapply(ranks, function(r){
  fit <- nmf( mx, rank=r, nrun=nruns, seed=123456, method='brunet', .options="vp")
  list(fit = fit, consensus = consensus(fit), coph = cophcor(fit))
})
names(estim) <- paste('rank', ranks)   

##### Generate the quality Plots:
createQualityPlots("estim", mx, ranks, transform.method, distribtype, output.dir);






#######################################################################################################
# Extract NMF results from a specific factorization rank (OPTIONAL): { only if 'fact.rank'>0 } ########
#######################################################################################################

if (fact.rank>0) {

customdistFUN <- function(x) {
				as.dist(1-cor(t(x), method="pearson"))  }
customhclustFUN <- function (x) {
				hclust(x, method="average") }

# A)  HEATMAP1: H-matrix
pidx <- match(fact.rank, ranks)
nmfres <- (estim[[pidx]])   # take the element of 'estim' Rlist related to the specific factorization rank.
binmap <- as.numeric(predict(nmfres[[1]], what="features"))
genomic.bins <- data.frame(genomic.bins,  name=paste(genomic.bins$chrom, paste(genomic.bins$chromStart, genomic.bins$chromEnd, sep="_"), sep=":"), score=binmap, strand=rep('*', nrow(genomic.bins)), stringsAsFactors=FALSE )
h <- coef(unlist(nmfres[[1]]))
# Plot the coefficients heatmap
dist.max <- unique(max(h))
dist.min <- unique(min(h))
breaks <- seq(0.001,dist.max, 0.01) 
colors.palette.6 <- colorRampPalette(c("dodgerblue4", "gold1", "red"))(length(breaks)-1)
pdf(file=file.path(output.dir,paste( transform.method, distribtype, paste('fr',fact.rank,sep=""),  "coeffMatrix", "pdf", sep=".")),height=6, width=10)
heatmap.2(h, Rowv=FALSE, distfun=customdistFUN, hclustfun=customhclustFUN, scale="none", col=colors.palette.6, dendrogram="col", hline=FALSE, vline=FALSE, margins=c(10,14), breaks=breaks,
sepwidth=c(0.03,0.03), labRow=paste("profile", c(1:nrow(h)), sep=" "), labCol=colnames(h), cexRow=1.6, cexCol=1.9, trace='none', key=TRUE, keysize=1.1, cex.main=.4, ,density.info=c("histogram","density","none"),
main="Coefficients matrix of chromatin profiles" )
dev.off() 



# B) HETMAP2: Mean Chromatin mark Signal per Profile
cclasses <- c(1:max(ranks))
signalstate.matrix <- matrix( rep(0, (length(cclasses) * length(chipseq.marks) )), nrow=length(cclasses), ncol=length(chipseq.marks) )
rownames(signalstate.matrix) <- cclasses
colnames(signalstate.matrix) <- chipseq.marks
for (cclass in cclasses) {
cat(cclass, "\n", sep=" ")
code.means <- colMeans(mx[which(genomic.bins$score==cclass),])
signalstate.matrix[cclass, ] <- code.means
}
breaks <- seq(0, 0.6, 0.01) 
colors.palette.6 <- colorRampPalette(c("dodgerblue4", "gold1", "red"))(length(breaks)-1)
pdf(file=file.path(output.dir, paste( transform.method, distribtype, paste('fr',fact.rank,sep=""),  "meanChromatinSignal", "pdf", sep=".")),height=6, width=10)
heatmap.2(signalstate.matrix[,chipseq.marks], Rowv=FALSE, Colv=FALSE, distfun=customdistFUN, hclustfun=customhclustFUN, scale="col", col=colors.palette.6, dendrogram="none", hline=FALSE, vline=FALSE, margins=c(10,14), #breaks=breaks,
sepwidth=c(0.03,0.03), labRow=paste('profile', rownames(signalstate.matrix), sep=" "), labCol=colnames(signalstate.matrix[,chipseq.marks]), cexRow=1.6, cexCol=1.9, trace='none', key=TRUE, keysize=1, cex.main=.3, ,density.info=c("histogram","density","none"),
main="Average signal in each profile" )
dev.off() 


## C) H-COEFFICIENT MATRX
write.table(h, file=file.path(output.dir, paste( transform.method, distribtype, paste('fr',fact.rank,sep=""),  "coeffMatrix", "txt", sep=".") ),
quote=FALSE, row.names=FALSE, sep="\t")


## D) PROFILES/BINS ASSIGNMENT in Tabular format:
write.table(genomic.bins, file=file.path(output.dir, paste(transform.method, distribtype, paste('fr',fact.rank,sep=""),  "profileAssignment", "bed", sep=".") ),
row.names=FALSE, quote=FALSE, dec=".", sep="\t" )        }

# Save 'estim' object:
save( list=c("estim"),file=file.path(output.dir, paste( transform.method, distribtype, "nmf.distributions.RData", sep=".") ) )





##########
## End ###
##########
