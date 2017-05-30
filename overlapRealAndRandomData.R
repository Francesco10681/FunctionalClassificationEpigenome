#!/usr/bin/Rscript

###################################################################################################################
###################################################################################################################
####### Function to overlap real and random results in NMF Quality plot for increasing factorization rank  ########
###################################################################################################################
####################################################################################################################
library(RColorBrewer);
library(NMF);
library(gplots);
library(ggplot2);


## 0) Multiplot function

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





####################################################
### 1) Define arguments and recall Rlibraries ######
####################################################

options("scipen"=0.1)

args <- commandArgs(trailingOnly = TRUE)
sample.matrix.file <- args[1];     # Define real original matrix  (Poisson or NegativeBinomial)
sample.nmfresults.file <- args[2];     # Define nmf results for real data
filename <- unlist(strsplit(sample.nmfresults.file, split="/")) [length(unlist(strsplit(sample.nmfresults.file, split="/")))]
transform.method <- unlist(strsplit(filename, split="[.]"))[1];
statdist <- unlist(strsplit(filename, split="[.]"))[2];
random.matrix.file <- args[3];			 # Define random original matrix
random.nmfresults.file <- args[4];       # Define nmf result file for random data
output.dir <- args[5];			# Define the output directory







chipseq.marks <- c("CTCF", "Pol2", "DNase-seq" , "H2A.Z",  "H3K4me1",  "H3K4me2",  "H3K4me3",  "H3K27ac",  "H3K27me3", "H3K36me3", "H3K9ac", "H3K9me3", "H3K79me1")


### Then, we have to import both the random and the real data-matrix in the same R-workspace and apply data-transformation on both if required.
## NB!!! Remember that compared to the random matrix, real matrix still have genomic bin coordinates in the first three columns!!







######################################
## 2) Import original matrixes #######
######################################
random.mx <- read.table(file=random.matrix.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(random.mx) <- chipseq.marks
real.mx <- read.table(file=sample.matrix.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)
real.mx <- as.matrix(real.mx[,c(4:ncol(real.mx))])
colnames(real.mx) <- chipseq.marks






##################################################################
## 3) Data transformation on original matrixes (if required) #####
##################################################################
for ( type in c("real", "random")) {

mx <- get(paste(type, "mx", sep="."))

if (transform.method == "sigmoid" ) {

for (j in c(1:ncol(mx) ) ) {
						cat("column", j, "\n", sep=" ")
						xj <- mx[,j]
						perc95th <- as.numeric(quantile(xj, probs=.95))
					    C = 2/(1+2.72^(-2*xj/perc95th)) - 1 
					    C[is.na(C)] <- 0
						mx[,j] <- C  }
rm(C);rm(xj);gc();					}  else

if (transform.method == "invhyperbolic" ) { mx = log2(mx + sqrt(mx^2+1)) } else

if (transform.method == "raw" ) { mx = mx }

assign(paste(type, "mx",sep="."), mx)  }





#############################################
## 4) Import real and random NMF results ####
#############################################
load( sample.nmfresults.file);
real.estim <- estim; rm(estim);
load( random.nmfresults.file);
random.estim <- estim; rm(estim);
gc();




#####################################################################
## 5) Collect NMF quality results for both random and real data #####
#####################################################################
combodat <- NULL
combosparseness <- NULL

for (type in c("real", "random")) {

estimates <- get(paste(type, "estim", sep="."))
mx <- get(paste(type, "mx", sep="."))

cophcoef.v <- NULL
rss.v <- NULL
basis.sparseness.v <- NULL
coef.sparseness.v <- NULL

# Define vector of Factorization ranks analysed;
ranks <- seq(3, (3+(length(estimates)-1)), 1);

for (k in c(1:length(ranks))) {

rank.factor <- ranks[k];
cat("extract", "nmf results", "for factorization rank", rank.factor, "in", type, "\n", sep=" ")
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



# Build the main dataframe
dat <- data.frame(ranks=ranks,
                 coph=cophcoef.v,
                 rss=rss.v,
                 class=rep(type, length(ranks)),
                 stringsAsFactors=F)
                  
# Build the 'sparseness' dataframe           
sparseness <- data.frame( ranks=rep(ranks,times=2),
						spar=c(coef.sparseness.v, basis.sparseness.v),
						class=paste( type, rep(c('coef', 'basis'), each=length(ranks)), sep=","),
						stringsAsFactors=F)        
						
combodat <- rbind(combodat, dat);
combosparseness <- rbind(combosparseness, sparseness);

						}     # all results from real and random data are now collected in the same data-frames             
					
					
			
			
			
			
					
##############################################			
### 6) Generate Quality Plots PDF  ###########
##############################################
pdf(file=file.path(output.dir, paste(transform.method, statdist, "Real+Random.quality.test.NMF.pdf", sep=".")), width=15, height=4)					     
# Cophenetic correlation coefficient
plot1 <- ggplot(combodat, aes(ranks, coph, group=class, color=factor(class))  ) + geom_line(size=1) + geom_point(size=3) + ggtitle("Cophenetic coeff.");
# Sparseness
plot2 <- ggplot(combosparseness, aes(ranks, spar, group=class, color=factor(class)) ) + geom_line(size=1) + geom_point(size=3) + ggtitle("Sparseness of W and H matrixes");
# RSS
plot3 <- ggplot(combodat, aes(ranks, rss, group=class, color=factor(class)) ) + geom_line(size=1) + geom_point(size=3) +  ggtitle("Residual sum of squares(RSS)");
# Call the 'multiplot' function;
multiplot(plot1,plot2,plot3,cols=3);
dev.off(); 


