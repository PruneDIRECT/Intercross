#------------------------------------------------------------------#
# Function   :  Main for 3d    - V6                                #
# Written by :  Behrang Mahjani                                    #
# Created on :  11/16/2015                                         #
# Purpose    :  Main										       #
#------------------------------------------------------------------#

#rm(list=ls(all=TRUE))
library(qtl);
library(qtlbook) # sample data
data(ch3b);datafile <- ch3b

hadoop <- 1;
filepath <- "C:/Users/Behrang/Documents/Revolution/Prunedirect-V7/" # (Change here)

source(paste(filepath,"calc_uboundB.R",sep=''))
source(paste(filepath,"calc_lboundB.R",sep=''))
source(paste(filepath,"CallObjFcnB.R",sep=''))
source(paste(filepath,"DIRdivideB.R",sep=''))
source(paste(filepath,"DIRiniB.R",sep=''))
source(paste(filepath,"find_poB.R",sep=''))
source(paste(filepath,"mydirectB.R",sep=''))
source(paste(filepath,"ObjectiveFB.R",sep=''))
source(paste(filepath,"probbcB.R",sep=''))
source(paste(filepath,"rssB.R",sep=''))
source(paste(filepath,"ProblemInt.R",sep=''))
source(paste(filepath,"FindSym.R",sep=''))
#datafile <- read.cross("csv",filepath, "hyper2.csv",genotypes=c("BB","BA"));  # your data in format of RQTL (Change here) if you want to change it to your own dataset
	
	
#change the objective function in the file objectiveFB.R, the name of the function is ObjectiveF3dInt 


Problem		 <- list("f"="ObjectiveF3dInt")
ideal_step   <- 0.1# The ideal discretization step used for the mesh 
chrom	     <- c(4,5,6) # Which chromosomes to use? (Change here)
Problem.data <- ProblemInt(datafile   = datafile,
						   chrom      = chrom,
						   ideal_step = ideal_step,
						   delta	  = 0.5)
results <- mydirect(Problem				= Problem,
					Problem.data		= Problem.data,
 					showits				= "none",
 					verbose				= T,
					permutation			= 0,
					calculated_rss_saved= matrix(0,1,1))
#output: X1:..., X2..., X3,... , Number of function evaluations:...