 ######################################################################################################
 #--------------------------------------------------------------------#
 # Function   :  FindSym                                              #
 # Written by :  Behrang Mahjani                                      #
 # Created on :  4/08/2015	                                          #
 # Purpose    :  Check if a box is out of the region                  #
 # Parts of this function is from symDIRECT, a Matlab code  writen    #
 # by Remigijus. Paulavi?ius, presented in: 					      #
 # https://github.com/alimaitevelina/matlabas/blob/master/symDIRECT.m #										
 #--------------------------------------------------------------------#

FindSym <- function(center_box,size_box){
	symPrune <- 0 # is the box out of the region 
	size_box <- size_box/2
	n		 <- length(size_box)
	V		 <- Problem.data$V # each column is a vertix
	real_V	 <- matrix(0,n,2^n)
	
	symPruneFlag <- 0
 	
	for (i in 1:2^n){		
		for(j in 1:n){
			if (V[j,i] == 0){
				real_V[j,i] <- center_box[j] - size_box[j];
            }else{
				real_V[j,i] <- center_box[j] + size_box[j];
            }
         }
		FlagVertex <- 0
        for(j in 1:(n-1)){               
			if (real_V[j,i] < real_V[j+1,i]){ # is the vertex above the line?
				FlagVertex <- FlagVertex + 1
			}
			else break			
		}
		if (FlagVertex == (n-1)) # vertex i above the line 
			{symPruneFlag <- symPruneFlag + 1}
	}
	if (symPruneFlag == (2^n)) symPrune <- 1 # all of the vertices above the line, then prune it

return(symPrune)
}
######################################################################################################
