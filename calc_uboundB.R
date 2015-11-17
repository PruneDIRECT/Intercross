####################################################################################################
 #------------------------------------------------------------------#
 # Function   :  calc_ubound                                        #
 # Written by :  Dan Finkel                                         #
 # Created on :  10/19/2002                                         #
 # Purpose    :  calculate the ubound used in determing potentially #
 #               optimal hrectangles                                #
 #------------------------------------------------------------------#
.calc_ubound<-function(lengths,fc,hull,szes){

 	ub <- vector()
 	hull_length  <- length(hull)
 	hull_lengths <- lengths[,hull, drop=FALSE]
	
 	for (i in 1:hull_length){
 	    tmp_rects = which(colSums(hull_lengths)< sum(lengths[,hull[i] ] ))
 	    if (length(tmp_rects) > 0){
 	        tmp_f     <- fc[hull[tmp_rects]]
 	        tmp_szes  <- szes[hull[tmp_rects]]
 	        tmp_ubs   <- (tmp_f - fc[hull[i]])/(tmp_szes - szes[hull[i]])
 	        ub[i]     <- max(tmp_ubs);
 	    }else{
 	        ub[i]     <- 1.976e14;
 	    }
 	}
 return(ub)
 }
######################################################################################################