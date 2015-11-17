######################################################################################################
 #------------------------------------------------------------------#
 # Function   :  calc_lbound                                        #
 # Written by :  Dan Finkel                                         #
 # Created on :  10/19/2002                                         #
 # Purpose    :  calculate the lbound used in determing potentially #
 #               optimal hrectangles                                #
 #------------------------------------------------------------------#
 .calc_lbound<-function(lengths,fc,hull,szes){
	
 	lb <- vector()
 	hull_length  <- length(hull)
 	hull_lengths <- lengths[,hull, drop=FALSE]
	
 	for (i in 1:hull_length){
 	    tmp_rects = which(colSums(hull_lengths)> sum(lengths[,hull[i] ] ))
 	    if (length(tmp_rects) > 0){
 	        tmp_f     <- fc[hull[tmp_rects]]
 	        tmp_szes  <- szes[hull[tmp_rects]]
 	        tmp_lbs   <- (fc[hull[i]]-tmp_f)/(szes[hull[i]]-tmp_szes)
 	        lb[i]     <- max(tmp_lbs);
 	    }else{
 	        lb[i]     <- -1.976e14;
 	    }
 	}
 	return(lb)
 }
######################################################################################################