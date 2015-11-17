###########################################################################################################      
 #------------------------------------------------------------------#
 # Function:   DIRini                                               #
 # Written by: Dan Finkel                                           #
 # Created on: 10/19/2002                                           #
 # Purpose   : Initialization of Direct                             #
 #             to eliminate storing floating points                 #
 # This function is modified by B. Mahjani on 3/26/2015             # 
 #------------------------------------------------------------------#
 .DIRini<-function(Problem,
					Problem.data,
					n,a,b,
 					param.names=c(1:length(a)),
					c, fc, szes,
					calculated_rss_saved,
					permutation,
					minval_prune,
					...){
 
 # DIRECT begins the optimization by transforming the domain of the problem into the unit hyper-cube.
 # The algorithm works in this normalized space, referring to the original space only when
 # making function calls. The center of this space is c1, and we begin by fnding f(c1).
 
 	# Start by calculating the thirds array
 	# Here we precalculate (1/3)^i which we will use frequently
 	l_thirds<-rep(NA, (Problem.data$maxdeep))
 	l_thirds[1] <- 1/3
 	for (i in 2:Problem.data$maxdeep){
 	   l_thirds[i] <- (1/3)*l_thirds[i-1];
 	}
 	
	# lengths=  length array stores number of slices in each dimension for each rectangle. 
	# Dimension will be rows; Each rectangle will be a column
 	
 	# First rectangle is the whole unit hyperrectangle
 	l_lengths			<- matrix(0,n,1);
 	rownames(l_lengths) <- param.names
 	
 	
 	# Store size of hyperrectangle in vector szes
 	szes		<- 1
 	names(szes) <- "start"
 	
 	# First element of c is the center of the unit hyperrectangle
 	l_c			  <- matrix(1/2,n,1)
 	colnames(l_c) <- "start"
 	rownames(l_c) <- param.names
 	 			
 	# First element of f is going to be the function evaluated at the center of the unit hyper-rectangle.
 	func.List<- .CallObjFcn(Problem,
							Problem.data,
 							point.x=l_c[,1, drop=FALSE],
							a, b,...)

 	l_fc			 <- func.List$fcn_value
 	fcncounter		 <- 1
	true_fcncounter  <- 1 # fcncounter changes when wer are pruning. So we need a new variable to save the number of function evaulations
 		
 	# Initialize minval and point.xatmin to be center of hyper-rectangle !  (NOT in the original intervals!!!!)
 	point.xatmin <- l_c[,1, drop=FALSE]      
 	minval		 <- l_fc[1]
	
    perror = 2

 	# Initialize History
 	History			 <- t(matrix(c(0,0,0,0)))  
    colnames(History)<- c("Iteration Nr", "Function Count - Pruned","Exact Funcion Count", "f_min")
	
	#if ((fc < calculated_rss_saved[1,1]) && (permutation))	prune <- 1 #!!!!! check this
	#else prune <- 0
    if ((l_fc[1] > calculated_rss_saved[1,1]) && (permutation))	prune <- 1 #!!!!! check this - The idea is to check if the first box should be pruned.
	else prune <- 0
	
	return(list(thirds		= l_thirds,
 				lengths		= l_lengths,
				c			= l_c,
 				fc			= l_fc,
				minval		= minval,
				point.xatmin= point.xatmin,
				perror		= perror,
 		        History		= History,
 				szes		= szes,
				fcncounter	= fcncounter,
				true_fcncounter = true_fcncounter,
				prune		= prune ))
 }
 
 ####################################################################################################