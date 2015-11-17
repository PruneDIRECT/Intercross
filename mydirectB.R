################################################################################################################## 
#------------------------------------------------------------------#
# Function   :  mydirect                                           #
# Written by :  Dan Finkel                                         #
# Created on :  10/19/2002                                         #
# Purpose    :  The main DIRECT algorithm		         		   #
# This function is modified by B. Mahjani on 11/27/2015             #
#------------------------------------------------------------------#

mydirect<-function(	Problem,
					Problem.data,	
	                # options
					maxits   =   10000,      #% maximum of iterations
					maxevals =   500000,     #% maximum # of function evaluations
					#maxdeep =    7,#1000,   #f% maximum number of side divisions
					ep		 =   1e-4,       #% global/local weight parameter.
					tol		 =   0.01,       #% allowable relative error if f_reach is set
					showits  =   c("none", "final", "all"), #%  plot iteration stats: none, final iteration, all iterations 
					verbose  =   TRUE, 	 #% print  iteration stats: none, final itertation, all iterations              
					# plot parameter
					pdf.name =   NULL, 
					pdf.width=   12, pdf.height=12,
					my.mfrow =   c(1,1), 
					permutation,
					per_minval= Inf,
					calculated_rss_saved,
					 ...){
#-- Initialize the variables --------------------------------------

lengths <- c_value <- fc <- prune <- vector()

szes		 <- vector()
om_lower     <- Problem.data$ends[,1]#Problem.data$bounds[,"lower", drop=FALSE]
om_upper     <- Problem.data$ends[,2]#Problem.data$bounds[,"upper", drop=FALSE]
fcncounter   <- 0
perror       <- 0
itctr        <- 1
done         <- 0
n            <-nrow(Problem.data$bounds)
true_fcncounter <-0

minval_prune <- per_minval;#input from main. if no premutation, then it is equal to Inf

#-- Pre-allocate memory for storage vectors
lengths    = matrix(0,n,c(maxevals + floor(.10*maxevals)))
c_value    = lengths;
fc         = matrix(0,1,c(maxevals + floor(.10*maxevals))) 
prune      = matrix(0,1,c(maxevals + floor(.10*maxevals))) 
szes       = fc


#-- Call DIRini ---------------------------------------------------
# define the multi-dim hypercube with center point c, fc=f(c) 

DIRini.list<-.DIRini(Problem,
					 Problem.data,
					 n,
 					 a=Problem.data$bounds[,"lower"],
 					 b=Problem.data$bounds[,"upper"],
					 param.names= rownames(Problem.data$bounds),
					 c=c_value,
					 fc=fc,
  					 szes=szes,
					 calculated_rss_saved,
					 permutation,
					 minval_prune,
					 ...)
          
thirds <- DIRini.list$thirds

#lengths is the  length array. It stores the number of slices in each dimension for each rectangle. 
#Dimensions are rows; Each rectangle is a column
lengths		 <- DIRini.list$lengths
minc_pruned <- c_value <-DIRini.list$c 
minfc_pruned  <- fc <-DIRini.list$fc

minval		 <-DIRini.list$minval
point.xatmin <-DIRini.list$point.xatmin
perror		 <-DIRini.list$perror
History		 <-DIRini.list$History
szes		 <-DIRini.list$szes                 # number of regions
fcncounter	 <-DIRini.list$fcncounter
prune        <- DIRini.list$prune;
true_fcncounter  <- DIRini.list$true_fcncounter
 
ret_point.xatmin <- point.xatmin

Hfcncounter <- 0

if(showits !="none" & !is.null(pdf.name)) { 
	pdf(pdf.name, pdf.width, pdf.height)
}
if(showits != "none")
	par(mfrow=my.mfrow)

#-- MAIN LOOP -----------------------------------------------------
minval		 <- fc[1] 
minfc_pruned <- fc[1] 
mintemp		 <- minval
 
while (perror > tol){ #check perror later
	#if(permutation && prune[1]) {if (verbose) print("First box Pruned");break;}# Check if the first boxed should to be pruned
	mintemp <- minval 
	
   	#-- Create list S of potentially optimal hyper-rectangles
    S <- .find_po(fc	  = fc,
				  lengths = lengths,
				  minval  = minval,
 				  ep	  = ep,
 				  szes	  = szes,
 				  prune   = prune)

	# if we don't find potentially hyper-rectangles --> break!
	if (ncol(S)==0) { if (verbose) {print("No potentially hyper-rectangles left")};break;}
	o <- order(S[1,])
	S <- rbind(S[1,o],S[2,o])
	
   #-- Loop through the potentially optimal hyper-rectangles ------
   for (i in 1:ncol(S)){
	
     # plot options: don't plot if not requested
     if ( (showits =="none") ) { 
     	showits.flag <-  FALSE
     } else {
	     # plot  last iteration's step  if requested
	     if ((showits =="final")) showits.flag <- ifelse ( (i == ncol(S)), TRUE, FALSE )
	     # plot all iterations' steps  if requested
	     if ( (showits =="all") ) showits.flag <-  TRUE
     }
	
	if(!permutation) minval_prune <- mintemp #if no permutation, use the last fmin for pruning rss
	
	# There are 3 different Prune
	# 1- Mesh - min should be saved if min < minval
	# 2- Symetirc - min should be saved 
	# 3- rss - min should be saved 
	
    tmp.list.divide<- .DIRdivide(a			 = Problem.data$bounds[,1],
 								 b	         = Problem.data$bounds[,2],
								 Problem	 = Problem,
								 Problem.data= Problem.data,
								 index		 = S[1,i],
 								 thirds		 = thirds,
 								 lengths	 = lengths,
								 fc			 = fc,
								 c			 = c_value,
								 p_fcncounter= fcncounter,
								 true_fcncounter=true_fcncounter,
 								 szes		 = szes,
 								 showits.flag= showits.flag,
 								 itctr		 = itctr,
 								 i			 = i,
								 prune		 = prune,
 								 minval_prune= minval_prune,
								 permutation =permutation,
								 calculated_rss_saved=calculated_rss_saved,
								 ...)

	
	  lengths	<- tmp.list.divide$lengths
	  fc		<- tmp.list.divide$fc
	  c_value   <- tmp.list.divide$c
	  szes		<-  tmp.list.divide$szes
	  fcncounter<-  tmp.list.divide$fcncounter
   	  prune		<- tmp.list.divide$prune
	  true_fcncounter <- tmp.list.divide$true_fcncounter


	  if (!ncol(lengths)){ 
         if (verbose) print("The whole space is pruned.")
         break
      }	
	  # if the min was pruned, we must save it - This can happned when we prune for mesh
		if(tmp.list.divide$flag_center_pruned){n = ncol(S);S[1,i:n] <- S[1,i:n] - 1 ;}
	    if(tmp.list.divide$flag_new_minval_prune)# find the min of the pruned boxes - if a min of a pruned box is smaller than fmin. save the min
		{
			if(mintemp >= tmp.list.divide$minval_prune){
				minfc_pruned  <- tmp.list.divide$minval_prune
				minc_pruned   <- (om_upper - om_lower)*tmp.list.divide$pruned_c + om_lower
				mintemp       <- minfc_pruned
				permutation   <- 0
			}
		}
  }
    
   #-- update minval, point.xatmin --------------------------------------

   if(length(fc)){
		minval	 <- min(fc)
		fminindex<- which.min(fc)
	}
	else{
		minval	<- Inf
		fminindex<- Inf
	}


if (minval < minfc_pruned)
 {
   fminindex	<- which.min(fc)
   point.xatmin <- (om_upper - om_lower)*c_value[,fminindex] + om_lower
 }
else
{
	minval		 <- minfc_pruned
	point.xatmin <- minc_pruned
}

   if((minval < per_minval)&&(permutation)) # the fmin from permutation is not useful anymore, because we found one smaller than that
	{
		per_minval  <- minval;
		permutation <- 0; # permutation = 1, as far as we want to use the old rss values based on the old fmin. New fmin, means new rss are required 
	}
   ret_minval		<- minval;
   ret_point.xatmin <- point.xatmin;

 
   #-- See if we are done ------------------------------------------

      #-- Have we exceeded the maxits?
      if (itctr >= maxits){
         if (verbose) print("Exceeded max iterations. Increase maxits")
         done <- 1
      }
	  # for permutation, when we are pruning, 
	  else if (!ncol(lengths)){ 
         if (verbose) print("The whole space is pruned.-")
         done <- 1
      }
      #-- Have we exceeded the maxevals?
      else if (fcncounter > maxevals){
         if (verbose) print("Exceeded max fcn evals. Increase maxevals")
         done <- 1
      }
	  else if (min(as.vector(lengths)) >= (0.5*Problem.data$maxdeep) ){ ###############Stoping criteria for PruneDIRECT
      #-- We've exceeded the max depth
		L_History <- length(History);
		if (abs(History[L_History]-History[L_History-1])<0.000001)
		{
			if (verbose) print("cover most of the space")
			perror = -1
		}
      }
	  if (done == 1) perror = -1
   
   #-- Store History
   History<-rbind(History,
                 c(itctr, fcncounter,true_fcncounter, minval))
  
  #-- show iteration stats
  if (verbose & ret_point.xatmin==2)  print(paste("Iter:", itctr, "f_min:", minval, "fn evals:","X1",ret_point.xatmin[1],"X2",ret_point.xatmin[2], fcncounter,true_fcncounter, sep="   "))
	else  print(paste("Iter:", itctr, "f_min:", minval, "fn evals:","X1",ret_point.xatmin[1],"X2",ret_point.xatmin[2],"X3",ret_point.xatmin[3], fcncounter,true_fcncounter, sep="   "))
  Hfcncounter <- c(Hfcncounter, fcncounter)
  itctr  = itctr + 1
} # end  of while (perror > tol)

if(showits !="none" & !is.null(pdf.name)) dev.off()

#-- Return values   
#-- return x
final_point.xatmin <- ret_point.xatmin;

#-- chop off (abschneiden) 1st row of History
History<-History[-1,]


return (list(final_point.xatmin=final_point.xatmin,
			 minval =minval, 
			 c=c_value,
 			 fc=fc, 
			 History=History,
 			 Hfcncounter= Hfcncounter))
}

###########################################################################################################
