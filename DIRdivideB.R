 ######################################################################################################
 #------------------------------------------------------------------#
 # Function   :  DIRdivide                                          #
 # Written by :  Dan Finkel                                         #
 # Created on :  10/19/2002                                         #
 # Purpose    :  Divides rectangle i that is passed in              #
 # This function is modified by B. Mahjani on 3/26/2015             #
 #------------------------------------------------------------------#
 
.DIRdivide<-function(a,
					 b,
					 Problem,
					 Problem.data,
					 index,
					 thirds,
					 lengths, fc, c,
 					 p_fcncounter,
					 true_fcncounter,
 					 szes,
					 showits.flag=TRUE,
 					 itctr="",
 					 i="",
 					 prune,
					 minval_prune,
					 permutation,
					 calculated_rss_saved,...){ 
								    						
# prune cases
# 0: no prune
# 1: prune based on rss
# 2: prune based on symmetry
# 3: prune based on mesh size
		
 fcncounter				<- p_fcncounter
 new_minval_prune		<- minval_prune
 pruned_c				<- 0
 flag_new_minval_prune	<- 0
 flag_center_pruned		<- 0
 
 starting_size <- length(fc)
 #---1. Determine which sides are the largest 
 li     <- lengths[,index]  
 biggy  <- min(li)
 ls     <- which(li==biggy)
 lssize <- length(ls)
 
 #---2. Evaluate function in directions of biggest size to determine which direction to make divisions
 oldc       <- c[,index,drop=FALSE]
 delta      <- thirds[biggy+1]   

 n <- length(Problem.data$pheno)
 

 # Add or create new centers  
 newc_left  <- newc_right <- matrix(rep(oldc, lssize), ncol=lssize) 
 
 # Initialize
 f_left  <- rep(NA, lssize)
 f_right <- rep(NA, lssize) 
 
 # For each dimention (parameter) in ls create new centers : left and right from the old center
 for (i in 1:lssize){
     lsi  <- ls[i]
     
     # c_i +/- delta*e_i, e_i has at the ith position 1, rest=0
     newc_left[lsi,i]  = newc_left[lsi,i] - delta;
     newc_right[lsi,i] = newc_right[lsi,i] + delta;

     # f(new_centers_left)
     func.left.list<- .CallObjFcn(Problem,
								  Problem.data,
  								  point.x=newc_left[,i,drop=FALSE],
								  a, b,...)
 	 f_left[i]<-  func.left.list$fcn_value
      
     # f(new_centers_right)
     func.right.list<- .CallObjFcn(Problem,
								   Problem.data,
 								   point.x=newc_right[,i,drop=FALSE],
								   a,b,...)
 	 f_right[i]<-  func.right.list$fcn_value
 
     fcncounter <- fcncounter + 2
	 true_fcncounter <- true_fcncounter + 2	
 }
   
 #--- 2.1 Calculate w - min of function in new centers 
 # like in DIRECTUserGuide: w_i := min( f_right, f_left ) for i in 1:N
 # Best function valueS ! in the largest space
 # This means for each dimention find separat min ! 
 w = apply(cbind(f_left, f_right), 1, min)
 
 
 #--- 3. Sort w for division order 
 tmp.sort <- sort(w,index.return=TRUE)
 V <- tmp.sort$x; 
 order <- tmp.sort$ix
 
 #--- 4. Make divisions in order specified by order 
 for (i in 1:length(order) ){

    newleftindex  <- p_fcncounter + 2*(i-1)+1
    newrightindex <- p_fcncounter + 2*(i-1)+2
	
    #--- 4.1 create new rectangles identical to the old one 
    oldrect <- lengths[,index, drop= FALSE]
    lengths <- cbind(lengths, oldrect, oldrect)
 
    # old, and new rectangles have been sliced in order(i) direction
    lengths[ls[order[i]],newleftindex ] <-  lengths[ls[order[i]],index] + 1
    lengths[ls[order[i]],newrightindex] <-  lengths[ls[order[i]],index]  + 1;
    lengths[ls[order[i]],index]         <-  lengths[ls[order[i]],index]  + 1;
 
    # add new columns to c
    c <- cbind(c, newc_left[,order[i]], newc_right[,order[i]] )
    colnames(c)[(ncol(c)-1):ncol(c)] <- c("left", "right")
		 
    # Add new values to f
    fc <- c(fc, f_left[order[i]], f_right[order[i]] )
 
	# maxdeep is the number of different sizes of lengths we have. 
	# how large is the center box? left and right boxes are the same size
	
	L		 <- rep(0,Problem.data$d) # this is the length of a box in each direction
	L_scaled <- rep(0,Problem.data$d) # ... scaled back to real values
	for (k in 1:Problem.data$d){
	
 		L[k]	   <- lengths[k,index]  
		
    	if(L[k]==0) delta <- 1
		else delta <- thirds[L[k]]	
		
		L[k] <- delta
		
		L_scaled[k]<-(b[k]-a[k])*delta #scale it back
	}
	
   	x <- sum(L_scaled)/(2*Problem.data$d) # explain why /4 - I think it is /2*dim
	
	# Prune for fmin
	leftprune <- rightprune <- 0
	
    if(!permutation){calculated_rss <- rss(n,x,minval_prune);} 
	else calculated_rss <- calculated_rss_saved[l1+1,l2+1] # this is wrong!!!!
	
	if (f_left[order[i]] >  calculated_rss)
		{leftprune <- 1;}
	
	if (f_right[order[i]] >  calculated_rss)
		{rightprune <- 1;}
		
	if (fc[index]  > calculated_rss) # center - careful with pruning the center, the index should be saved
		{prune[index] <- 1;flag_center_pruned <- 1;}
	#-- Prune for symmetry
	# value : 2
	if (Problem.data$is.symmetric)
	{
		if (FindSym(center_box=newc_left[,order[i]],size_box=L)) # left
			{leftprune <- 2;}
		if (FindSym(center_box=newc_right[,order[i]],size_box=L)) # right
			{rightprune <- 2;}
		if (FindSym(center_box=c[,index],size_box=L)) # center, careful with pruning the center, the index should be saved
			{prune[index] <- 2;flag_center_pruned <- 1;}
	}
    # Prune for maxdeep/mesh
	if ((lengths[ls[order[i]],newleftindex ]+1) >= Problem.data$maxdeep)
		{leftprune <- 3;}
	
	if ((lengths[ls[order[i]],newrightindex ]+1) >= Problem.data$maxdeep)
		{rightprune <- 3; flag_center_pruned <- 1;}
		
	if ((lengths[ls[order[i]],index ]+1) >= Problem.data$maxdeep) # center
		{prune[index] <- 3;flag_center_pruned <- 1;} # careful with pruning the center, the index should be saved		
	
	prune <- c(prune,leftprune,rightprune)	
	#------   		
	tmp.szes.l<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newleftindex, drop=FALSE])
 	tmp.szes.l<- 1/2*max(svd(tmp.szes.l)$d)	
 		
 	tmp.szes.r<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newrightindex, drop=FALSE])
 	tmp.szes.r<- 1/2*max(svd(tmp.szes.r)$d)	
 		   
    szes<-c(szes, tmp.szes.l, tmp.szes.r )   
    names(szes)[(length(szes)-1):length(szes)]<- c("left", "right")
 } #end of for
   
 
 ## plot old and new centers ####
 if (showits.flag ){
 	my.col <- rep("black", ncol(c))
 	my.col[grep("left",colnames(c))] <- "blue"
 	my.col[grep("right",colnames(c))] <- "red"
 	
 	if (nrow(c)==1){
 		toPlot<- rbind(c, fc) 
 		#plot(as.data.frame(t(toPlot)), pch=20, type="p", xlim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  ) )
		plot(as.data.frame(t(toPlot)), pch=20, type="p",xlim=c(0,1),ylim=c(0,1))# xlim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  ) )
 		
 		#text(x=toPlot[1,],y=toPlot[2,], round(fc,1), pos=1, cex=0.3)		
 		#text(x=toPlot[1,],y=toPlot[2,], c(1:ncol(c)), pos=3, cex=0.3, col="#008080")	
 		#legend("topright", c("old center", "new left", "new right", "center's number"), fill=c("black", "blue", "red", "#008080"), cex=0.5)
 	} # end for 1 tuning parameter
 	
 	# for 2 tuning parameters
 	if (nrow(c)==2){
 	toPlot <- c
	plot(as.data.frame(t(toPlot)), cex=0.6,pch=20, type="p",xlim=c(0,1),ylim=c(0,1))#, xlim=c(0,1), ylim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  )  )
 	#plot(as.data.frame(t(toPlot)), pch=20, type="p", xlim=c(0,1), ylim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  )  )
  	#text(x=toPlot[1,],y=toPlot[2,], round(fc,1), pos=1, cex=0.3)		
 	#text(x=toPlot[1,],y=toPlot[2,], c(1:ncol(c)), pos=3, cex=0.3, col="#008080")	
 	#legend("topright", c("old center", "new left", "new right", "center's number"), fill=c("black", "blue", "red", "#008080"), cex=0.5)
 	} # end for 2 tuning parameters
 }
 ## end of plot old and new centers #### 		

 # How many boxed are created
 tmp.szes.ind <-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,index, drop=FALSE])
 tmp.szes.ind <- 1/2*max(svd(tmp.szes.ind)$d)	
 szes[index]  <- tmp.szes.ind

 #!!!!!! this part is not efficent. Choose the indexes more wisely! 
 prune_index <- c()
 updated_size <- length(fc)
 for (i in (1):(updated_size))
 {	
	if(prune[i]>0)
	{	
		if(prune[i] == 3)# prune for mesh. If we pruned a fmin, then we should save its value.
		{
			if (new_minval_prune >= fc[i])
			{
				new_minval_prune <- fc[i]
				pruned_c		 <-	c[,i]
				flag_new_minval_prune <- 1 				
			}
		}
		fcncounter  <- fcncounter - 1
		prune_index <- c(prune_index,i)
	}
 }

if(length(prune_index))
{
	lengths		<- lengths[,-prune_index,drop=F] # important to have drop=F, otherwise the dimention drops to a 1 dmin matrix
	fc			<- fc[-prune_index]
	c			<- c[,-prune_index,drop=F]
	szes		<- szes[-prune_index]
	prune		<- prune[-prune_index]		
}

 return(list(lengths=lengths,
			 fc=fc,c=c,
 			 szes=szes,
			 fcncounter=fcncounter,
			 true_fcncounter=true_fcncounter,
			 prune=prune,
			 minval_prune=new_minval_prune,
			 pruned_c=pruned_c,
			 flag_new_minval_prune=flag_new_minval_prune,
			 flag_center_pruned=flag_center_pruned))
}
 
 ####################################################################################################
