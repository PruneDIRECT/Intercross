##################################################################################################################
#------------------------------------------------------------------#
# Function   :  ProblemInt                                         #
# Written by :  Behrang Mahjani                                    #
# Created on :  3/26/2015                                          #
# Purpose    :  Extract the information from the chromosomes       #
#------------------------------------------------------------------#
#???? What is needed: 1. add cumulative effect of delta for symetric cases
ProblemInt<-function(datafile,chrom,ideal_step,delta = 0)
{
# Extract genotypes + phenotypes
# Save the probabilites in a list
# Creat a mesh and save it in a list
# pheo:datafile$pheno contains both genotypes and locations,
# Choose only genotypes for chrom i: datafile$geno[[i]][[1]]
# Choose locations: datafile$geno[[i]][[2]]
# bounds are for the mesh, ends are for the chrom. We don't use bounds because our mesh never starts from zero. because we work with the midel points

# delta is the bound for the symmetric case
chrom_length <- length(chrom)

#Check if the problem is symetric. In case of multidimentional scan on one chromosome
is.symmetric <- 0  

if (diff(chrom)[chrom_length-1] ==0){
	is.symmetric <- 1	
}

step3p <- rep(0,chrom_length)
step   <- rep(0,chrom_length)
bounds <- matrix(0,chrom_length,2)
#bounds variable is based on two ends of mesh. ends is the varibale for the first and the last marker for scaling
ends    <- matrix(0,chrom_length,2)

if(is.symmetric){# if symmetry
	#cut the search space to avoid between markers QTL - 
	locations		<- datafile$geno[[chrom[1]]][[2]]
	locations <- locations - locations[1] # start the locations from 0 because of prob function from RQTL
	numberOfMarkers <- length(locations)
		
	## 1- Calculate the ideal step size :
	step3p[1:chrom_length] <- round(log((locations[numberOfMarkers] - locations[1])/ideal_step,base=3))
		
	## 2- calculate n, where (n*deltax) /(n*3^i deltax) = delta/L
	n <- floor(3^step3p[1]*delta/(locations[numberOfMarkers] - locations[1]))
	if(n==0) { n <- 1;}
		
	## 3- calculate deltax
	deltax	<- (locations[numberOfMarkers] - locations[1])/(3^step3p[1]+n)
	step[1:chrom_length]	<- deltax
	
	## 5- pre- Calculate the probabilites
	prob     <- probbc(datafile, chrom[1], step = 0.5*step[1], err=0.001)
	prob_removed_markers <- prob[,!(names(prob[1,,1]) %in% names(locations[-1])),1 ]# remove the markers - keep the first one 			
	
	# 6.1:(chrom_length-1)
	startCh	<- n*deltax + locations[1]
	endCh	<- locations[numberOfMarkers]
	meshSize <- (endCh - startCh)/step[1] 
	mesh	 <- c(startCh + step[1]*0.5, startCh + step[1]*0.5 + step[1]*(1:(meshSize-1)))
	

	prob_removed_markers_2  <- prob_removed_markers[,(2*n+1+locations[1]):dim(prob_removed_markers)[2]] # ???add the dimention later
	prob_removed_extras <- prob_removed_markers_2[, seq(2,2*length(mesh),by=2)]# remove the extra points
	for(i in 1:(chrom_length-1)){
		if (i == 1){
			P <- list(prob_removed_extras)  
			meshSave <-list(mesh)
		}else{ 
			P[[i]] <- prob_removed_extras   
			meshSave[[i]] <- mesh
		}
		
		bounds[i,1] <- mesh[1]
		bounds[i,2] <- mesh[length(mesh)]
		
		ends[i,1]   <- startCh
		ends[i,2]   <- endCh
	}
	
	# 6.2:(chrom_length-1)
	startCh		<- locations[1]
	endCh		<- locations[numberOfMarkers]-n*deltax
	meshSize <- (endCh - startCh)/step[1] 	
	mesh	 <- c(startCh + step[1]*0.5, startCh + step[1]*0.5 + step[1]*(1:(meshSize-1)))
	
	prob_removed_markers_2 <- prob_removed_markers[,(1):(dim(prob_removed_markers)[2]-2*n)] 
	prob_removed_extras <- prob_removed_markers_2[, seq(2,2*length(mesh),by=2)]# remove the extra point
	P[[chrom_length]] <- prob_removed_extras    
	meshSave[[chrom_length]] <- mesh
	
	bounds[chrom_length,1] <- mesh[1]
	bounds[chrom_length,2] <- mesh[length(mesh)]
	
	ends[chrom_length,1]   <- startCh
	ends[chrom_length,2]   <- endCh
	
}else{ # non-symemteric 
	for(i in 1:chrom_length){
		locations		<- datafile$geno[[chrom[i]]][[2]]
		locations		<- locations - locations[1] # start the locations from 0 because of prob function from RQTL
		numberOfMarkers <- length(locations)
		step3p[i]		<- round(log((locations[numberOfMarkers] - locations[1])/ideal_step,base=3))# Transform the ideal_Step size in form of 3^-i
		
		# The mesh points should be in form of 0.5*3^-i*(2k-1). 
		#	With how many 3^-i we can cover the chromosome?
		#   So we are finding i where 3^-i is almost equal to ideal_step
		step[i]			<- (locations[numberOfMarkers] - locations[1])*3^(-step3p[i])
		prob		    <- probbc(datafile, chrom[i], step = step[i]*0.5 , err=0.001)# there are some extra point because the final messh is not uniform

		meshSize		<- (locations[numberOfMarkers] - locations[1])/step[i] 	
		mesh			<- c(locations[1] + step[i]*0.5, locations[1] + step[i]*0.5 + step[i]*(1:(meshSize-1)))

		prob_removed_markers <- prob[,!(names(prob[1,,1]) %in% names(locations[-1])),1 ]# Remove the markers
		prob_removed_extras <- prob_removed_markers[, seq(2,2*length(mesh),by=2)]# remove the extra points
	    if (dim(prob)[3] == 2) # backcross - Enough with P_AA
		{
			if (i == 1){
				P <- list(prob_removed_extras)
				meshSave <-list(mesh)
			}else{ 
				P[[i]] <- prob_removed_extras 
				meshSave[[i]] <-(mesh)
			}
		}else{# inter cross
			P_AA <- prob_removed_extras
			prob_removed_markers <- prob[,!(names(prob[1,,1]) %in% names(locations[-1])),2 ]# Remove the markers
			P_AB <- prob_removed_markers[, seq(2,2*length(mesh),by=2)]# remove the extra points
			prob_removed_markers <- prob[,!(names(prob[1,,1]) %in% names(locations[-1])),3 ]# Remove the markers
			P_BB <- prob_removed_markers[, seq(2,2*length(mesh),by=2)]# remove the extra points
		}
			if (i == 1){
				PAA <- list(P_AA)
				PAB <- list(P_AB)
				PBB <- list(P_BB)
				meshSave <-list(mesh)
			}else{ 
				PAA[[i]] <- P_AA
				PAB[[i]] <- P_AB
				PBB[[i]] <- P_BB
				meshSave[[i]] <-(mesh)
			}
		bounds[i,1] <- mesh[1]
		bounds[i,2] <- mesh[length(mesh)]
		ends[i,1]   <- locations[1]
		ends[i,2]   <- locations[numberOfMarkers]
	}
}#end of else
if (dim(prob)[3] == 3) P <- list(PAA=PAA,PAB=PAB,PBB=PBB) 
#- binary representation of vertices for dim n, used in case of symmetry
n <- chrom_length 

if(is.symmetric){
	b	<- matrix(0,n,1)       
    k	<- 1
	V	<- matrix(0,n,2^n)
    V[,k] <- b
    k	  <- k + 1
	flag  <- 1
    while(flag){
		j <- n
        while ((j>0) && (b[j] == 1)){
			b[j] <- 0
            j	 <- j - 1
        }
        if(j > 0){
			b[j] <- 1
            V[,k]<- t(b)
            k	 <- k+1
        }else{
            flag <- 0 
        }
    }

}else V <- 0
#---
pheno	 <- datafile$pheno[[1]]
pheno <- pheno + rnorm(length(pheno),mean=mean(pheno)+0.1*mean(pheno),sd=sd(pheno)+0.5*sd(pheno));

mu		 <- mean(pheno)
variance <- var(pheno)
maxdeep  <- min(step3p) + 1
colnames(bounds)<-c("lower", "upper")

Problem.data = list(pheno		 = pheno,
					mu			 = mu,
					variance	 = variance,
					prob		 = P,
					mesh		 = meshSave,
					bounds		 = bounds,
					is.symmetric = is.symmetric,
					maxdeep		 = maxdeep,
					step		 = step,
					d			 = chrom_length,
					V			 = V,
					ends		 = ends)
				
return(Problem.data)
}				 
##################################################################################################################
