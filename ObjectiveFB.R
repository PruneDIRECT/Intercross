##################################################################################################################
#------------------------------------------------------------------#
# Function   :  Objective function                                 #
# Written by :  Behrang Mahjani                                    #
# Created on :  3/26/2015                                          #
# Purpose    :												       #
#------------------------------------------------------------------#

# Change this to multidimentional 

ObjectiveF <- function(x,Problem.data){ # 2 dimentional
	
	probAA1_local <- Problem.data$prob[[1]]
	probAA2_local <- Problem.data$prob[[2]]
	pheno_local   <- Problem.data$pheno
	mesh1_local   <- Problem.data$mesh[[1]]
	mesh2_local   <- Problem.data$mesh[[2]]
	variance_local<- Problem.data$variance 
	step1		  <- Problem.data$step[1]
    step2		  <- Problem.data$step[2]		
	
	x1 <- x[1];
	x2 <- x[2];

	if (x[1] == mesh1_local[1]){
		location1 <- 1
	} else
		location1 <- (x1-mesh1_local[1])/step1 + 1
	
	if (x[2] == mesh2_local[1]){
		location2 <- 1
	} else
		location2 <- (x2-mesh2_local[1])/step2 + 1

	Px1 <- probAA1_local[, as.integer(as.character(location1))]-0.5;
	Px2 <- probAA2_local[, as.integer(as.character(location2))]-0.5;
	
	Px1x2 <- Px1*Px2; # If they are independet. 
	
	lm.x <- lm(pheno_local~Px1+Px2+Px1x2); 
 
	rss <- summary(lm.x)$sigma;
	sigma2 <- rss^2;
	f <-  -log(abs(-sigma2+variance_local)/variance_local)# make it positive

	return(f);
}
##################################################################################################################
#3 - Dim, backcross

ObjectiveF3d <- function(x,Problem.data){ 
	
	probAA1_local <- Problem.data$prob[[1]]
	probAA2_local <- Problem.data$prob[[2]]
	probAA3_local <- Problem.data$prob[[3]]
	
	pheno_local   <- Problem.data$pheno
	
	mesh1_local   <- Problem.data$mesh[[1]]
	mesh2_local   <- Problem.data$mesh[[2]]
	mesh3_local   <- Problem.data$mesh[[3]]
	
	variance_local<- Problem.data$variance 
	
	step1		  <- Problem.data$step[1]
    step2		  <- Problem.data$step[2]		
	step3		  <- Problem.data$step[3]		
	
	x1 <- x[1];
	x2 <- x[2];
	x3 <- x[3];

	if (x[1] == mesh1_local[1]){
		location1 <- 1
	} else
		location1 <- (x1-mesh1_local[1])/step1 + 1
	
	if (x[2] == mesh2_local[1]){
		location2 <- 1
	} else
		location2 <- (x2-mesh2_local[1])/step2 + 1
		
	if (x[3] == mesh3_local[1]){
		location3 <- 1
	} else
		location3 <- (x3-mesh3_local[1])/step3 + 1

	Px1 <- probAA1_local[, as.integer(as.character(location1))]-0.5;
	Px2 <- probAA2_local[, as.integer(as.character(location2))]-0.5;
	Px3 <- probAA3_local[, as.integer(as.character(location3))]-0.5;
	
	Px1x2 <- Px1*Px2; # If they are independet. 
	Px1x3 <- Px1*Px3; # If they are independet. 
	Px2x3 <- Px2*Px3; # If they are independet. 
	
	lm.x <- lm(pheno_local~Px1+Px2+Px3+Px1x2+Px1x2+Px1x3+Px2x3); 
 
	rss <- summary(lm.x)$sigma;
	sigma2 <- rss^2;
	f <-  -log(abs(-sigma2+variance_local)/variance_local)# is this correct? do we need dimention?

	return(f);
}
##################################################################################################################
##################################################################################################################
#3 - Dim, intercross

ObjectiveF3dInt <- function(x,Problem.data){ 
	# we need P(AA), P(AB), P(BB)
	probAA_X1_local <- Problem.data$prob$PAA[[1]] 
	probAA_X2_local <- Problem.data$prob$PAA[[2]] 
	probAA_X3_local <- Problem.data$prob$PAA[[3]] 
	
	probAB_X1_local <- Problem.data$prob$PAA[[1]] 
	probAB_X2_local <- Problem.data$prob$PAB[[2]] 
	probAB_X3_local <- Problem.data$prob$PAB[[3]] 
	
	probBB_X1_local <- Problem.data$prob$PBB[[1]] 
	probBB_X2_local <- Problem.data$prob$PBB[[2]] 
	probBB_X3_local <- Problem.data$prob$PBB[[3]] 
	
	pheno_local   <- Problem.data$pheno
	
	mesh1_local   <- Problem.data$mesh[[1]]
	mesh2_local   <- Problem.data$mesh[[2]]
	mesh3_local   <- Problem.data$mesh[[3]]
	
	variance_local<- Problem.data$variance 
	
	step1		  <- Problem.data$step[1]
    step2		  <- Problem.data$step[2]		
	step3		  <- Problem.data$step[3]		
	
	x1 <- x[1];
	x2 <- x[2];
	x3 <- x[3];

	if (x[1] == mesh1_local[1]){
		location1 <- 1
	} else
		location1 <- (x1-mesh1_local[1])/step1 + 1
	
	if (x[2] == mesh2_local[1]){
		location2 <- 1
	} else
		location2 <- (x2-mesh2_local[1])/step2 + 1
		
	if (x[3] == mesh3_local[1]){
		location3 <- 1
	} else
		location3 <- (x3-mesh3_local[1])/step3 + 1

	PAAx1 <- probAA_X1_local[, as.integer(as.character(location1))]; # prob for AA
	PAAx2 <- probAA_X2_local[, as.integer(as.character(location2))]; # prob for BB
	PAAx3 <- probAA_X3_local[, as.integer(as.character(location3))]; # prob for CC
	
	PABx1 <- probAB_X1_local[, as.integer(as.character(location1))]; # prob for AA
	PABx2 <- probAB_X2_local[, as.integer(as.character(location2))]; # prob for BB
	PABx3 <- probAB_X3_local[, as.integer(as.character(location3))]; # prob for CC
	
	PBBx1 <- probBB_X1_local[, as.integer(as.character(location1))]; # prob for AA
	PBBx2 <- probBB_X2_local[, as.integer(as.character(location2))]; # prob for BB
	PBBx3 <- probBB_X3_local[, as.integer(as.character(location3))]; # prob for CC
	
	w1a <- PAAx1-PBBx1 # addative effect for x1
	w1d <- 2*PAAx1;		#dominant effect for x1
	
	w2a <- PAAx2-PBBx2
	w2d <- 2*PAAx2;
	
	w3a <- PAAx2-PBBx2
	w3d <- 2*PAAx2;
	
	#interactions:
	w12a <- w1a*w2a
	w13a <- w1a*w3a
	w23a <- w2a*w3a
	
	w12d <- w1d*w2d
	w13d <- w1d*w3d
	w23d <- w2d*w3d
	
	lm.x <- lm(pheno_local~w1a+w1d+w2a+w2d+w3a+w3d+w12a+w13a+w23a+w12d+w13d+w23d); # (Change here), choose which effects you want
 
	rss <- summary(lm.x)$sigma;
	sigma2 <- rss^2;
	f <-  -log(abs(-sigma2+variance_local)/variance_local)# is this correct? do we need dimention?

	return(f);
}
##################################################################################################################
