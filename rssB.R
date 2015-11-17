##################################################################################################################
#------------------------------------------------------------------#
# Function   :			                                           # 
# Written by :  Behrang Mahjani                                    #
# Created on :  3/26/2015                                          #
# Purpose    :												       #
#------------------------------------------------------------------#

# These functions should be improved
# Thesis functions are only for 2d

 rss <- function(n,x,fmin)
{
	# distance after x  !!! be carefull what x means 
	p <- 1-(0.5 + 0.5*(exp(1))^(-2*x/100))

	n0_0 <- n/2;
	n1_0 <- n/2;	

	mu_0 <- 0;
	a_0  <-	sqrt((1.0 / (1.0 - exp(-fmin)) - 1)) * 2;
	RSSt <- (a_0 * 0.5 * a_0* 0.5 + 1) * (n);

	mp01 <- floor(n*p/2+8*sqrt(n*p*(1-p)/2))
	mp10 <- floor(n*p/2+8*sqrt(n*p*(1-p)/2))
	
	l<- floor(RSSt/10);
	for(i in 1:l)
	{	
		rssx_grid	   <-seq(RSSt-i*l,(RSSt-(i-1)*l),by=0.1)	
		quantile_index <- min(which(FRSS(rssx_grid,mp01,n,a_0,mu_0,RSSt,p)>0.999999))
		
		if (quantile_index != Inf) {
			f <- -log(abs(-rssx_grid[quantile_index]+RSSt)/RSSt)
			return(f)
		}
	}
    
    return(Inf)# just for safty *
}

##################################################################################################################

 FRSS <- function(rssx,mp01,n,a_0,mu_0,RSSt,p)
 {
	F1 <- 0;
	F2 <- 0;
	
	mp10 <- mp01;
	n0_0 <- n/2;
	n1_0 <- n/2;
	
	for (m01 in 0:mp01){
		F1 <- 0;
		for (m10 in 0:mp10){
			n0_x <- n0_0 + m10 - m01
			n1_x <- n1_0 + m01 - m10
			zbar_x <- n1_x/n

			a2  <- n0_x * zbar_x^2 +n1_x*(1-zbar_x)^2

		    a12 <- (1-zbar_x)*a_0*n1_0
			
			mu_a11	<- mu_0*(m01-m10)-a_0*m10
			s0		<- 1; 
			s1		<- 1;
			sigma_a11 <- sqrt(s0*m01*(1-(m01-1)/(n0_0-1))+s1*m10*(1-(m10-1)/(n1_0-1)))
			FCDF	  <- 1 - pnorm((sqrt((RSSt-rssx)*a2)-a12-mu_a11)/sigma_a11) + pnorm((-sqrt((RSSt-rssx)*a2)-a12-mu_a11)/sigma_a11)
			F1		  <- F1+ dbinom(m10,n1_0, p)*FCDF
		}
		F2 <- F2 + dbinom(m01,n0_0, p)*F1
	}
	return(F2);
 }
##################################################################################################################

