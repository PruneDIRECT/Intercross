###########################################################################################################
 #------------------------------------------------------------------#
 # Function   :  CallObjFcn                                         #
 # Written by :  Dan Finkel                                         #   
 # Created on :  06/07/2004                                         #
 # Purpose    :  Evaluate ObjFcn at pointed specified               #
 # This function is modified by B. Mahjani on 3/26/2015             #
 #------------------------------------------------------------------#
 .CallObjFcn<-function(Problem,
					   Problem.data,
					   point.x,
					   a,b,...){
 
# point.x = vector of values for tuning parametr(s)  
# Scale variable back to original space
 a <- Problem.data$ends[,1] # corect for non-sym
 b <- Problem.data$ends[,2]

 point		<- abs(b - a)*point.x + a
 fcn_value  <- eval(parse(text=Problem$f))(point,Problem.data, ...)  
 	   
 return(data.frame("fcn_value"=fcn_value))
 }
##################################################################################################################