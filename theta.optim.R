

#variance sigma
get.sigma <- function(resp, # reponse/dependent variable vector of n * 1
	 					  u, #centers vector of n * k
						  gamma # posterior the matrix of n(instances) * k(components)
						  ){
						  	
							  den.sum <- length(resp)
							  
							  #browser()
							  
							  num <-  gamma * (resp - u)^2 
							  #resp <- c(1,2,3)
							  #u <- cbind(c(0.9,1.8, 2.8), c(0.8,1.9, 2.7))
							  #(resp - u)^2
						      #[,1] [,2]
						 	  #[1,] 0.01 0.04
						 	  #[2,] 0.04 0.01
						 	  #[3,] 0.04 0.09
							  #num
				              #V1          V2
				 			  #[1,] 0.00453240 0.021870398
				 			  #[2,] 0.01712745 0.005718137
				 			  #[3,] 0.01822326 0.048997663
							  
							  num.sum <- sum(num)
							  
							  
							  val <- num.sum/den.sum
							  
							  val <- sqrt(val)
							  val <- c(val, val)
							  
							  return(val)
						  }
	
					 
#coeff
get.coeff <- function(resp, # reponse/dependent variable vector of n * 1
					  pred, # predictor variable matrix of n * (p+1)
					  gamma, # posterior the matrix of n(instances) * k(components)
					  coeff # regression coefficients matrix of (p+1) * k
					  ){
	
						  #browser()
   
    					  pred.itcp <- pred[, 1]
						  #head(pred.itcp, 2)
						  #1 1
						  coeff.itcp <- coeff[1, ]
	
	
						  pred.other <- pred[, -1]
						  coeff.other <- coeff[-1,]
						  
						 
						  ##estimate Intercept Coefficients for two components
	
						  den.sum.itcp <- colSums(gamma)
						  num.other <- gen.center(pred.other, coeff.other)
						  #   [,1] [,2]
						  #[1,]    7 11.4
						  #[2,]    8 13.2
						  #[3,]    9 15.0
						  num.sum.itcp <- gamma * (resp - num.other)
						  num.sum.itcp <- colSums(num.sum.itcp)
						  beta.itcp <- num.sum.itcp/den.sum.itcp
	


    					  ##estimate other coefficients for two components,
						  #note here they are literally same
	
						  coeff.other.common <- coeff.other[,1]
						  for(i in 1:length(coeff.other.common)){
							  
							  den.sum.other <- sum(pred.other[,i]^2)
							  #pred <- matrix((1:12), 3, 4)
							  #     [,1] [,2] [,3] [,4]
							  #[1,]    1    4    7   10
							  #[2,]    2    5    8   11
							  #[3,]    3    6    9   12
	

	    					  #select the the random column (all columns are same)
						  
							  coeff.ss <- gen.coeff.pred.prod.sum(pred.other, coeff.other.common, i)
	    					  #[,1] [,2] [,3] 
							  #[1,]  6.9  6.2  4.9  
							  #[2,]  7.8  7.0  5.6  
							  #[3,]  8.7  7.8  6.3  
	
							  coeff.itcp.ss <- as.vector(gamma %*% beta.itcp) * pred.other[, i]
	
							  coeff.ss <- coeff.ss + coeff.itcp.ss
   
	    					  #browser()	   
		
							  num.sum.other <- sum(resp * pred.other[, i] - coeff.ss) 

							  beta.other <- num.sum.other/den.sum.other
							  coeff.other.common[i] <- beta.other
						  }
						  
						  
	
						  coeff.other.common <- cbind(coeff.other.common, coeff.other.common)
	
						  val <- rbind(beta.itcp, coeff.other.common)


						  return(val)
}				  

#beta
get.beta <- function(w, # n(number of instances) * V(vocabulary size)
	 				gamma # n(number of instances) * k(number of component)
					){
					
					    epis <- 10^(-10)
						#denominator vecotr of k 
						den <- gamma * rowSums(w)
						den.sum <- colSums(den)
						
						#numerator V * k
						num.sum <- t(w) %*% gamma
					    
						# k * V
						val <- (t(num.sum) + epis)/(den.sum + epis)
						
						# V * k
						val <- t(val)
						
						return(val)
						
						#gamma <- cbind(c(0.3, 0.4, 0.5), c(0.7, 0.6, 0.5))
						#w <- rbind(c(1,2), c(1,1), c(1,0))
}




theta.optim <- function(resp, pred, mult, ini.Theta, curr.X0){
	
	
	##convert theta to the required format
	coeff0 <- ini.Theta[c(1, (2+1):(2+ncol(pred[,-1])))]
	sigma0 <- ini.Theta[2 + ncol(pred[,-1]) + 1]
	beta0 <- ini.Theta[(2 + ncol(pred[,-1]) + 1 + 1): (2 + ncol(pred[,-1]) + 1 + ncol(mult))]
	
	coeff1 <- ini.Theta[2:(2+ncol(pred[,-1]))]
	sigma1 <- ini.Theta[2 + ncol(pred[,-1]) + 1]
	beta1 <- ini.Theta[(2 + ncol(pred[,-1]) + 1 + ncol(mult) + 1):(2 + ncol(pred[,-1]) + 1 + ncol(mult) + ncol(mult))]
	
	coeff <- cbind(coeff0, coeff1)
	sigma <- c(sigma0, sigma1)
	beta <- cbind(beta0, beta1)
	
	gamma <- cbind(1 - curr.X0, curr.X0)
	
	
	## update the new parameter
	coeff <- get.coeff(resp, pred, gamma, coeff)
	
	u <- gen.center(pred, coeff)
	sigma <- get.sigma(resp, u, gamma)
	
	beta <- get.beta(mult, gamma)
	
	## convert the obtained format to the Theta-based format
	val <- c(coeff[1,1], coeff[1,2], coeff[-1,1], sigma[1], beta[,1], beta[,2])
	names(val) <- c("alpha0_0","alpha1_0",paste("alpha", 1:ncol(pred[,-1]), sep = "_"), "sigma", paste("beta0", 1:ncol(mult), sep="_"), paste("beta1", 1:ncol(mult), sep="_"))
	
	return(val)
}

gs.theta.optim <- function(resp, pred, mult, ini.Theta, curr.X0){
	
	
	##convert theta to the required format
	coeff0 <- ini.Theta[c(1, (2+1):(2+ncol(pred[,-1])))]
	sigma0 <- ini.Theta[2 + ncol(pred[,-1]) + 1]
	beta0 <- ini.Theta[(2 + ncol(pred[,-1]) + 1 + 1): (2 + ncol(pred[,-1]) + 1 + ncol(mult))]
	
	coeff1 <- ini.Theta[2:(2+ncol(pred[,-1]))]
	sigma1 <- ini.Theta[2 + ncol(pred[,-1]) + 1]
	beta1 <- ini.Theta[(2 + ncol(pred[,-1]) + 1 + ncol(mult) + 1):(2 + ncol(pred[,-1]) + 1 + ncol(mult) + ncol(mult))]
	
	coeff <- cbind(coeff0, coeff1)
	sigma <- c(sigma0, sigma1)
	beta <- cbind(beta0, beta1)
	
	gamma <- cbind(1 - curr.X0, curr.X0)
	
	
	## update the new parameter
	coeff <- get.coeff(resp, pred, gamma, coeff)
	
	u <- gen.center(pred, coeff)
	sigma <- get.sigma(resp, u, gamma)
	
	#beta <- get.beta(mult, gamma)
	
	## convert the obtained format to the Theta-based format
	val <- c(coeff[1,1], coeff[1,2], coeff[-1,1], sigma[1], beta[,1], beta[,2])
	names(val) <- c("alpha0_0","alpha1_0",paste("alpha", 1:ncol(pred[,-1]), sep = "_"), "sigma", paste("beta0", 1:ncol(mult), sep="_"), paste("beta1", 1:ncol(mult), sep="_"))
	
	return(val)
}

mult.theta.optim <- function(resp, pred, mult, ini.Theta, curr.X0){
	
	
	##convert theta to the required format
	coeff0 <- ini.Theta[c(1, (2+1):(2+ncol(pred[,-1])))]
	sigma0 <- ini.Theta[2 + ncol(pred[,-1]) + 1]
	beta0 <- ini.Theta[(2 + ncol(pred[,-1]) + 1 + 1): (2 + ncol(pred[,-1]) + 1 + ncol(mult))]
	
	coeff1 <- ini.Theta[2:(2+ncol(pred[,-1]))]
	sigma1 <- ini.Theta[2 + ncol(pred[,-1]) + 1]
	beta1 <- ini.Theta[(2 + ncol(pred[,-1]) + 1 + ncol(mult) + 1):(2 + ncol(pred[,-1]) + 1 + ncol(mult) + ncol(mult))]
	
	coeff <- cbind(coeff0, coeff1)
	sigma <- c(sigma0, sigma1)
	beta <- cbind(beta0, beta1)
	
	gamma <- cbind(1 - curr.X0, curr.X0)
	

	
	beta <- get.beta(mult, gamma)
	
	## convert the obtained format to the Theta-based format
	val <- c(coeff[1,1], coeff[1,2], coeff[-1,1], sigma[1], beta[,1], beta[,2])
	names(val) <- c("alpha0_0","alpha1_0",paste("alpha", 1:ncol(pred[,-1]), sep = "_"), "sigma", paste("beta0", 1:ncol(mult), sep="_"), paste("beta1", 1:ncol(mult), sep="_"))
	
	return(val)
}






