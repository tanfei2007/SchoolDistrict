

library(matrixStats)

#calculate the gaussian probability in log space
gs.prob <- function(x, u, sigma){
	val <- sapply(x, dnorm, mean = u, sd = sigma, log = T)
	return(val)
}
#gp(c(3,4), c(1,0.5), c(0.2,0.3))
#          [,1]       [,2]      [,3]
#[1,] -49.30950 -111.80950 -199.3095
#[2,] -34.43719  -67.77052 -112.2150


#calculate the multinomial probability
mult.prob <- function(w, beta){
 val <- dmultinom(w, prob = beta, log = T)
 return(val)  # the log of probability 
}


#calculate the log sum
# log(sum(exp(logx)))
log.sum <- function(x){
	val <- logSumExp(x)
	return(val)
}
#x <- c(3,4)
#log.sum(log(x))
#1.94591
#log(sum(x))
#1.94591


#generate the center
gen.center <- function(pred, coeff){
	val <- pred %*% coeff
	return(val)
}
#coeff <- cbind(c(0.1,0.2, 0.3, 0.4), c(0.3, 0.4, 0.5, 0.6))
#     [,1] [,2]
#[1,]  0.1  0.3
#[2,]  0.2  0.4
#[3,]  0.3  0.5
#[4,]  0.4  0.6

#pred <- matrix((1:12), 3, 4)
#     [,1] [,2] [,3] [,4]
#[1,]    1    4    7   10
#[2,]    2    5    8   11
#[3,]    3    6    9   12

#gen.center(pred, coeff)
#     [,1] [,2]
#[1,]    7 11.4
#[2,]    8 13.2
#[3,]    9 15.0


#calculate the sufficient statistics coeff.pred.sum
gen.coeff.pred.prod.sum <- function(pred, coeff.part, z){
	center <- pred %*% coeff.part
	center <- as.vector(center)
	#coeff.part
	#[1] 0.1 0.2 0.3 0.4
	#center
	#[1] 7 8 9
	
	center.part <- t(t(pred) * coeff.part)
	#center.part
    # [,1] [,2] [,3] [,4]
	#[1,]  0.1  0.8  2.1  4.0
	#[2,]  0.2  1.0  2.4  4.4
	#[3,]  0.3  1.2  2.7  4.8
	val <- center - center.part
	#val
    #[,1] [,2] [,3] [,4]
	#[1,]  6.9  6.2  4.9  3.0
	#[2,]  7.8  7.0  5.6  3.6
	#[3,]  8.7  7.8  6.3  4.2
	
	val <- pred * val
	
	val <- val[,z]
	
	return(val)
}
#


#calculate the posterior for regression with gaussian error
gs.gamma.post <- function(resp, # reponse/dependent variable vector of n * 1
						u, #centers vector of n * 2
						sigma, # standard deviation here is vector of 2 * 1
						Pi # the prior of hidden variable here is vector of 2 * 1
						){
								
							u.frame <- as.data.frame(u)
							log_prob <- mapply(function(x,y){
								mapply(gs.prob, resp, x, y)
							}, u.frame, sigma)
							
							#resp <- c(1,2,3)
							#u <- cbind(c(0.9,1.8, 2.8), c(0.8,1.9, 2.7))
							#u.frame
						    #V1  V2
						 	#0.9 0.8
						 	#1.8 1.9
							#2.8 2.7
							#sigma <- c(0.5, 0.6)
							#log_prob
			                #V1         V2
			   			 	#[1,] -0.2457914 -0.4636685
			   			 	#[2,] -0.3057914 -0.4220018
			   			 	#[3,] -0.3057914 -0.5331129
							
							log_prob_Pi <- log(Pi) + t(log_prob)
							log_prob_Pi <- t(log_prob_Pi)
							
							#Pi <- c(0.4, 0.6)
							#log_prob_Pi
				            #V1         V2
							#[1,] -1.162082 -0.9744941
							#[2,] -1.222082 -0.9328274
							#[3,] -1.222082 -1.0439385
							
							log_prob_Pi_sum <- apply(log_prob_Pi, 1, log.sum)
							#log_prob_Pi_sum
							#[1] -0.3707487 -0.3738853 -0.4359015
							
							log_gamma <- log_prob_Pi - log_prob_Pi_sum
							#log_gamma
			                #V1         V2
			   			 	#[1,] -0.7913334 -0.6037454
			   			 	#[2,] -0.8481968 -0.5589421
			   			 	#[3,] -0.7861806 -0.6080371
							gamma <- exp(log_gamma)
							#gamma
							#            V1        V2
							#[1,] 0.4532400 0.5467600
							#[2,] 0.4281863 0.5718137
							#[3,] 0.4555815 0.5444185
							return(gamma)
						}
	
#calculate the likelihood for gaussian regression						
gs.log.likelihood <- function(resp, # reponse/dependent variable vector of n * 1
						u, #centers vector of n * 2
						sigma, # standard deviation here is vector of 2 * 1
						Pi # the prior of hidden variable here is vector of 2 * 1
						){
							
							u.frame <- as.data.frame(u)
							log_prob <- mapply(function(x,y){
								mapply(gs.prob, resp, x, y)
							}, u.frame, sigma)
							
							#resp <- c(1,2,3)
							#u <- cbind(c(0.9,1.8, 2.8), c(0.8,1.9, 2.7))
							#u.frame
						    #V1  V2
						 	#0.9 0.8
						 	#1.8 1.9
							#2.8 2.7
							#sigma <- c(0.5, 0.6)
							#log_prob
			                #V1         V2
			   			 	#[1,] -0.2457914 -0.4636685
			   			 	#[2,] -0.3057914 -0.4220018
			   			 	#[3,] -0.3057914 -0.5331129
							
							log_prob_Pi <- log(Pi) + t(log_prob)
							log_prob_Pi <- t(log_prob_Pi)
							
							#Pi <- c(0.4, 0.6)
							#log_prob_Pi
				            #V1         V2
							#[1,] -1.162082 -0.9744941
							#[2,] -1.222082 -0.9328274
							#[3,] -1.222082 -1.0439385
							
							log_prob_Pi_sum <- apply(log_prob_Pi, 1, log.sum)
							#log_prob_Pi_sum
							#[1] -0.3707487 -0.3738853 -0.4359015
							
							log_likelihood <- sum(log_prob_Pi_sum)
							#log_likelihood
							#[1] 0.7055473
							return(log_likelihood)
						}
						

					
#calculate the posterior for multinomial distribution						
mult.gamma.post <- function(w, # the matrix of n (the number of instances) * V
							beta, # the matrix of  V (the vocabulary size) * k(the number of component)
							Pi #priori
							){
								log_prob <- apply(beta, 2, function(x){
									val <- apply(w, 1, mult.prob, beta = x)
								})
								#beta <- cbind(c(0.3,0.7), c(0.4,0.6))
								#w <- rbind(c(1,2), c(1,1), c(1,0))
							    #exp(log_prob)
							    #      [,1]  [,2]
								#[1,] 0.441 0.432
								#[2,] 0.420 0.480
								#[3,] 0.300 0.400
								
								log_prob_Pi <- log(Pi) + t(log_prob)
								log_prob_Pi <- t(log_prob_Pi)
								 
								log_prob_Pi_sum <- apply(log_prob_Pi, 1, log.sum)
							
								log_gamma <- log_prob_Pi - log_prob_Pi_sum
						
								gamma <- exp(log_gamma)
							
								return(gamma)
}			


#calculate the likelihood for gaussian regression		
mult.log.likelihood <- function(w, beta, Pi){
	log_prob <- apply(beta, 2, function(x){
		val <- apply(w, 1, mult.prob, beta = x)
	})
	#beta <- cbind(c(0.3,0.7), c(0.4,0.6))
	#w <- rbind(c(1,2), c(1,1), c(1,0))
    #exp(log_prob)
    #      [,1]  [,2]
	#[1,] 0.441 0.432
	#[2,] 0.420 0.480
	#[3,] 0.300 0.400
	
	log_prob_Pi <- log(Pi) + t(log_prob)
	log_prob_Pi <- t(log_prob_Pi)
	 
	log_prob_Pi_sum <- apply(log_prob_Pi, 1, log.sum)
	
	log_likelihood <- sum(log_prob_Pi_sum)
	
	return(log_likelihood)	
}


#calculate the posterior for regression with gaussian error and multinomial distribution
gsmult.gamma.post <- function(resp, # reponse/dependent variable vector of n * 1
						u, #centers vector of n * 2
						sigma, # standard deviation here is vector of 2 * 1
						w, # the matrix of n (the number of instances) * V
						beta, # the matrix of  V (the vocabulary size) * k(the number of component)
						Pi # the prior of hidden variable here is vector of 2 * 1
						){
								
								
							#gaussian part	
							u.frame <- as.data.frame(u)
							log_prob_gs <- mapply(function(x,y){
								mapply(gs.prob, resp, x, y)
							}, u.frame, sigma)
							
							#multinomial part
							log_prob_mult <- apply(beta, 2, function(x){
								val <- apply(w, 1, mult.prob, beta = x)
							})
							
							
							log_prob <- log_prob_gs + log_prob_mult
							
							log_prob_Pi <- log(Pi) + t(log_prob)
							log_prob_Pi <- t(log_prob_Pi)
							
							#Pi <- c(0.4, 0.6)
							#log_prob_Pi
				            #V1         V2
							#[1,] -1.162082 -0.9744941
							#[2,] -1.222082 -0.9328274
							#[3,] -1.222082 -1.0439385
							
							log_prob_Pi_sum <- apply(log_prob_Pi, 1, log.sum)
							#log_prob_Pi_sum
							#[1] -0.3707487 -0.3738853 -0.4359015
							
							log_gamma <- log_prob_Pi - log_prob_Pi_sum
							#log_gamma
			                #V1         V2
			   			 	#[1,] -0.7913334 -0.6037454
			   			 	#[2,] -0.8481968 -0.5589421
			   			 	#[3,] -0.7861806 -0.6080371
							gamma <- exp(log_gamma)
							#gamma
							#            V1        V2
							#[1,] 0.4532400 0.5467600
							#[2,] 0.4281863 0.5718137
							#[3,] 0.4555815 0.5444185
							return(gamma)
						}



#calculate the log likelihood for regression with gaussian error and multinomial distribution
gsmult.log.likelihood <- function(resp, # reponse/dependent variable vector of n * 1
							  u, #centers vector of n * 2
							  sigma, # standard deviation here is vector of 2 * 1
							  w, # the matrix of n (the number of instances) * V
							  beta, # the matrix of  V (the vocabulary size) * k(the number of component)
							  Pi # the prior of hidden variable here is vector of 2 * 1
						){	
							#gaussian part	
							u.frame <- as.data.frame(u)
							log_prob_gs <- mapply(function(x,y){
									mapply(gs.prob, resp, x, y)
							}, u.frame, sigma)
							
							#multinomial part
							log_prob_mult <- apply(beta, 2, function(x){
									val <- apply(w, 1, mult.prob, beta = x)
							})
							
							log_prob <- log_prob_gs + log_prob_mult
							
							log_prob_Pi <- log(Pi) + t(log_prob)
							log_prob_Pi <- t(log_prob_Pi)
						
							log_prob_Pi_sum <- apply(log_prob_Pi, 1, log.sum)
													
							log_likelihood <- sum(log_prob_Pi_sum)
	
							return(log_likelihood)
					}


#calculate the loglikelihood for regression with gaussian error and multinomial distribution
gsmult.mrf.log.likelihood <- function(resp, # reponse/dependent variable vector of n * 1
							  coeff, #centers vector of p
							  pred, # n * p
							  sigma, # standard deviation here is scalar value
							  w, # the matrix of n (the number of instances) * V
							  beta # the vector of  V (the vocabulary size) 
							  
						){	
							
							
							#gaussian part	
							u <- gen.center(pred, coeff)
							
							log_prob_gs <- mapply(gs.prob, resp, u, sigma)
						
		
							#multinomial part
							log_prob_mult <- apply(w, 1, mult.prob, beta = beta)
							
							log_prob <- log_prob_gs + log_prob_mult
							
							#browser()
															
							log_likelihood <- log_prob

							return(log_likelihood)
					}
					



#calculate the loglikelihood for regression with gaussian error and multinomial distribution
gs.mrf.log.likelihood <- function(resp, # reponse/dependent variable vector of n * 1
							  coeff, #centers vector of p
							  pred, # n * p
							  sigma, # standard deviation here is scalar value
							  w, # the matrix of n (the number of instances) * V
							  beta # the vector of  V (the vocabulary size) 
							  
						){	
							
							#gaussian part	
							u <- gen.center(pred, coeff)
							
							log_prob_gs <- mapply(gs.prob, resp, u, sigma)
						
							
							#browser()
															
							log_likelihood <- log_prob_gs

							return(log_likelihood)
					}


#calculate the loglikelihood for regression with gaussian error and multinomial distribution
mult.mrf.log.likelihood <- function(resp, # reponse/dependent variable vector of n * 1
							  coeff, #centers vector of p
							  pred, # n * p
							  sigma, # standard deviation here is scalar value
							  w, # the matrix of n (the number of instances) * V
							  beta # the vector of  V (the vocabulary size) 
							  
						){	
							
							#multinomial part
							log_prob_mult <- apply(w, 1, mult.prob, beta = beta)
							
							
							#browser()							
							log_likelihood <- log_prob_mult

							return(log_likelihood)
					}

	
beta.init <- function(w){
	beta1 <- colSums(w)/sum(w)
	beta2 <- beta1 + runif(length(beta1), 0, min(beta1))
	beta2 <- beta2/sum(beta2)
	beta <- cbind(beta1, beta2)	
}




