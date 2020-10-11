rm(list = ls())
source("Utility.R")

#read gaussian data
all <- read.csv("../Data/gs.data.csv", string = F)


ground.truth <- all$Key

#read multinomial data
w <- read.csv("../Data/mult.data.csv", string = F)
w <- as.matrix(w)



ground.truth <- all$Key

#################################################################

all <- all[, -c(2, 3,5, 9,19)]

all <- cbind(price = all[,1], Intercept = rep(1, nrow(all)), all[,-1])




pred <- as.matrix(all[,-1])
resp <- all[,1]

pred <- as.matrix(all[,-1])
resp <- all[,1]




##initialize the parameters
#set coeff to generate the initial u
model <- lm(resp~pred[, -1])
init.coeff <- cbind(model$coefficients, model$coefficients)
init.coeff[1,] <- init.coeff[1,] * c(1.01, 0.99)



init.sigma <- c(summary(model)$sigma,summary(model)$sigma)

init.Pi <- c(0.5,0.5)


init.u <- gen.center(pred, init.coeff)



init.beta <- beta.init(w)

###################################
#Pi
get.Pi <- function(gamma # posterior the matrix of n(instances) * k(components)
	){
		val <- colSums(gamma)/(nrow(gamma))
		#gamma
		#            V1        V2
		#[1,] 0.4532400 0.5467600
		#[2,] 0.4281863 0.5718137
		#[3,] 0.4555815 0.5444185
		#val
		#       V1        V2 
		#0.4456693 0.5543307
		
		return(val)
	}
	
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


#########################
###EM Process

#browser()
Pi <- init.Pi
coeff <- init.coeff
u <- init.u
sigma <- init.sigma

beta <- init.beta

loglike <- gsmult.log.likelihood(resp, u, sigma, w, beta, Pi)
eta = 10^(-10)

#browser()
maxItr = 200
count <- 0

gs.mult.start <- proc.time()

repeat{
	
	
	#browser()
	##E step
	gamma <- gsmult.gamma.post(resp, u, sigma, w, beta, Pi)

	##M step
	Pi <- get.Pi(gamma)
	old.coeff <- coeff
	coeff <- get.coeff(resp, pred, gamma, old.coeff)
	
	old.u <- u
	u <- gen.center(pred, coeff)
	sigma <- get.sigma(resp, u, gamma)
	
	beta <- get.beta(w, gamma)
	
	
	old.loglike <- loglike
	loglike <- gsmult.log.likelihood(resp, u, sigma, w, beta, Pi)
	
	#check convergence
	bound <- (loglike - old.loglike)/abs(old.loglike)
	#print(bound)
	#print(coeff[1,])
	count <- count + 1
	#print(count)
	if(abs(bound) < eta || count > maxItr){
		break
	}
}

gs.mult.time <- proc.time() - gs.mult.start
browser()

save(gamma, beta, coeff, ground.truth, file = "gamma.gs.mult.RData")
rslt <- cbind(ground.truth, gamma, resp, pred, w)
write.csv(rslt, "rslt.gs.mult.csv", row.names = F, quote = T)






