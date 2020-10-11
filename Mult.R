rm(list = ls())
source("Utility.R")

w <- read.csv("../Data/mult.data.csv", string = F)
all <- read.csv("../Data/gs.data.csv", string = F)
ground.truth <- all$Key

w <- as.matrix(w)


init.Pi <- c(0.5, 0.5)

init.beta <- beta.init(w)

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

###EM Process
#browser()
Pi <- init.Pi
beta <- init.beta

loglike <- mult.log.likelihood(w, beta, Pi)
eta = 10^(-10)

maxItr = 200
count <- 0

mult.start <- proc.time()

repeat{
	
	
	#browser()
	##E step
	gamma <- mult.gamma.post(w, beta, Pi)

	##M step
	Pi <- get.Pi(gamma)
	beta <- get.beta(w, gamma)
	
	
	old.loglike <- loglike
	loglike <- mult.log.likelihood(w, beta, Pi)
	
	#check convergence
	bound <- (loglike - old.loglike)/abs(old.loglike)
	#print(bound)
	#print(beta)
	count <- count + 1
	if(bound < eta || count > maxItr){
		break
	}
}

mult.time <- proc.time() - mult.start
browser()

save(gamma, beta, ground.truth, file = "gamma.mult.RData")

rslt <- cbind(ground.truth, gamma, w)
write.csv(rslt, "rslt.mult.csv", row.names = F, quote = T)





