
thetaInit <- function(pred, mult){
	

	
	model <- lm(resp~pred[, -1])
	coeff <- cbind(model$coefficients, model$coefficients)
	coeff[1,] <- coeff[1,] * c(0.99, 1.01)
	sigma <- c(summary(model)$sigma,summary(model)$sigma)
	
	beta <- beta.init(mult)
	

	val <- c(coeff[1,1], coeff[1,2], coeff[-1,1], sigma[1], beta[,1], beta[,2])
	names(val) <- c("alpha0_0","alpha1_0",paste("alpha", 1: ncol(pred[,-1]), sep = "_"), "sigma", 
	paste("beta0", 1:ncol(mult), sep="_"), paste("beta1", 1:ncol(mult), sep="_"))
	
	return(val)
}


