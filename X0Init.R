
X0Init <- function(currTheta, resp, pred, mult){
	

	coeff0 <- currTheta[c(1, (2+1):(2+ncol(pred[,-1])))]
	sigma0 <- currTheta[2 + ncol(pred[,-1]) + 1]
	beta0 <- currTheta[(2 + ncol(pred[,-1]) + 1 + 1): (2 + ncol(pred[,-1]) + 1 + ncol(mult))]
    like.theta0 <- gsmult.mrf.log.likelihood(resp = resp, coeff = coeff0, pred = pred, sigma = sigma0, w = mult, beta = beta0)
	
	coeff1 <- currTheta[2:(2+ncol(pred[,-1]))]
	sigma1 <- currTheta[2 + ncol(pred[,-1]) + 1]
	beta1 <- currTheta[(2 + ncol(pred[,-1]) + 1 + ncol(mult) + 1):(2 + ncol(pred[,-1]) + 1 + ncol(mult) + ncol(mult))]
    like.theta1 <- gsmult.mrf.log.likelihood(resp = resp, coeff = coeff1, pred = pred, sigma = sigma1, w = mult, beta = beta0)
	
	#browser()
	
	val <- as.numeric((like.theta1-like.theta0)>0)
	
	
	
	return(val)
}


