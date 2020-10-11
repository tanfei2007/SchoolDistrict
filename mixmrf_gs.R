
mxmrf.gs <-
function(resp, pred, mult, nb, maxIter=200) {

  ###initilize the parameters
  ini.Theta <- thetaInit(pred, mult)
  ini.Phi=c(1,1,1)
  Theta.trace=ini.Theta;  Phi.trace=ini.Phi
  Theta.trace=cbind(Theta.trace,ini.Theta);  Phi.trace=cbind(Phi.trace,ini.Phi)
  

  curr.X0 <- X0Init(ini.Theta, resp, pred, mult)
  

  
  
  for(run in 1:maxIter){
	  
	#print(run)

	 # browser()
    u=t(apply(nb,1, function(x) c(sum((curr.X0[x])==0,na.rm=T),sum((curr.X0[x])==1,na.rm=T))))
	#u <- u/rowSums(u)
    
	##derivative
	newTheta <- gs.theta.optim(resp, pred, mult, ini.Theta, curr.X0)
	
	newPhi = optim(par = ini.Phi,fn = loglikPhi, gr = grrLoglikPhi, X0=curr.X0, u=u, method = c("BFGS"), control = list(fnscale=-1))

    #parameter update
	currTheta=newTheta
    #currTheta=newTheta$par   
	currPhi=newPhi$par
    
	#check for optimization
	ini.Theta=currTheta
    #if(newTheta$convergence==0) {
    #   ini.Theta=currTheta
    #}
    if(newPhi$convergence==0) {
       ini.Phi=currPhi
    }
    
    Theta.trace=cbind(Theta.trace,currTheta); Phi.trace=cbind(Phi.trace, currPhi)
    
	
    #multi regression
	coeff0 <- currTheta[c(1, (2+1):(2+ncol(pred[,-1])))]
	sigma0 <- currTheta[2 + ncol(pred[,-1]) + 1]
	beta0 <- currTheta[(2 + ncol(pred[,-1]) + 1 + 1): (2 + ncol(pred[,-1]) + 1 + ncol(mult))]
    like.theta0 <- gs.mrf.log.likelihood(resp = resp, coeff = coeff0, pred = pred, sigma = sigma0, w = mult, beta = beta0)
	
	coeff1 <- currTheta[2:(2+ncol(pred[,-1]))]
	sigma1 <- currTheta[2 + ncol(pred[,-1]) + 1]
	beta1 <- currTheta[(2 + ncol(pred[,-1]) + 1 + ncol(mult) + 1):(2 + ncol(pred[,-1]) + 1 + ncol(mult) + ncol(mult))]
    like.theta1 <- gs.mrf.log.likelihood(resp = resp, coeff = coeff1, pred = pred, sigma = sigma1, w = mult, beta = beta1)
    
	#mrf
    like.phi0 <- apply(u,1, logPhi0, Phi.alpha0=currPhi[1], Phi.alpha1=currPhi[2], Phi.beta=currPhi[3])
    like.phi1 <- apply(u,1, logPhi1, Phi.alpha0=currPhi[1], Phi.alpha1=currPhi[2], Phi.beta=currPhi[3]) 
	
	#browser()
    
	#total conditional prob
    prob0 <- like.theta0 + like.phi0;
    prob1 <- like.theta1 + like.phi1;
	
	#assigmnment of state 
    curr.X0 <- as.numeric((prob1-prob0)>0)

    #convergence criterion
    stop1 <- max(abs(Phi.trace[,ncol(Phi.trace)]-Phi.trace[,ncol(Phi.trace)-1])/Phi.trace[,ncol(Phi.trace)-1])
    stop1 <- min(stop1, max(abs(Phi.trace[,ncol(Phi.trace)]-Phi.trace[,ncol(Phi.trace)-2])/Phi.trace[,ncol(Phi.trace)-2]))
    
	stop2 <- max(abs(Theta.trace[,ncol(Theta.trace)]-Theta.trace[,ncol(Theta.trace)-1])/Theta.trace[,ncol(Theta.trace)-1])
    stop2 <-  min(stop2, max(abs(Theta.trace[,ncol(Theta.trace)]-Theta.trace[,ncol(Theta.trace)-2])/Theta.trace[,ncol(Theta.trace)-2]))
         
    #print(c(stop1, stop2))	
	#print(run) 
    if (stop1<0.00001 & stop2<0.00001 | all(curr.X0==0)) {
    	break;
    }
  }
  
  postProb0 <- 1/(1+exp(prob1-prob0))
 
  #Theta
  Theta.trace <- Theta.trace[,-1]
  colnames(Theta.trace)=paste("run",1:ncol(Theta.trace),sep="")
  rownames(Theta.trace) <- c("alpha0_0","alpha1_0",paste("alpha", 1:ncol(pred[,-1]), sep = "_"), "sigma", paste("beta0", 1:ncol(mult), sep="_"), paste("beta1", 1:ncol(mult), sep="_"))
  
  #Phi
  Phi.trace <- Phi.trace[,-1]
  colnames(Phi.trace)=paste("run",1:ncol(Phi.trace),sep="")
  rownames(Phi.trace)=c("gamma0","gamma1","beta")
  
  return(list(X0=curr.X0, postProb0=postProb0,Theta.trace=Theta.trace, Phi.trace=Phi.trace))

}

