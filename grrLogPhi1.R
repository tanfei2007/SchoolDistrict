grrLogPhi1 <-
function(x,Phi.alpha0,Phi.alpha1,Phi.beta){
    u0=x[1]; u1=x[2];
    alpha0Prime= - exp(Phi.alpha0-Phi.beta*u1)/( exp(Phi.alpha0-Phi.beta*u1)+exp(Phi.alpha1-Phi.beta*u0) )
    alpha1Prime=1- exp(Phi.alpha1-Phi.beta*u0)/( exp(Phi.alpha0-Phi.beta*u1)+exp(Phi.alpha1-Phi.beta*u0) )
    betaPrime= -u0+ ( u1*exp(Phi.alpha0-Phi.beta*u1)+u0*exp(Phi.alpha1-Phi.beta*u0) )/( exp(Phi.alpha0-Phi.beta*u1)+exp(Phi.alpha1-Phi.beta*u0) )
    c(alpha0Prime, alpha1Prime, betaPrime)
}

