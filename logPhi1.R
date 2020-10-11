logPhi1 <-
function(x,Phi.alpha0,Phi.alpha1,Phi.beta){
  u0=x[1]; u1=x[2]
  (Phi.alpha1-Phi.beta*u0)-log(exp(Phi.alpha0-Phi.beta*u1)+exp(Phi.alpha1-Phi.beta*u0))
}

