grrLoglikPhi <-
function(phi,u,X0) {
    p.alpha0=phi[1]; p.alpha1=phi[2]; p.beta=phi[3];
    lphi0=apply(matrix(u[(X0==0),],ncol=ncol(u)),1, grrLogPhi0, Phi.alpha0=p.alpha0, Phi.alpha1=p.alpha1, Phi.beta=p.beta)
    lphi1=apply(matrix(u[(X0==1),],ncol=ncol(u)),1, grrLogPhi1, Phi.alpha0=p.alpha0, Phi.alpha1=p.alpha1, Phi.beta=p.beta)

    if (length(lphi0)>0) {
    	ans<-rowSums(lphi0)
    }
    if (length(lphi1)>0) {
    	ans<-rowSums(lphi1)
    }
    if ((length(lphi0)>0) && (length(lphi1)>0)){
    	ans<-rowSums(lphi0)+rowSums(lphi1)
    }
    ans
}

