negloglike_egss_remle <- function(yt,tt,fguess_egss){
  
  #theta         <- fguess_egss[1];
  sigmasq       <- exp(fguess_egss[2]);
  tausq         <- exp(fguess_egss[3]);
  #xo            <- fguess_egss[4];
  q             <- length(yt) - 1;
  qp1           <- q+1;
  vx            <- matrix(0,qp1,qp1);
  for(i in 1:q){ 
    vx[(i+1):qp1,(i+1):qp1] <- matrix(1,(qp1-i),(qp1-i))*tt1[i+1];
  }
  Sigma.mat     <- sigmasq*vx;
  Itausq        <- matrix(rep(0,(qp1*qp1)), nrow=qp1, ncol=qp1);
  diag(Itausq)  <- rep(tausq,qp1);
  V             <- Sigma.mat + Itausq;
  ss=tt[2:qp1]-tt[1:q]; 
  D1mat=cbind(-diag(1/ss),matrix(0,q,1))+cbind(matrix(0,q,1),diag(1/ss));
  D2mat=cbind(-diag(1,q-1),matrix(0,q-1,1)) + cbind(matrix(0,q-1,1),diag(1,q-1));
  Phi.mat=D2mat%*%D1mat%*%V%*%t(D1mat)%*%t(D2mat); 
  wt=(yt[2:qp1]-yt[1:q])/ss; 
  ut=wt[2:q]-wt[1:q-1]; 
  ofn=((qp1)/2)*log(2*pi)+(0.5*log(det(Phi.mat))) + (0.5*(ut%*%ginv(Phi.mat)%*%ut)); 
  
  return(ofn)
}