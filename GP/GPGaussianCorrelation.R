#General Gaussian process regression with MLE estimation for Gaussian Correlation.
likelihood<-function(dat,mu,sigma,rr){
  GCF<-function(h,theta){
    if(theta>0){
      return(exp(-h^2*theta))
    }else{
      message("Parameter theta must >0!")
      return(-Inf)
    }
  }
  if(sigma<=0){return(-Inf)}
  n<-dim(dat)[1];
  x<-dat[,1:2];
  R<-matrix(NA,nrow=dim(dat)[1],ncol=dim(dat)[1]);
  for(k in 1:dim(dat)[1]){
    for(l in 1:dim(dat)[1]){
      h<-t(x[k,]-x[l,])%*%(x[k,]-x[l,]);
      R[k,l]<-GCF(sqrt(h),theta=rr);
    }
  }
  y<-dat[,3];
  res<-det(R)^(-.5);
  res<-res*(2*pi*sigma^2)^(-n/2);
  res<-res*exp( -(1/sigma^2)*t(y-mu)%*%solve(R)%*%(y-mu) ) 
  print(res)
  return(log(res));
}

dat1<-cbind(1:5,rnorm(5,mean=11));

MLE<-function(dattmp){
  cost<-function(x){
    return(likelihood(dattmp,mu=x[1],sigma=x[2],rr=x[3]))
  }
  res<-optim(c(0,10,100), fn=cost, method="BFGS",control = list(fnscale = -1) )$par
  return(res)
}
