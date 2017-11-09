#Sparse Peudo-input Gaussian Process(SPGPs)
#Source: Snelson, Edward, and Zoubin Ghahramani. "Sparse Gaussian processes using pseudo-inputs." Advances in neural information processing systems. 2006.
#Author: Hengrui Luo(luo.619@osu.edu)
library(mvtnorm)
K<-function(point1,point2){
  #This is the covariance function we specify for GP regression model.
  norm<-as.numeric( t(point1-point2)%*%(point1-point2) );
  #norm for stationary kernel
  return( exp(-norm/0.01) )
}
makeCovariance<-function(vector1,vector2){
  #This is used for generating the covariance matrix using the covariance function provided.
  cov<-matrix(NA,nrow=length(vector1),ncol=length(vector2));
  for(i in 1:dim(cov)[1]){
    for(j in 1:dim(cov)[2]){
      cov[i,j]<-K(vector1[i],vector2[j]);
    }
  }
  return(cov);
}
p<-function(y,x,barx,barf,sigma=1){
  #Calculating the single data point likelihood.
  if(length(barx)!=length(barf)){stop("p: barx and barf must have conformal dimension!")}
  K_xx<-makeCovariance(x,x);
  K_M<-makeCovariance(barx,barx);
  k_x<-makeCovariance(barx,x);
  
  mu<-t(k_x)%*%chol2inv(chol(K_M))%*%barf;
  sigmamat<-K_xx-t(k_x)%*%chol2inv(chol(K_M))%*%k_x+sigma^2;
  result<-dmvnorm(x=y, mean=mu, sigma=sigmamat, log=FALSE)
  return(result)
  #Test sample: p(y=c(1,2),x=c(1,2.5),barx=c(1.25,1.75,2.25),barf=c(3,4,5))
}
prior_pseudoTarget<-function(barf,barx){
  #Calculating the prior on pseudo targets barf;
  if(length(barx)!=length(barf)){stop("p: barx and barf must have conformal dimension!")}
  K_M<-makeCovariance(barx,barx);
  
  mu<-rep(0,length(barf));
  sigmamat<-K_M;
  result<-dmvnorm(x=barf, mean=mu, sigma=sigmamat, log=FALSE)
  return(result)
}
posterior_pseudoTarget<-function(x,y,barf,barx,sigma=1){
  #Calculating the posterior on pseudo targets barf;
  #x,y and real data,
  #barx,barf are pseudo input and pseudo data respectively.
  
  if(length(barx)!=length(barf)){stop("p: barx and barf must have conformal dimension!")}
  k_n<-makeCovariance(barx,x);
  K_M<-makeCovariance(barx,barx);
  K_MN<-makeCovariance(barx,x);
  K_NM<-makeCovariance(x,barx);
  K_nn<-makeCovariance(x,x);
  lambda<-diag(
    diag(K_nn-t(k_n)%*%chol2inv(chol(K_M))%*%k_n)
  )
  
  Q_M<-K_M+K_MN%*%chol2inv(chol( lambda+sigma^2*diag(dim(lambda)[1]) ))%*%K_NM;

  
  mu<-K_M%*%chol2inv(chol(Q_M))%*%K_MN%*%chol2inv( chol( lambda+sigma^2*diag(dim(lambda)[1]) ) )%*%y
  sigmamat<-K_M%*%chol2inv(chol(Q_M))%*%K_M;
  result<-dmvnorm(x=barf, mean=mu, sigma=sigmamat, log=FALSE)
  return(result)
  #Test sample: posterior_pseudoTarget(x=c(1,2,3),y=c(1,4,9),barf=c(3.5,7.5),barx = c(1.5,2.5))
}
predictive_pseudoTarget<-function(ynew,xnew,x,y,barx,barf,sigma=1){
  k_n<-makeCovariance(barx,x);
  K_M<-makeCovariance(barx,barx);
  K_MN<-makeCovariance(barx,x);
  K_NM<-makeCovariance(x,barx);
  K_nn<-makeCovariance(x,x);
  K_starstar<-makeCovariance(xnew,xnew);
  kstar<-makeCovariance(xnew,x);
  lambda<-diag(
    diag(K_nn-t(k_n)%*%chol2inv(chol(K_M))%*%k_n)
  )
  
  Q_M<-K_M+K_MN%*%chol2inv( chol( lambda+sigma^2*diag(dim(lambda)[1]) ))%*%K_NM;
  
  mustar<-t(kstar)%*%chol2inv(chol(Q_M))%*%K_MN%*%chol2inv( chol( lambda+sigma^2*diag(dim(lambda)[1]) ))%*%y
  sigmastar<-K_starstar-t(kstar)%*%( chol2inv(chol(K_M))-chol2inv(chol(Q_M)) )%*%kstar+sigma^2
  result<-dmvnorm(x=ynew, mean=mustar, sigma=sigmastar, log=FALSE)
  return(result)
}
marginal_likelihood<-function(y,x,barx,sigma=1){
  k_n<-makeCovariance(barx,x);
  K_M<-makeCovariance(barx,barx);
  K_MN<-makeCovariance(barx,x);
  K_NM<-makeCovariance(x,barx);
  K_nn<-makeCovariance(x,x);
  lambda<-diag(
    diag(K_nn-t(k_n)%*%chol2inv(chol(K_M))%*%k_n)
  )
  d<-dim(lambda)[1];
  mu<-rep(0,d);
  sigmamat<-K_NM%*%chol2inv(chol(K_M))%*%K_MN+as.numeric(sigma^2)*diag(d);
  result<-dmvnorm(x=y, mean=mu, sigma=sigmamat, log=TRUE)
  return(result)
}
#You shall specify the number of pseudo input N before you optimize the marginal likelihood to yield parameter estimations.
M=50;
N=5;
set.seed(11)
#set.seed(1) will yield the local extrema and a very bad optimal point as claimed by [Titsias].
#Source: Titsias, Michalis K. "Variational learning of inducing variables in sparse Gaussian processes." International Conference on Artificial Intelligence and Statistics. 2009.
fulldatay<-rnorm(M);
fulldatax<-rnorm(M);
wrapper1<-function(x){
  if(x[N+1]<=0){return(-Inf)}
  return( marginal_likelihood(y=fulldatay,x=fulldatax,barx=x[1:N],sigma=x[N+1]) 
          )
}

random_init<-c(rnorm(N),0.01)
wrapper1_res<-optim(par=random_init, fn=wrapper1, method="BFGS",control = list(fnscale = -1) )$par

random_init
wrapper1_res

plot(fulldatax,fulldatay)
points(x=random_init[1:N],y=rep(0,N),col="blue",pch=20)
points(x=wrapper1_res[1:N],y=rep(0,N),col="red",pch=20)
#we can see that the optimized pseudo locations are somehow more scattered around.