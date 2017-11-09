#
# Bayesian Gaussian Process Regression
# For STAT8810
# Fall, 2017
# M.T. Pratola


# Bayesian Implementation of GP Regression Model
#
# Data:   y(xi) = eta(xi), i=1,...,n ; xi=1xK vector
#
#
# Likelihood: L(Y|.) propto. 1/(2pi*det(E))exp(-(Y')(Einv)(Y)) ; Y=[y(x1),...,y(xn)]
# where,
#         E=1/lambdaz*[R(.,.;rho)]
#
# Priors:
#         pi(lambdaz) propto. lambdaz^(az-1)exp(-bz*lambdaz)     [Gamma prior]
#         pi(rhoz_k) propto. rhoz_k^(rhoa-1)(1-rhoz_k)^(rhob-1)  [Beta prior]
#
# Posterior:
#         pi(.|Y) propto. L(Y|.)pi(lambdaz)pi(rhoz_1)*...*pi(rhoz_K)
#
# We fit the model using MCMC with adaptive stepwidths for the uniform proposal distributions.
#
#

# Note: requires dace.sim.r for the correlation function (rhogeodacecormat).

logp<-function(y,betaz,lambdaz,rhoz,l.dez,pi,perturb=1e-5)
{
  #lambdaz is the parameter for the mean regression covariance matrix.
  #rhoz is the parameter for each dimension of the tensored covariance function.
  az=pi$az; bz=pi$bz
  #pi is used to store all apriori hyperparameters.
  #az bz are hyperparameters for the Beta priors on covariance coefficients rho_d , d=1,2,...k=D
  
  Fz=pi$Fz; mu=pi$mu; lambdamu=pi$lambdamu;
  #Fz is set to be unital columns for simplicity. 
  #mu lambdamu are hyperparameters for the Noraml prior (conjugate) on beta parameter.
  #betaz=as.vector(betaz)
  #read those hyper parameters.
  if(any(rhoz<=0) || any(rhoz>=1))
    logposterior=-Inf
  else
  {
    n=length(y)
    #The observation is n*1 vector.
    
    k=length(rhoz)
    #The locations are n d*1 vectors, in this list rhoz, each dimension  rhoz_d (d=1,2,...,D) is stored as a component of the list.
    
    In=diag(n)
    #This is an identity matrix of dimension n.
    
    Rz=rhogeodacecormat(l.dez,rhoz)$R
    Ez=(1/lambdaz)*(Rz+In*perturb)
    #Generate the covariance matrix for the observation locations 
    cholEz=chol(Ez)
    Ezinv=chol2inv(cholEz)
    logdetEz=2*sum(log(diag(cholEz)))
    logpz=-1/2*logdetEz-1/2*t(y-Fz%*%betaz)%*%Ezinv%*%(y-Fz%*%betaz)
    #log density for the likelihood of observation.
    
    logplambdaz=(az-1)*log(lambdaz)-bz*lambdaz
    #log density for Gamma prior of lambdaz.
    
    logprhoz=0
    for(i in 1:k) logprhoz=logprhoz+(pi$rhoa[i]-1)*log(rhoz[i])+(pi$rhob[i]-1)*log(1-rhoz[i])
    #log density for Beta prior of rhoz.
    
    COVbetaz=(1/lambdamu)*diag(length(mu))
    COVbetazinv=chol2inv(chol(COVbetaz))
    cholCOVbetaz=chol(COVbetaz)
    logdetCOVbetaz=2*sum(log(diag(cholCOVbetaz)))
    logpbeta=-1/2*logdetCOVbetaz-1/2*t(betaz-mu)%*%COVbetazinv%*%(betaz-mu)
    #log density for prior of betaz
    
    logposterior=logpz+logplambdaz+logprhoz+logpbeta
  }
  logposterior
}


# Gibbs/MCMC Sampler
# y - observed response
# l.dez - distances in my format
# N - number of Gibbs/MCMC iterations for the sampler
# pi=list(az=5,bz=5,rhoa=rep(1,k),rhob=rep(5,k))
# mh=list(rr=.07)
regress<-function(y,l.dez,N,pi,mh,last=1000,adapt=TRUE)
{
  if(N<2000) stop("N>=2000 please\n")
  az=pi$az; bz=pi$bz
  Fz=pi$Fz; mu=pi$mu; lambdamu=pi$lambdamu;
  
  last=last-1
  n=length(y)
  k=length(l.dez)
  p=length(mu)
  
  #We shall check the dimension matches.
  cat("Observation dimension:",n,"\n Regressor dimension:",p,"\n Location dimension",k)
  
  draw.lambdaz=rep(NA,N)
  draw.betaz=matrix(NA,nrow=N,ncol=p)
  draw.rhoz=matrix(NA,nrow=N,ncol=k)
  
  rr=rep(mh$rr,k)
  #This is used for perturbation of the matrix before inversion.
  
  # initial guesses
  draw.lambdaz[1]=1/var(y)
  draw.betaz[1,]=rep(mean(y),p)
  draw.rhoz[1,]=rep(.5,k)
  
  #accept/reject ratio trackers for MH step for non-conjugate rho.
  accept.rhoz=rep(0,k)
  
  lastadapt=0
  
  In=diag(n)
  one=rep(1,n)
  rhoz.new=rep(NA,k)
  
  cat("\n Bayesian Gaussian Process Interpolation model")
  cat("\n The last ",last," samples from the posterior will be reported")
  cat("\n The stepwidth for uniform corr. param proposal distn is rr=",rr)
  cat("\n Prior params:  az=",pi$az," bz=",pi$bz)
  cat("\n ----------------------------------------------------------------\n\n\n")
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for(i in 2:N)
  {
    
    # Draw the correlation parameters (Metropolis-Hastings steps)
    draw.rhoz[i,]=draw.rhoz[i-1,]
    for(j in 1:k)
    {
      #This loop is just for proposing new move on each of D dimensions for the rhoz parameters used in covariance matrix.
      rhoz.new=draw.rhoz[i,]
      rhoz.new[j]=runif(1,draw.rhoz[i-1,j]-rr[j],draw.rhoz[i-1,j]+rr[j])
      #cat("\n")
      #print(draw.betaz)
      a=min(0,logp(y=y,betaz=draw.betaz[i-1,],lambdaz = draw.lambdaz[i-1],rhoz = rhoz.new,l.dez = l.dez,pi = pi)
             -logp(y=y,betaz=draw.betaz[i-1,],lambdaz = draw.lambdaz[i-1],rhoz = draw.rhoz[i,],l.dez = l.dez,pi = pi)
             )
      if(log(runif(1))<a)
      {
        draw.rhoz[i,j]=rhoz.new[j]
        accept.rhoz[j]=accept.rhoz[j]+1
      }
    }
    
    # Draw the marginal precision of the data (Gibbs step)
    Rz=rhogeodacecormat(l.dez,draw.rhoz[i,])$R
    Rz=Rz+In*1e-14  #cheat
    cholRz=chol(Rz)
    Rzinv=chol2inv(cholRz)
    draw.lambdaz[i]=rgamma(1,az+n/2,bz+0.5*t(y)%*%Rzinv%*%y)
    
    # Draw the marginal mean coefficient of the data (Gibbs step)
    Ip=diag(p)
    Rtmp=(draw.lambdaz[i])*Rzinv
    sigma_mu_tmp=(t(Fz)%*%Rtmp%*%Fz+lambdamu*Ip)
    sigma_mu=chol2inv(chol(sigma_mu_tmp))
    mean_mu=sigma_mu%*%(t( t(y)%*%Rtmp%*%Fz )+lambdamu*mu )
    draw.betaz[i,]=rmvnorm(n=1,mean=mean_mu,sigma=sigma_mu)
    
    #cat(i/N*100," percent complete\r")
    
    # Adaptive MCMC stuff:
    # adapt the proposal step every N/20 iterations, but only for the first 50% of the iterations
    if(adapt && i%%(N/20)==0 && i<(N*.5+1))
    {
      rate.rhoz=accept.rhoz/(i-lastadapt)
      cat("Adapting rates from ",rate.rhoz,"\n");
      for(j in 1:k)
        if(rate.rhoz[j]>.49 || rate.rhoz[j]<.39) rr[j]=rr[j]*rate.rhoz[j]/.44
      lastadapt=i
      accept.rhoz=rep(0,k)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  #cat("\n Complete.")
  rate.rhoz=accept.rhoz/(i-lastadapt)
  cat("[rate.rhoz=",rate.rhoz,"]\n")
  
  return(list(y=y,
              betaz=draw.betaz[(N-last):N],
              lambdaz=draw.lambdaz[(N-last):N],
              rhoz=as.matrix(draw.rhoz[(N-last):N,],pi=pi)))
}

plotfit<-function(fit,simtrue){
  par(mfrow=c(2,3))
  plot(fit$lambdaz,type='l',xlab="draw",ylab="lambdaz")
  abline(h=simtrue$lambdaytrue,col="red")
  plot(fit$rhoz,type='l',xlab="draw",ylab="rhoz")
  abline(h=simtrue$rhotrue,col="red")
  plot(fit$betaz,type='l',xlab="draw",ylab="betaz")
  abline(h=simtrue$betatrue,col="red")
  
  boxplot(fit$lambdaz,type='l',xlab="draw",ylab="lambdaz")
  abline(h=simtrue$lambdaytrue,col="red")
  boxplot(fit$rhoz,type='l',xlab="draw",ylab="rhoz")
  abline(h=simtrue$rhotrue,col="red")
  boxplot(fit$betaz,type='l',xlab="draw",ylab="betaz")
  abline(h=simtrue$betatrue,col="red")
  
}

# Draw realizations from posterior predictive distribution
# l.v - distances in my format where first n entries are for the observed data and 
#       the remaning are the prediction locations
# fit - output from running regress().
# eps.yy - fudge factor for inverting the observation correlation matrix
# eps - fudge factor for inverting the conditional covariance matrix
predict<-function(l.v,fit,eps.yy=0,eps=1e-6)
{
  n=length(fit$y)  # num observations
  m=nrow(l.v[[1]][[1]])-n  # num prediction points
  N=length(fit$lambdaz)
  
  draw.preds=matrix(0,nrow=N,ncol=m)
  
  for(i in 1:N)
  {
    Rall=rhogeodacecormat(l.v,fit$rhoz[i,])$R
    
    # Extract the sub-matrices we need
    Ryy=Rall[1:n,1:n]
    Rgg=Rall[(n+1):(n+m),(n+1):(n+m)]
    Rgy=Rall[(n+1):(n+m),1:n]
    Ryy.inv=chol2inv(chol(Ryy+diag(n)*eps.yy))
    
    # Mean of conditional distribution:
    m.cond=Rgy%*%Ryy.inv%*%(y)
    
    # Covariance of conditional distribution:
    E.cond=fit$lambdaz[i]^(-1)*(Rgg-Rgy%*%Ryy.inv%*%t(Rgy))
    
    # Let's generate a realization!
    L=t(chol(E.cond+diag(m)*eps))
    u=rnorm(m)
    draw.preds[i,]=m.cond + L%*%u
  }
  
  return(list(preds=draw.preds))
}

