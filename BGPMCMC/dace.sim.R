#
# Gaussian Processes
# For STAT8810
# Fall, 2017
# M.T. Pratola


# Simulate an unconditional draw of a GP with the power exponential correlation model.
# alpha=2 gives us draws from the Gaussian correlation model.
#
# More info:
# design.distmat is the per-dimension distance matrices already in the usual
# list format, eg:
#
# l1=list(m1=abs(outer(design[,1],design[,1],"-")))
# l2=list(m2=abs(outer(design[,2],design[,2],"-")))
# l3=list(m3=abs(outer(design[,3],design[,3],"-")))
# l.design=list(l1=l1,l2=l2,l3=l3)
#
# rho is the vector of (0,1) correlation parameters
# s2 is the process variance
# se2 is the measurement error variance
# alpha is the correlation smoothness parameter
#
# conditioning may help the cholesky decomposition if we are generating
# a field at >=100 locations.  Can try conditioning=.Machine$double.eps
# or a reasonbly small number such as 1e-14.
#
sim.field<-function(design.distmat,rho,s2,se2=0,mu=0,alpha=2,conditioning=0)
{
  R=rhogeodacecormat(design.distmat,rho,alpha)$R
  nd=nrow(R)
  K=sqrt(s2)*chol(R+diag(nd)*conditioning)
  u=rnorm(nd)
  z=mu + t(K)%*%u + rnorm(nd,sd=se2)
  z
}


# Calculate the correlation matrix for the power exponential/Gaussian model.
#
# The correlation parameters are rho_1,...,rho_p for a p-dim space, each in [0,1] and
# we will have p distance (geoD) matrices.  We construct these in a list of lists, see
# for example the following:
# l1=list(m1=design.distmat.dim1)
# l2=list(m2=design.distmat.dim2)
# l=list(l1=l1,l2=l2)
#
rhogeodacecormat<-function(geoD,rho,alpha=2)
{
  rho=as.vector(rho)
  if(!is.vector(rho)) stop("non-vector rho!")
  if(any(rho<0)) stop("rho<0!")
  if(any(rho>1)) stop("rho>1!")
  if(any(alpha<1) || any(alpha>2)) stop("alpha out of bounds!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(rho)) stop("rho vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=nrow(geoD$l1$m1))
  for(i in 1:length(rho))
    R=R*rho[i]^(geoD[[i]][[1]]^alpha)
  
  return(list(R=R))
}

matern32<-function(geoD,theta)
{
  theta=as.vector(theta)
  if(!is.vector(theta)) stop("non-vector theta!")
  if(any(theta<0)) stop("theta<0!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(theta)) stop("theta vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=nrow(geoD$l1$m1))
  for(i in 1:length(theta))
  {
    D=(1+sqrt(3)*geoD[[i]][[1]]/theta[i])*exp(-sqrt(3)*geoD[[i]][[1]]/theta[i])
    R=R*D
  }
  
  return(list(R=R))
}

matern52<-function(geoD,theta)
{
  theta=as.vector(theta)
  if(!is.vector(theta)) stop("non-vector theta!")
  if(any(theta<0)) stop("theta<0!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(theta)) stop("theta vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=nrow(geoD$l1$m1))
  for(i in 1:length(theta))
  {
    D=(1+sqrt(3)*geoD[[i]][[1]]/theta[i]+5*geoD[[i]][[1]]^2/(5*theta[i]^2))*exp(-sqrt(5)*geoD[[i]][[1]]/theta[i])
    R=R*D
  }
  
  return(list(R=R))
}

wendland1<-function(geoD,theta)
{
  theta=as.vector(theta)
  if(!is.vector(theta)) stop("non-vector theta!")
  if(any(theta<0)) stop("theta<0!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(theta)) stop("theta vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=nrow(geoD$l1$m1))
  for(i in 1:length(theta))
  {
    D=(1-geoD[[i]][[1]]/theta[i])
    D[D<0]=0
    D=D^4*(1+4*geoD[[i]][[1]]/theta[i])
    R=R*D
  }
  
  return(list(R=R))	
}

wendland2<-function(geoD,theta)
{
  theta=as.vector(theta)
  if(!is.vector(theta)) stop("non-vector theta!")
  if(any(theta<0)) stop("theta<0!")
  if(!is.list(geoD)) stop("wrong format for distance matrices list!")
  if(length(geoD)!=length(theta)) stop("theta vector doesn't match distance list")
  
  R=matrix(1,nrow=nrow(geoD$l1$m1),ncol=nrow(geoD$l1$m1))
  for(i in 1:length(theta))
  {
    D=(1-geoD[[i]][[1]]/theta[i])
    D[D<0]=0
    D=D^6*(1+6*geoD[[i]][[1]]/theta[i]+35*geoD[[i]][[1]]^2/(3*theta[i]^2))
    R=R*D
  }
  
  return(list(R=R))	
}

generalized.wendland<-function(geoD,theta,kap)
{
  d=length(geoD)
  mu=(d+1)/2  # strictly speaking we need mu>=(d+1)/2 but we can change kappa also so
  # this is fair to assume.
  
  if(length(theta)>1) stop("theta is incorrect dimension\n")
  if(length(kap)>1) stop("kappa is incorrect dimensions\n")
  if(mu<((d+1)/2) ) stop("mu does not satisfy constraints\n")
  if(kap<0) stop("kappa > 0 required\n")
  
  D=matrix(0,nrow=nrow(geoD$l1$m1),ncol=nrow(geoD$l1$m1))
  for(i in 1:length(geoD))
    D=D+(geoD[[i]][[1]]^2)
  D=sqrt(D)
  
  if(kap==0) {
    # kap=0 is essential the Askey correlation
    D=D/theta
    R=1-D
    R[R<0]=0
    R=R^(mu+kap)
  }
  else
  {
    # library(fields) implements the general case
    R=fields::Wendland(D,theta=theta,dimension=d,k=kap)
  }
  
  rm(D)
  return(list(R=R))
}

logl<-function(rho,y,design.distmat,alpha=2,conditioning=0)
{
  n=length(y)
  R=rhogeodacecormat(design.distmat,rho,alpha)$R
  cR=chol(R+diag(n)*conditioning)
  Rinv=chol2inv(cR)
  
  s2.hat=(1/n)*t(y)%*%Rinv%*%y
  logdetR.div2=sum(log(diag(cR)))
  l=-n/2*log(s2.hat)-logdetR.div2-n/2
  
  return(l)
}


