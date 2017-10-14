# Example of Predictive Distribution
source("dace.sim.r")
set.seed(99)
n=4
x1=runif(n)
l1=list(m1=abs(outer(x1,x1,"-")))
l.dez=list(l1=l1)
s2=1
# sim.field defaults to Gaussian correlation
R=matern32(l.dez,theta=3)$R
L=t(chol(R+diag(n)*.Machine$double.eps*100));Z=rnorm(n);
z=L%*%Z;
####################################This line is the generating covariance.
# Now we assume that the z's are observed error-free:
y=z
# Let's write down the predictive distribution 
# over a fine grid
ng=100;xg=seq(0,1,length=ng);X=c(x1,xg);l1=list(m1=abs(outer(X,X,"-")));l.dez=list(l1=l1)
Rall=generalized.wendland(l.dez,theta=0.05,kap=1)$R
####################################This is the model covariance function.

# Extract the sub-matrices we need
Ryy=Rall[1:n,1:n]
Rgg=Rall[(n+1):(n+ng),(n+1):(n+ng)]
Rgy=Rall[(n+1):(n+ng),1:n]
Ryy.inv=chol2inv(chol(Ryy))

# Mean of conditional distribution:
m.cond=Rgy%*%Ryy.inv%*%y

# Covariance of conditional distribution:
E.cond=s2*(Rgg-Rgy%*%Ryy.inv%*%t(Rgy))


# Now the predictive distn is N(m.cond,E.cond).
# Let's generate a realization!
L=t(chol(E.cond+diag(ng)*1e-5))
u=rnorm(ng)
z.cond=m.cond + L%*%u

# And make a plot
plot(x1,y,pch=20,col="red",cex=2,xlim=c(0,1),ylim=c(-2,2))
m.cond1<-m.cond
E.cond1<-E.cond
z.cond1<-z.cond



# Generate some more realizations from the generating
plot(x1,y,pch=20,col="red",cex=2,xlim=c(0,1),ylim=c(-2,2))
for(i in 1:100){
  ng=100;xg=seq(0,1,length=ng);X=c(x1,xg);l1=list(m1=abs(outer(X,X,"-")));l.dez=list(l1=l1)
  Rall=matern32(l.dez,theta=3)$R
  ####################################This line is the generating covariance.
  Ryy=Rall[1:n,1:n]
  Rgg=Rall[(n+1):(n+ng),(n+1):(n+ng)]
  Rgy=Rall[(n+1):(n+ng),1:n]
  Ryy.inv=chol2inv(chol(Ryy))
  m.cond=Rgy%*%Ryy.inv%*%y
  E.cond=s2*(Rgg-Rgy%*%Ryy.inv%*%t(Rgy))
  L=t(chol(E.cond+diag(ng)*1e-5))
  u=rnorm(ng)
  z.cond=m.cond + L%*%u
  lines(xg,z.cond,col="grey")
  
}

lines(xg,m.cond1,lwd=5,col="black")
lines(xg,m.cond1-1.96*sqrt(diag(E.cond1)),lwd=2,col="black")
lines(xg,m.cond1+1.96*sqrt(diag(E.cond1)),lwd=2,col="black")
#lines(xg,z.cond1,col="grey")