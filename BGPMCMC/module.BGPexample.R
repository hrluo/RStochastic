#
# Bayesian Gaussian Process Example
# For STAT8810
# Fall, 2017
# M.T. Pratola
library(fOptions)
library(mvtnorm
        )
source("dace.sim.R")
source("regress.original.R")
##########################################################################################################################
# First, generate the response.
##########################################################################################################################
cat("\nGenerating response\n")
n=12
k=2
p=1
simtrue=list(rhotrue=0.2,lambdaytrue=1,betatrue=rep(3.14,p))
pi=list(az=5,bz=5,rhoa=rep(1,k),rhob=rep(5,k),Fz=matrix(1,nrow=n,ncol=p),mu=rep(1,p),lambdamu=1 )
design=as.matrix(runif.sobol(n,1))

l1=list(m1=outer(design[,1],design[,1],"-"))
l.dez=list(l1=l1)
R=(1/simtrue$lambdaytrue)*rhogeodacecormat(l.dez,c(simtrue$rhotrue))$R
L=chol(R)
u=rnorm(nrow(R))
z=t(L)%*%u+pi$Fz%*%simtrue$betatrue

# now set up our observed data:
y=z
l1=list(m1=outer(design[,1],design[,1],"-"))
l.dez=list(l1=l1)


##########################################################################################################################
# Now fit the GP model
##########################################################################################################################
#source("regression.r")



# Run MCMC using a proposal width of 1e-5 for the rho parameters.
mh=list(rr=1e-5)
fit=regress(y,l.dez,5000,pi,mh,last=2000,adapt=FALSE)
plotfit(fit,simtrue)

mh=list(rr=0.05)
fit2=regress(y,l.dez,5000,pi,mh,last=2000,adapt=TRUE)
plotfit(fit2,simtrue)

# Predict the last fit on a grid of 100 equally spaced points
grid=as.matrix(seq(0,1,length=11))
design.all=rbind(design,grid)
l1=list(m1=outer(design.all[,1],design.all[,1],"-"))
l.v=list(l1=l1)
fitp=predict(l.v,fit,eps.yy=1e-4)


par(mfrow=c(1,2))
# Plot the posterior predictions using mean and sd to quantify uncertainty (pointwise)
plot(design,y,pch=20,col="red",cex=2,xlim=c(0,1),
     ylim=c(0,5),xlab="x",
     main="Predicted mean response +/- 2s.d.")
for(i in 1:nrow(fitp$preds))
  lines(grid,fitp$preds[i,],col="grey",lwd=0.25)
mean=apply(fitp$preds,2,mean)
sd=apply(fitp$preds,2,sd)
lines(grid,mean-1.96*sd,lwd=0.75,col="black")
lines(grid,mean+1.96*sd,lwd=0.75,col="black")
lines(grid,mean,lwd=2,col="blue")

# Plot the posterior predictions using median and quantiles to quantify uncertainty (pointwise)
plot(design,y,pch=20,col="red",cex=2,xlim=c(0,1),ylim=c(0,5),
     xlab="x",main="Predicted median, q.025 and q.975")
for(i in 1:nrow(fitp$preds))
  lines(grid,fitp$preds[i,],col="grey",lwd=0.25)
med=apply(fitp$preds,2,quantile,0.5)
q.025=apply(fitp$preds,2,quantile,0.025)
q.975=apply(fitp$preds,2,quantile,0.975)
lines(grid,q.025,lwd=0.75,col="black")
lines(grid,q.975,lwd=0.75,col="black")
lines(grid,mean,lwd=2,col="blue")

