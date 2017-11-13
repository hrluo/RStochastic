library(rjmcmc)
## Comparing two binomial models -- see Barker & Link (2013) for further details.
y=c(8,16); 
sumy=sum(y)
n=c(20,30); 
sumn=sum(n)
#Binomial success rates. Suppose yi~Binom(ni, pi ) for i=1,2 and we have observations 
#y1 = 8, n1 = 20
#y2 = 16,n2 = 30
#To evaluate the evidence for p1 != p2 versus p1 = p2.

L1=function(p){if((all(p>=0))&&(all(p<=1))) sum(dbinom(y,n,p,log=TRUE)) else -Inf}
#This is the log likelihood function of model 1, that is to say, a multi-nomial model with parameter \Vec(p).
#dim model 1=dim of vector p.

L2=function(p){if((p[1]>=0)&&(p[1]<=1)) sum(dbinom(y,n,p[1],log=TRUE)) else -Inf}
#This is the log likelihood function of model 2, that is to say, a binomial model with parameter p.
#dim model 2=1.

p.prior1=function(p){sum(dbeta(p,1,1,log=TRUE))}
#This is the prior function for parameter of model 1, i.e. pi( theta^(1) )

p.prior2=function(p){dbeta(p[1],1,1,log=TRUE)+dbeta(p[2],17,15,log=TRUE)}
#This is the prior function for parameter of model 2, i.e. pi( theta^(2) )


## The bijection we used in this example is psi=(p1,p2) for a three dimensional multinomial model (p1+p2+p3=1)
## g(psi)=g(p1,p2)=(  (n1*p1+n2*p2)/(n1+n2), u ).
## u is the supplementary proposal random variable.

g1=function(psi){p=psi}
# This is the first component of the bijection.

g2=function(psi){w=n[1]/sum(n); p=c(w*psi[1]+(1-w)*psi[2],psi[2])}
# This is the second component of the bijection.

ginv1=function(p){p}
ginv2=function(p){c(sum(n)/n[1]*p[1]-n[2]/n[1]*p[2],p[2])}
# These two are the inverse function of the bijection.
# We must calculate Jacobian corresponding to the bijection we construct using these two functions.

draw1=function(){rbeta(2,y+1,n-y+1)}
draw2=function(){c(rbeta(1,sumy+1,sumn-sumy+1),rbeta(1,17,15))}
#These two functions are posterior draw from model 1 and 2 respectively. We know their closed form because of the Beta conjugate priors.
#If we do not know their closed form, a Metropolis-Hasting algorithm shall be implemented.

coda1=matrix(rbeta(2000,y+1,n-y+1), nrow=1000, ncol=2, byrow=TRUE) 
#This is the closed form of full conditional posterior due to conjugacy of (separable) Beta-multinomial.

coda2=matrix(c(rbeta(1000,sumy+1,sumn-sumy+1),rbeta(1000,17,15)), nrow=1000, ncol=2)
#This is the closed form of full conditional posterior due to conjugacy of (separable) Beta-binomial.

out_DB=defaultpost(coda=list(coda1,coda2), 
                #If we already have the posterior draws, then we feed them to this function.
                likelihood=list(L1,L2),
                #A list of likelihood of model 1,2.
                param.prior=list(p.prior1,p.prior2), 
                #A list of N(=2) functions specifying the prior distributions for each model-specific parameter vector. We are dealing with two models here.
                model.prior=c(1,1), 
                #A numeric vector of the prior model probabilities. Note that this argument is not required to sum to one as it is automatically normalised.
                chainlength=1000)
# Generate synthetic 'MCMC output' for a model with 3 parameters. 
# There is one column per parameter, and 1000 iterations.
# In the returned value we have used 'Default Bijections'  
# This function uses a 'Default Bijection' scheme based on approximating each posterior by a multivariate normal distribution.

out_RJ=rjmcmcpost(post.draw=list(draw1,draw2),
                  #If we DO NOT have the posterior draws, then we feed them to this function.
                  g=list(g1,g2), 
                  ginv=list(ginv1,ginv2),
                  likelihood=list(L1,L2), 
                  param.prior=list(p.prior1,p.prior2),
                  model.prior=c(1,1), 
                  chainlength=1000)
# In the returned value we have used 'g Bijections' supplied in the parameter g=list(...), ginv=list(...),
# This function uses a 'g Bijection' scheme based on user-specified bijection that is used in reversible jump MCMC.

# rjmethods refer to a collection of overloaded functions that analyzes the rj class
densities(out_RJ)
densities(out_DB)

probplot(out_RJ)
probplot(out_DB)

#Posteror predictive distribution from fitted models containing two possible models 1 and 2.
predPost<-function(x,rj){
  probs<-c(rj$result$`Posterior Model Probabilities`)
  u<-runif(1);
  sample<-NULL;

    if(u<probs[1]){
      sample<-exp(L1(x))
    }else{
      sample<-exp(L2(x))
    }

  return(sample)
}

#Plot the fitted sample density.
predPostSample_RJ<-c()
predPostSample_DB<-c()
xseq=seq(-5,5,.01)
for(x in xseq){
  predPostSample_RJ<-c( predPostSample_RJ,predPost(x=x,rj=out_RJ) )
  predPostSample_DB<-c( predPostSample_DB,predPost(x=x,rj=out_DB) )
}
plot(xseq,predPostSample_RJ,type="l",col="red")
lines(xseq,predPostSample_DB+.01,type="l",col="blue")
#These two methods yield almost identical posterior.


posteriorRJ = getsampler(modelfit=out_RJ$psidraws, sampler.name="posteriorRJ") 
posteriorDB = getsampler(modelfit=out_DB$psidraws, sampler.name="posteriorDB") 
# This successfully defines a function named posterior1b but also defines an identical function corresponding to the value
# of sampler.name, i.e. the default "post.draw" in this case.
# To construct such a function, we need to obtain a posterior sample and supply it to modelfit parameter in getsampler().
# Usually we want 
set.seed(100)

posteriorRJ()
posteriorDB()
# The output is a vector consisting of values of 2/4 parameters. That is to say we draw from the posterior of the parameters psi.
