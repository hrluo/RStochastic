#Sparse Peudo-input Gaussian Process(SPGPs)
#Source: Snelson, Edward, and Zoubin Ghahramani. "Sparse Gaussian processes using pseudo-inputs." Advances in neural information processing systems. 2006.
#Author: Hengrui Luo(luo.619@osu.edu)
library(mvtnorm)
K<-function(point1,point2,c=1,ell=2){
  #This is the covariance function we specify for GP regression model.
  norm<-as.numeric( t(point1-point2)%*%(point1-point2) );
  #norm for stationary kernel
  return( c*exp(-0.5*norm/ell^2) )
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
#We adopt following notations from
#Source: Quiñonero-Candela, Joaquin, and Carl Edward Rasmussen. "A unifying view of sparse approximate Gaussian process regression." Journal of Machine Learning Research 6.Dec (2005): 1939-1959.
Kmatrix<-function(a,b){
  return(makeCovariance(a,b) )
}
Qmatrix<-function(a,b,u){
  result<-Kmatrix(a,u)%*%chol2inv(chol(Kmatrix(u,u)))%*%Kmatrix(u,b);
  #Will return a matrix in conformal with a,b.
}
Diagmatrix<-function(mat){
  #return a matrix retaining only diagonal entries.
  return(diag(diag(mat)));
}
marginal_likelihood<-function(y,u,f,sigma=1){
  library(mvtnorm)
  #y is the full data  of size M.
  #u is the pseudo-inputs of size N, or latent inducing variables that we think we can reduce to, which are what we want to come up with by maximizing the marginal.
  #f is the full inputs of size M.
  
  d<-length(f);
 
  Lambda<-Diagmatrix( Kmatrix(a=f,b=f)-Qmatrix(a=f,b=f,u=u) );
  #Kmatrix(a=f,b=f) is a M*M matrix.
  #Qmatrix(a=f,b=f,u=u) is a M*M matrix.
  #Lambda is a M*M matrix.
  #print(Lambda)
  
  K1<-Kmatrix(a=f,b=u);
  #print(K1)
  #K1 is a M*N matrix.
  K2<-Kmatrix(a=u,b=f);
  #print(K2)
  #K2 is a N*M matrix.
  Ku<-Kmatrix(a=u,b=u);
  #K_M is a N*N matrix.
  #print(Ku)
  Kuinv<-chol2inv(chol(Ku))
  

  mu<-rep(0,d);
  sigmaMat<-K1%*%Ku%*%K2+Lambda+sigma^2*diag(d);
  y_projected<-K1%*%Kuinv%*%u;
  y_projected<-as.vector(y_projected);
  #print(mu)
  #print(sigmaMat)
  #print(y_projected)
  
  result<-dmvnorm(x=y_projected, mean=mu, sigma=sigmaMat, log=TRUE)
  return(result)
}
#You shall specify the number of pseudo input N before you optimize the marginal likelihood to yield parameter estimations.
set.seed(11)
#set.seed(1) will yield the local extrema and a very bad optimal point as claimed by [Titsias].
#Source: Titsias, Michalis K. "Variational learning of inducing variables in sparse Gaussian processes." International Conference on Artificial Intelligence and Statistics. 2009.
M=50;#full data set contains M observations.
N=5;#We only want 6 pseudo inputs.
fulldatax<-seq(1,M,1);
fulldatay<-2*fulldatax+rnorm(M)
wrapper1<-function(x){
  sigma_param<-tail(x, n=1)
  if(sigma_param<=0){return(-Inf)}
  #the last entry is used for representing sigma in the Gaussian model.
  return( marginal_likelihood(y=fulldatay,f=fulldatax,u=x[1:N],sigma=sigma_param) 
          )
}
random_init<-c(rnorm(N),0.01)
wrapper1_res<-optim(par=random_init, fn=wrapper1, method="Nelder-Mead",control = list(fnscale = -1) )$par

random_init
wrapper1_res

plot(fulldatax,fulldatay)
points(x=random_init[1:N],y=rep(0,N),col="blue",pch=20,xlim=range(fulldatax),ylim=range(fulldatay))
points(x=wrapper1_res[1:N],y=rep(0,N),col="red",pch=20,xlim=range(fulldatax),ylim=range(fulldatay))
#we can see that the optimized pseudo locations are somehow more scattered around..
u1<-wrapper1_res[1:N];
sig1<-tail(wrapper1_res,1);
#We yield these parameter estimations from the optimization result of the marginal likelihood as shown above.
predictive_distribution<-function(ynew,y,u=u1,f,sigma=sig1){
  if(length(ynew)!=1){stop("New response must be a length 1 value!")}
  Lambda<-Diagmatrix( Kmatrix(a=f,b=f)-Qmatrix(a=f,b=f,u=u) );
  Qstarf<-Qmatrix(a=ynew,b=f,u=u);
  Q1<-Qmatrix(a=f,b=f,u=u)+Lambda;
  Q2<-chol2inv(chol(Q1));
  
  munew<-Qstarf%*%Q2%*%y;
  
  Kstarstar<-Kmatrix(a=ynew,b=ynew);
  Qfstar<-Qmatrix(a=f,b=ynew,u=u);
  Q3<-Qstarf%*%Q2%*%Qfstar
  sigmaMatnew<-Kstarstar-Q3;
  
  result<-dnorm(ynew, mean=munew, sd=sqrt(sigmaMatnew));
  return(result);
}
#Let us plot the predictive distribution to see if this works well enough.
probs<-c();
probs_real<-c();
iter<-1;
for(k in fulldatay){
  probs[iter]<-predictive_distribution(ynew=k,y=fulldatay,f=fulldatax);
  probs_real[iter]<-dnorm(fulldatay,mean=2*fulldatax[iter]);
  iter<-iter+1;
}

plot(fulldatay,log(probs),type="l",xlab="Y value",xlim=range(fulldatay))
lines(fulldatay,log(probs_real),col="red",xlab="Y value",xlim=range(fulldatay))
#The red line is the true probability such a Y should have under the model Y=2*X+N(0,1).
#The black line is the predictive distribution probability.
#They turn out to be approximately the same,which mean the sparse GP is working pretty good!

#Now we encase the procedure into one function
SPGP<-function(full.x=fulldatax,full.y=fulldatay,
               N.pi=3,cov.fun=K){
  wrapper<-function(x){
    sigma_param<-tail(x, n=1)
    if(sigma_param<=0){return(-Inf)}
    #the last entry is used for representing sigma in the Gaussian model.
    return( marginal_likelihood(y=full.y,f=full.x,u=x[1:N.pi],sigma=sigma_param) )
  }
  random_init<-c(rnorm(N.pi),0.01);
  wrapper_res<-optim(par=random_init, fn=wrapper, method="Nelder-Mead",control = list(fnscale = -1) )$par;
  u_opt<-wrapper_res[1:N.pi];
  sigma_opt<-tail(wrapper_res,1);
  result_fun<-function(x){
    return(  predictive_distribution(ynew=x,y=full.y,u=u_opt,f=full.x,sigma=sigma_opt)  )
  }
  return(result_fun)
  #return a predictive distribution.
}
message("SPGP initialized!")

fun1<-SPGP()
fun2<-SPGP(N.pi=length(fulldatay))

#fun2 is the exact posterior because it uses all the data without any "pseudo-ness"