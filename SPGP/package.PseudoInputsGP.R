#Sparse Peudo-input Gaussian Process(SPGPs)
#Source: Snelson, Edward, and Zoubin Ghahramani. "Sparse Gaussian processes using pseudo-inputs." Advances in neural information processing systems. 2006.
#Author: Hengrui Luo (luo.619@osu.edu) and Giovanni Nattino (nattino.1@osu.edu)


#Function that simulates a realization of a GP
simulateGP_GaussianCorr <- function(n, p, theta, sigma2) {
  
  X <- matrix(NA, nrow = n, ncol = p)
  #This is the data matrix with original data inputs of size n and p regressors.
  #####
  #This could be replaced with real data. Ideally,
  for(i in 1:p) {
    
    x_i <- seq(0, 1, length = n)
    X[,i] <- x_i
    #Generate a uniform grid, when put into use, this can be substituted with real data inputs.
    
    diff_x_i <- abs(outer(x_i, x_i, "-"))
    #distance between any two of these two locations.
    
    R_i <- exp(- theta[i] * diff_x_i^2 )
    #Correlation function for direction i, i=1,...,p
    
    L_i <- t(chol(R_i + diag(n)*.Machine$double.eps*100))
    
    if(i == 1) {
      L <- L_i
    } else {
      L <- L %x% L_i
    }
    
  }
  L <- sqrt(sigma2) * L
  
  Z <- rnorm(n^p, mean = 0, sd = 1)
  #Normal random sample
  
  Yvect <- L %*% Z
  #Realization (in vector form)
  
  Y <- aperm(array(Yvect, dim = rep(n,p))) 
  #Grid with realizations of GP
  
  return(list(Y = Y, 
              Yvect = Yvect,
              X = X))  
}

#Covariance kernel (default Gaussian). 
#Inputs: scalars a, b
#Outputs: scalar of covariance value.
K <- function(a, b, sigma_K, rho, type="Gaussian"){
  # You can choose from different covariance kernels to fit the Gaussian process, however this does not affect the Gaussian-Gaussian hierarchy.
  # Moreover, you need package "fields" and pay attention to the parametrization.
  if(type=="Gaussian"){
    output<-sigma_K^2 * rho ^ (a-b)^2;
  }
  if(type=="Matern"){
    require(fields)
    output<-sigma_K^2 * Matern(d=a-b,smoothness = rho);
  }
  if(type=="Exponential"){
    require(fields)
    output<-sigma_K^2 * Exponential(d=a-b);
  }
  return(output) 
}

#Covariance matrix between two inputs. 
#Inputs: vectors x1, x2
#Outputs: matrix of covariance between x1, x2.
Kmatrix <- function(x1, x2, sigma_K, rho) {
  output <- outer(x1, x2, FUN = K, 
                  sigma_K = sigma_K, rho = rho)
  return(output)
}

# Marginal likelihood function of pseudo-inputs xbar given observations (y,x).
# Input:  vectors 
#         y (N*1, response), 
#         x (N*1, inputs), 
#         xbar (M*1, pseudo-inputs), 
#         sigma_noise (scalar,sd noise on responses),
#         sigma_K (scalar, sd GP), 
#         rho (scalar, correlation for covariance kernel of GP)
#         cond (the perturbation coefficient used for addressing the numerical issue in inverting the covariance matrix)
# Output: the value of the likelihood function.
marginal_likelihood <- function(y, x, xbar, sigma_noise, sigma_K, rho, cond = 0){
  
  N <- length(y)
  M <- length(xbar)
  
  Knm <- Kmatrix(x,xbar,sigma_K, rho)
  Kmm <- Kmatrix(xbar,xbar,sigma_K, rho)
  KmmInv <- chol2inv(chol(Kmm + cond*diag(M)))
  
  Knm_KmmInv_Knm <- Knm %*% KmmInv %*% t(Knm)
  
  Lambda <- diag( diag( Kmatrix(x,x,sigma_K, rho) - Knm_KmmInv_Knm))
  
  muMargLik <- rep(0,N)
  covMatrixMargLik <- Knm_KmmInv_Knm + Lambda + sigma_noise^2*diag(N)

  y_projected <- as.vector(Knm %*% KmmInv %*% xbar)
  
  result<- (-.5 * log(det(2 * pi * covMatrixMargLik)) +
             (-.5 * (y_projected - muMargLik) %*% chol2inv(chol(covMatrixMargLik + cond*diag(N))) %*% (y_projected - muMargLik)))
  return(result)
}

wrapperMargLik <- function(input, cond = 0){
  #This function is used for wrapping the marginal likelihood function, for optimization over pseudo inputs.
  #Source: Titsias, Michalis K. "Variational learning of inducing variables in sparse Gaussian processes." International Conference on Artificial Intelligence and Statistics. 2009.
  sigma_noise <- rev(input)[3]
  sigma_K <- rev(input)[2]
  rho <- rev(input)[1]
  
  
  M<-length(input)-3;
  if(sigma_noise<=0){return(-Inf)}
  if(sigma_K<=0){return(-Inf)}
  if(rho<=0 | rho>=1){return(-Inf)}
  
  margLikValue <- marginal_likelihood(y = y, x = x,
                                     xbar=input[1:M], 
                                     sigma_noise = sigma_noise, 
                                     sigma_K = sigma_K, rho = rho,
                                     cond = cond)
  return(margLikValue)
}

posterior_distribution <- function(y, x, xbar, sigma_noise, sigma_K, rho, cond = 0,mu=rep(0,M), A=diag(M) ) {
  #The prior distribution for the pseudo-inputs is defined to be N(mu,A) by deafult to specify the prior selection of pseudo-inputs
  #This result is given in Titsias(8-10)
  N <- length(y)
  M <- length(xbar)
  
  Knn <- Kmatrix(x,x,sigma_K, rho)
  Knm <- Kmatrix(x,xbar,sigma_K, rho)
  Kmm <- Kmatrix(xbar,xbar,sigma_K, rho)
  KmmInv <- chol2inv(chol(Kmm + cond*diag(M)))
  
  Knm_KmmInv_Knm <- Knm %*% KmmInv %*% t(Knm)
  B<-KmmInv%*%A%*%KmmInv
  
  post_mean<-Knm%*% KmmInv %*%mu
  post_cov<-Knn-Knm_KmmInv_Knm+Knm%*%B%*%Kmn
  
  require(mvtnorm)
  K_tilde<-Knn-Knm_KmmInv_Knm;
  marginal_likelihood_bound<-dmvnorm(x=y,
                                     mu=rep(0,N),
                                     sigma=sigma_noise*diag(N)+Knm_KmmInv_Knm  )-1/(2*sigma_noise)*tr(K_tilde)
  FUN_prior<-function(xbar){
    #This is the optimal prior on the peuso inputs that allows us to obtain optimal selections of the pseudo-input points.
    Sigma<-Kmm+sigma_noise^(-2)%*%Kmn%*%Knm;
    Sigma<-chol2inv( chol(Sigma) );
    A1<-Kmm%*%Sigma%*%Kmm;
    mu1<-sigma_noise^(-2)%*%Kmm%*%Sigma%*%Kmn%*%y;
    output <- dmvnorm(x=xbar,mu=mu1,sigma = A1)
    return(output)
  }
  return(list(mu = post_mean,
              sigma2 = post_cov,
              bound=marginal_likelihood_bound,
              optimal_prior=FUN_prior))
}


predictive_distribution <- function(xNew, y, x, xbar, sigma_noise, sigma_K, rho, cond = 0) {
  #This is the predictive posterior given by Snelson-Zoubin
  N <- length(y)
  M <- length(xbar)
    
  kNewM <- Kmatrix(xNew, xbar, sigma_K, rho)
  kNewNew <- Kmatrix(xNew, xNew, sigma_K, rho)
  Knm <- Kmatrix(x,xbar,sigma_K, rho)
  Kmm <- Kmatrix(xbar,xbar,sigma_K, rho)
  KmmInv <- chol2inv(chol(Kmm + cond*diag(M)))
  
  Knm_KmmInv_Knm <- Knm %*% KmmInv %*% t(Knm)
  
  Lambda <- diag( diag( Kmatrix(x,x,sigma_K, rho) - Knm_KmmInv_Knm))
  
  Lambda_sigma2I_inv <- chol2inv(chol(Lambda + sigma_noise^2*diag(N) + cond*diag(N)))
  
  Qmm_inv <- chol2inv(chol(Kmm + t(Knm) %*% Lambda_sigma2I_inv %*% Knm + cond*diag(M)))
  
  muPredDistr <- kNewM %*% Qmm_inv %*% t(Knm) %*% Lambda_sigma2I_inv %*% y
  sigma2PredDistr <- kNewNew - kNewM %*% (KmmInv - Qmm_inv) %*% t(kNewM) + sigma_noise^2
  
  return(list(mu = muPredDistr,
              sigma2 = sigma2PredDistr))
}

wrapperPredDistr <- Vectorize(FUN = predictive_distribution,
                              vectorize.args = "xNew")


fitGpPseudoInputs <- function(y, x, M, 
                              xbarStart = runif(M), sigma_noiseStart = 1,
                              sigma_KStart = 1, rhoStart = .1,
                              cond = 0) {
  
  N <- length(y)
  
  parStart <- c(xbarStart, sigma_noiseStart, sigma_KStart, rhoStart)
  objOpt <- optim(par = parStart, 
                  fn = wrapperMargLik, cond = 1e-7, 
                  method="Nelder-Mead", control = list(fnscale = -1))
  parOpt<-objOpt$par;
  valOpt<-objOpt$value;
  xbarOpt <- parOpt[1:M]
  sigma_noiseOpt <- parOpt[M+1]
  sigma_KOpt <- parOpt[M+2]
  rhoOpt <- parOpt[M+3]
  
  return(list(xbarOpt = xbarOpt,
              sigma_noiseOpt = sigma_noiseOpt,
              sigma_KOpt = sigma_KOpt,
              rhoOpt = rhoOpt,
              likelihoodOpt=valOpt,
              x=x,
              y=y,
              M=M))
}
#The following function allowing an adaptive choice of number of pseudo-inputs M.
fitGpPseudoInputsAdapt<-function(y,x,M.range=seq(1,10)){
  if(min(M.range)<1){stop("The number of pseudo inputs must be greater than 1!")}
  result<-matrix(NA,ncol=2,nrow=length(M.range));
  colnames(result)<-c("M","maximalLL");
  iter<-1;
  for(m in M.range){
    value_m<-fitGpPseudoInputs(y,x,M=m)$likelihoodOpt;
    value_m<-as.numeric(value_m)
    result[iter,]<-c(m,value_m)
    iter<-iter+1;
  }
  par(mfrow=c(1,1))
  result<-as.data.frame(result)
  plot(x=result$M,y=result$maximalLL,
       xlab="M",ylab="maximal marginal likelihood",main="Adaptive Choice of number of pseudo inputs M.")
  return(results)
}

####################
# Parameters of GP #
####################

set.seed(123)

# Number of observed points
N <- 10000 

#SD kernel
sigma_K <- 1

#Correlation
rho <- .1

#Draw of a GP
simulation <- simulateGP_GaussianCorr(p = 1, n = N, theta = -log(rho), sigma2 = sigma_K)
x <- as.vector(simulation$X)
y <- as.vector(simulation$Y)


#Single fit fixed number of pseudo inputs
M <- 10
xbarStart <- runif(M)
fit <- fitGpPseudoInputs(y, x, M, 
                         xbarStart = xbarStart, sigma_noiseStart = 1,
                         sigma_KStart = 1, rhoStart = .1, cond = 0)


#Grid to display points
par(mfrow = c(1,1))
plot(x, y, xlim = range(c(x,fit$xbarOpt, xbarStart)))
points(xbarStart, rep(max(y), M), col="blue", pch=20)
points(fit$xbarOpt, rep(min(y), M), col="green3", pch=20)

xPred <- seq(0,1, length = 100)
pred <- wrapperPredDistr(xPred, y, x, 
                         fit$xbarOpt, fit$sigma_noiseOpt, fit$sigma_KOpt, fit$rhoOpt, cond = 1e-6)
lines(xPred, pred[1,], col = "red")

#Different choices of number of pseudo-inputs

plotGpPseudoInputs<-function(M=4,xPred=seq(0,1,length=100)) {
  xbarStart <- runif(M)
  fit <- fitGpPseudoInputs(y, x, M, 
                           xbarStart = xbarStart, sigma_noiseStart = 1,
                           sigma_KStart = 1, rhoStart = .1, cond = 0)
  
  #Grid to display points
  plot(x, y, xlim = range(c(x,fit$xbarOpt, xbarStart,xPred)), 
       main = paste0("M =", M))
  points(xbarStart, rep(max(y), M), col="blue", pch=20)
  #The blue points are randomly chosen startup pseudo input points.
  points(fit$xbarOpt, rep(min(y), M), col="green3", pch=20)
  #The green points are chosen pseudo inputs points.
  #xPred <- seq(0,1, length = 100)
  pred <-    wrapperPredDistr(xPred, y, x, 
                           fit$xbarOpt, fit$sigma_noiseOpt, fit$sigma_KOpt, fit$rhoOpt, cond = 1e-6)
  
  fullpred<- wrapperPredDistr(x, y, x, 
                              fit$xbarOpt, fit$sigma_noiseOpt, fit$sigma_KOpt, fit$rhoOpt, cond = 1e-6)
  fullpredError<-y-as.numeric(fullpred[1,])
  fullpredError<-sum(fullpredError^2)
  #This is the prediction error of this specific fit for the whole full data set.
  
  #Here we can replace the "xbar" option with the optimal design chosen pseudo input locations and therefore we can actually see how is the performance of Snelson-Zoubin compared with our solution.
  
  multiplier<-1.96;
  
  lines(xPred, pred[1,], col = "red", lwd = 2)
  lines(xPred, as.numeric(pred[1,])-sqrt( as.numeric(pred[2,]) )*multiplier, col = "magenta", lty = 2)
  lines(xPred, as.numeric(pred[1,])+sqrt( as.numeric(pred[2,]) )*multiplier, col = "magenta", lty = 2)
  
  #The read lines are fit under the "sparsified" Gaussian process with only pseudo inputs and thus fasten up the fitting process.
  return(list(fit=fit,fullpredError=fullpredError) )
}
par(mfrow = c(2,2))
fit2<-plotGpPseudoInputs(M=2)
fit3<-plotGpPseudoInputs(M=3)
fit4<-plotGpPseudoInputs(M=4)
fit5<-plotGpPseudoInputs(M=5)
#This summary function allows you to analyze the fit of GP with pseudo inputs, closely.
summaryPseudoInput<-function(fit,echo=T,return.data=F){
  if(echo==T){
    print("Gaussian Process with Sparse Pesudo-inputs ");
    print(paste(fit$M," pseudo-inputs were chosen.",sep=""))
    print(paste("The sigma_noise is ",fit$sigma_noiseOpt,sep=""))
    print(paste("The marginal likelihood of pseudo-inputs is ",fit$likelihoodOpt,sep=""))
    print("The optimal pseudo-inputs based on the observed data are")
    print(as.data.frame(fit$xbarOpt))
  }
  if(return.data==T){
    return(list(x=fit$x,y=fit$y)
           )
  }
}
summaryPseudoInput(fit2)
summaryPseudoInput(fit3)
summaryPseudoInput(fit4)
summaryPseudoInput(fit5)

#Minimal Example
plotGpPseudoInputs(M=4,xPred=seq(0,2,length=100))