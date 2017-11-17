#Covariance kernel (default Gaussian)
#-------------------------------------
#Inputs: scalars a, b, sigma_K (SD Kernel), 
#                rho (correlation for obs. 1 unit apart), 
#                p (power of exponent, p=2 Gaussian, p=1 exponential)
#Outputs: scalar of covariance value.

K <- function(a, b, sigma_K, rho, p = 2){
  output <- sigma_K^2 * rho^((a-b)^p);
  return(output) 
}

#Covariance matrix between two inputs. 
#-------------------------------------
#Inputs: vectors x1, x2
#        scalar sigma_K, rho
#Outputs: covariance matrix between x1 and x2.

Kmatrix <- function(x1, x2, sigma_K, rho) {
  output <- outer(x1, x2, FUN = K, 
                  sigma_K = sigma_K, rho = rho)
  return(output)
}

#Marginal likelihood function of pseudo-inputs xbar given observations (y,x).
#----------------------------------------------------------------------------
#Input:  
# - N*1 vector y, response 
# - N*1 vector x, inputs 
# - M*1 vector xbar, pseudo-inputs 
# - scalar sigma_noise, sd noise on responses
# - scalar sigma_K, SD covariance kernel GP 
# - scalar rho, correlation for covariance kernel GP
# - scalar cond, perturbation coefficient used to address the numerical issue in inverting matrices
#Output: the value of the likelihood function.
#Source: Snelson, Edward, and Zoubin Ghahramani. "Sparse Gaussian processes using pseudo-inputs." In Advances in neural information processing systems, pp. 1257-1264. 2006.

marginalLikelihood <- function(y, x, xbar, sigma_noise, sigma_K, rho, cond = 0){
  #browser()
  N <- length(y)
  M <- length(xbar)
  
  Knm <- Kmatrix(x,xbar,sigma_K, rho)
  Kmm <- Kmatrix(xbar,xbar,sigma_K, rho)
  KmmInv <- chol2inv(chol(Kmm + cond*diag(M)))
  
  Knm_KmmInv_Knm <- Knm %*% KmmInv %*% t(Knm)
  
  Lambda <- diag( diag( Kmatrix(x,x,sigma_K, rho) - Knm_KmmInv_Knm))
  
  muMargLik <- rep(0,N)
  covMatrixMargLik <- Knm_KmmInv_Knm + Lambda + (sigma_noise^2)*diag(N)
  
  # y_projected <- as.vector(Knm %*% KmmInv %*% xbar)
  # 
  # result<- (-.5 * (N * log(2 * pi) + as.numeric(determinant(covMatrixMargLik, logarithm = T)$modulus) ) + 
  #             (-.5 * (y_projected - muMargLik) %*% chol2inv(chol(covMatrixMargLik + cond*diag(N))) %*% (y_projected - muMargLik)))
  result<- (-.5 * (N * log(2 * pi) + as.numeric(determinant(covMatrixMargLik, logarithm = T)$modulus) ) + 
              (-.5 * (y - muMargLik) %*% chol2inv(chol(covMatrixMargLik + cond*diag(N))) %*% (y - muMargLik)))
  
  return(result)
}

#Wrapper marginal likelihood for optimization over pseudo-inputs (xbar), 
# sd error (sigma_noise) and covariance kernel parameters (sigma_K and rho).
#----------------------------------------------------------------------------
wrapperMargLik <- function(input, y, x, cond = 0){
  
  M<-length(input)-3;
  
  xbar <- input[1:M]
  sigma_noise <- rev(input)[3]
  sigma_K <- rev(input)[2]
  rho <- rev(input)[1]

  if(sigma_noise<=0){return(-1e20)}
  if(sigma_K<=0){return(-1e20)}
  if(rho<=0 | rho>=1){return(-1e20)}
  
  trial <- try(marginalLikelihood(y = y, x = x,
                     xbar=xbar, 
                     sigma_noise = sigma_noise, 
                     sigma_K = sigma_K, 
                     rho = rho,
                     cond = cond), silent = T)
  
  if(class(trial) == "try-error") {
    margLikValue <- NA 
  } else {
    margLikValue <- trial
  }
  
  return(margLikValue)
}

wrapperMargLikMax <- function(input, y, x, cond = 0){
  
  minusMargLik <- (-1) * wrapperMargLik(input, y, x, cond)
  
  return(minusMargLik)
}


#Wrapper marginal likelihood for optimization ON sd error (sigma_noise) and 
# covariance kernel parameters (sigma_K and rho) ONLY (the pseudo-inputs xbar are FIXED).
#----------------------------------------------------------------------------
wrapperMargLikGivenPseudoInputs <- function(input, y, x, xbar, cond = 0){
  
  M <- length(xbar)
  
  sigma_noise <- input[1]
  sigma_K <- input[2]
  rho <- input[3]
  
  if(sigma_noise<=0){return(-1e20)}
  if(sigma_K<=0){return(-1e20)}
  if(rho<=0 | rho>=1){return(-1e20)}
  
  trial <- try( marginalLikelihood(y = y, x = x,
                                     xbar=xbar, 
                                     sigma_noise = sigma_noise, 
                                     sigma_K = sigma_K, 
                                     rho = rho,
                                     cond = cond), silent = T)
  
  if(class(trial) == "try-error") {
    margLikValue <- NA 
  } else {
    margLikValue <- trial
  }
  
  return(margLikValue)
}

wrapperMargLikGivenPseudoInputsMax <- function(input, y, x, xbar, cond = 0){

  minusMargLik <- (-1) * wrapperMargLikGivenPseudoInputs(input, y, x, xbar, cond)

  return(minusMargLik)
}


#Predictive distribution
#-------------------------
#Input:  
# - scalar xNew, new input value where perform predictions
# - N*1 vector y, response 
# - N*1 vector x, inputs 
# - M*1 vector xbar, pseudo-inputs 
# - scalar sigma_noise, sd noise on responses
# - scalar sigma_K, SD covariance kernel GP 
# - scalar rho, correlation for covariance kernel GP
# - scalar cond, perturbation coefficient used to address the numerical issue in inverting matrices
#Output: mu (predicted mean), sigma2 (variance o)
#Source: Snelson, Edward, and Zoubin Ghahramani. "Sparse Gaussian processes using pseudo-inputs." In Advances in neural information processing systems, pp. 1257-1264. 2006.

predictiveDistribution <- function(xNew, y, x, xbar, sigma_noise, sigma_K, rho, cond = 0) {

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

#Wrapper predictive distribution that allows the vectorial evaluation
#---------------------------------------------------------------------
wrapperPredDistr <- Vectorize(FUN = predictiveDistribution,
                              vectorize.args = "xNew")


# Lower bound corresponding to KL divergence between approximate posterior and real posterior (Titsias, 2009)
#------------------------------------------------------------------------------------------------------------
#Input:  
# - N*1 vector y, response 
# - N*1 vector x, inputs 
# - M*1 vector xbar, pseudo-inputs 
# - scalar sigma_noise, sd noise on responses
# - scalar sigma_K, SD covariance kernel GP 
# - scalar rho, correlation for covariance kernel GP
# - scalar cond, perturbation coefficient used to address the numerical issue in inverting matrices
#Output: the value of the likelihood function.
#Source: Titsias, Michalis K. "Variational learning of inducing variables in sparse Gaussian processes." International Conference on Artificial Intelligence and Statistics. 2009.

F_V <- function(y, x, xbar, sigma_noise, sigma_K, rho, cond = 0) {

  N <- length(y)
  M <- length(xbar)
  
  Knn <- Kmatrix(x,x,sigma_K, rho)
  Knm <- Kmatrix(x,xbar,sigma_K, rho)
  Kmm <- Kmatrix(xbar,xbar,sigma_K, rho)
  KmmInv <- chol2inv(chol(Kmm + cond*diag(M)))
  
  Knm_KmmInv_Knm <- Knm %*% KmmInv %*% t(Knm)
  
  K_tilde <- Knn-Knm_KmmInv_Knm
  
  covMatrixMargLik <- sigma_noise^2 * diag(N) + Knm_KmmInv_Knm
  result <- ( (-.5 * (N * log(2 * pi) + as.numeric(determinant(covMatrixMargLik, logarithm = T)$modulus)) + 
                (-.5 * t(y) %*% chol2inv(chol(covMatrixMargLik + cond*diag(N))) %*% (y))) 
               - 1/(2*sigma_noise^2) * sum(diag(K_tilde)) )
  
  
  return(result)
}

#Wrapper of F_V for optimization over pseudo-inputs (xbar), 
# sd error (sigma_noise) and covariance kernel parameters (sigma_K and rho).
#----------------------------------------------------------------------------
wrapperF_V <- function(input, y, x, cond = 0){
  
  M <- length(input)-3;
  
  xbar <- input[1:M]
  sigma_noise <- rev(input)[3]
  sigma_K <- rev(input)[2]
  rho <- rev(input)[1]
  
  if(sigma_noise<=0){return(-1e20)}
  if(sigma_K<=0){return(-1e20)}
  if(rho<=0 | rho>=1){return(-1e20)}
  
  trial <- try(F_V(y = y, x = x,
                        xbar=xbar, 
                        sigma_noise = sigma_noise, 
                        sigma_K = sigma_K, 
                        rho = rho,
                        cond = cond), silent = T)
  
  if(class(trial) == "try-error") {
    margLikValue <- NA 
  } else {
    margLikValue <- trial
  }
  
  return(margLikValue)
}

wrapperF_V_Max <- function(input, y, x, cond = 0){
  
  minusMargLik <- (-1) * wrapperF_V(input, y, x, cond)
  
  return(minusMargLik)
}

#Approximate posterior distribution to make predictions in new points
#--------------------------------------------------------------------
#Input:  
# - scalar xNew, new input value where perform predictions
# - N*1 vector y, response 
# - N*1 vector x, inputs 
# - M*1 vector xbar, pseudo-inputs 
# - scalar sigma_noise, sd noise on responses
# - scalar sigma_K, SD covariance kernel GP 
# - scalar rho, correlation for covariance kernel GP
# - scalar cond, perturbation coefficient used to address the numerical issue in inverting matrices
#Output: mu (predicted mean), sigma2 (variance o)
#Source: Titsias, Michalis K. "Variational learning of inducing variables in sparse Gaussian processes." International Conference on Artificial Intelligence and Statistics. 2009.

approxPosterior <- function(xNew, y, x, xbar, sigma_noise, sigma_K, rho, cond = 0) {
  
  N <- length(y)
  M <- length(xbar)
  
  kNewM <- Kmatrix(xNew, xbar, sigma_K, rho)
  kNewNew <- Kmatrix(xNew, xNew, sigma_K, rho)
  Knm <- Kmatrix(x,xbar,sigma_K, rho)
  Kmm <- Kmatrix(xbar,xbar,sigma_K, rho)
  KmmInv <- chol2inv(chol(Kmm + cond*diag(M)))
  
  Knm_KmmInv_Knm <- Knm %*% KmmInv %*% t(Knm)

  K_tilde <- Knn-Knm_KmmInv_Knm
  
  Sigma <- chol2inv(chol(Kmm +  1/sigma_noise^2 * t(Knm) %*% Knm + cond * diag(M)))
    
  muPost <- 1/sigma_noise^2 * Kmm %*% Sigma %*%  t(Knm) %*% y
  
  A <- Kmm %*% Sigma %*% Kmm
  B <- KmmInv %*% A %*% KmmInv
    
  muPredDistr <- kNewM %*% KmmInv %*% muPost
  sigma2PredDistr <- kNewNew - kNewM %*% KmmInv %*% t(kNewM) + kNewM %*% B %*% t(kNewM) 
  
  return(list(mu = muPredDistr,
              sigma2 = sigma2PredDistr))
}

#Wrapper approximate posterior distribution that allows the vectorial evaluation
#---------------------------------------------------------------------
wrapperApproxPosterior <- Vectorize(FUN = predictiveDistribution,
                                    vectorize.args = "xNew")





#Fit GP with pseudo-inputs 
#--------------------------
#Input:  
# - scalar xNew, new input value where perform predictions
# - N*1 vector y, response 
# - N*1 vector x, inputs 
# - M*1 vector xbar, pseudo-inputs. If provided, "pseudoInputSelection" must be "none"  
# - character pseudoInputSelection, method used to select pseudo-inputs. Must be one of: "none", "ml", "kl".
# - scalar sigma_noiseStart, starting value in optimization for sd noise on responses
# - scalar sigma_KStart, starting value in optimization for SD covariance kernel GP 
# - scalar rhoStart, starting value in optimization for correlation for covariance kernel GP
# - scalar cond, perturbation coefficient used to address the numerical issue in inverting matrices
#Output: mu (predicted mean), sigma2 (variance o)

fitGpPseudoInputs <- function(y, x, M = NULL,
                              pseudoInputSelection,
                              xbar = NULL,
                              xbarStart = seq(0,1,length = M), 
                              sigma_noiseStart = .05,
                              sigma_KStart = 1, 
                              rhoStart = .1,
                              sigma_noiseUpper = 10,
                              sigma_KUpper = 10,
                              cond = 0,
                              maxTime = 24*60*60) {
  
  if(!is.null(xbar) & pseudoInputSelection != "none") {
    stop("If you specify the pseudo-inputs, 'pseudoInputSelection' must be 'none'")
  }
  
  if(is.null(xbar) & pseudoInputSelection == "none") {
    stop("If 'pseudoInputSelection' is 'none', you must specify the pseudo-inputs")
  }
  
  if(is.null(M) & !is.null(xbar)) {
    M <- length(xbar)
  }
  
  N <- length(y)
  
  #Rescale input to be in [0,1]:
  xOrig <- x
  x <- (xOrig-min(xOrig))/max((xOrig-min(xOrig)))
  
  if(pseudoInputSelection == "none") {
    
    parStart <- c(sigma_noiseStart, sigma_KStart, rhoStart)
    
    lower <- c(0,0,0)
    upper <- c(sigma_noiseUpper,sigma_KUpper,1)
    out.GenSA <- GenSA(par = parStart, 
                       lower = lower, 
                       upper = upper, 
                       fn = wrapperMargLikGivenPseudoInputsMax,
                       y = y, x = x, xbar= xbar, cond = cond, 
                       control=list(verbose=T, max.time = maxTime)) 
    
    xbarOpt <- xbar
    sigma_noiseOpt <- out.GenSA$par[1]
    sigma_KOpt <- out.GenSA$par[2]
    rhoOpt <- out.GenSA$par[3]  
    
    optimalValue <- out.GenSA$value
    convergence <- NA
    
    # objOpt <- optim(par = parStart, 
    #                 fn = wrapperMargLikGivenPseudoInputs, 
    #                 y = y, x = x, xbar= xbar, cond = cond, 
    #                 method="L-BFGS-B", 
    #                 lower = c(0,0,0), upper = c(1,10,1),
    #                 control = list(fnscale = -1,trace = 1, maxit = 5000))
    # parOpt<-objOpt$par
    # 
    
    # xbarOpt <- xbar
    # sigma_noiseOpt <- parOpt[1]
    # sigma_KOpt <- parOpt[2]
    # rhoOpt <- parOpt[3]  
    # 
    # optimalValue <- objOpt$value
    # convergence <- objOpt$convergence
  }
  
  if(pseudoInputSelection == "ml") {
    
    parStart <- c(xbarStart, sigma_noiseStart, sigma_KStart, rhoStart)

    #1:
    # objOpt <- optim(par = parStart,
    #                 fn = wrapperMargLik,
    #                 y = y, x = x, cond = cond,
    #                 method="L-BFGS-B",
    #                 lower = c(rep(0,M),0,0,0), upper = c(rep(1,M),1,10,1),
    #                 control = list(fnscale = -1, trace = 1, maxit = 5000))
    # parOpt<-objOpt$par
    # 
    # xbarOpt <- parOpt[1:M]
    # sigma_noiseOpt <- parOpt[M+1]
    # sigma_KOpt <- parOpt[M+2]
    # rhoOpt <- parOpt[M+3]
    # 
    # optimalValue <- objOpt$value
    # convergence <- objOpt$convergence

    # #2:Genetic algorithms
    # GA <- ga(type = "real-valued", fitness =  wrapperMargLik,
    #          min = c(rep(0,M), 0,0,0), max = c(rep(1,M), 10,10,1),
    #          y = y, x = x, cond = cond,
    #          popSize = 200, maxiter = 500, parallel=T)
    # 
    # objOpt <- optim(par = as.numeric(GA@solution),
    #                 fn = wrapperMargLik,
    #                 y = y, x = x, cond = cond,
    #                 method="Nelder-Mead",
    #                 control = list(fnscale = -1, trace = 1))
    # parOpt<-as.numeric(objOpt)
    
    # #3: Conjugate gradient
    # objOpt <- optim(par = parStart,
    #                 fn = wrapperMargLik,
    #                 y = y, x = x, cond = cond,
    #                 method="CG",
    #                 control = list(fnscale = -1, trace = 1, parscale = c(rep(1,M),1e-1,1,1), type = 3,maxit = 7))
    # parOpt<-objOpt$par
    # 
    # xbarOpt <- parOpt[1:M]
    # sigma_noiseOpt <- parOpt[M+1]
    # sigma_KOpt <- parOpt[M+2]
    # rhoOpt <- parOpt[M+3]
    # 
    # optimalValue <- objOpt$value
    # convergence <- objOpt$convergence

    lower <- c(rep(0,M),0,0,0)
    upper <- c(rep(1,M),sigma_noiseUpper,sigma_KUpper,1)
    out.GenSA <- GenSA(par = parStart,
                       lower = lower,
                       upper = upper,
                       fn = wrapperMargLikMax,
                       y = y, x = x, cond = cond,
                       control=list(verbose=T, max.time = maxTime))

    xbarOpt <- out.GenSA$par[1:M]
    sigma_noiseOpt <- out.GenSA$par[M+1]
    sigma_KOpt <- out.GenSA$par[M+2]
    rhoOpt <- out.GenSA$par[M+3]

    optimalValue <- out.GenSA$value
    convergence <- NA
    

  }
  
  if(pseudoInputSelection == "kl") {
    
    parStart <- c(xbarStart, sigma_noiseStart, sigma_KStart, rhoStart)
    
    # objOpt <- optim(par = parStart,
    #                 fn = wrapperF_V,
    #                 y = y, x = x, cond = cond,
    #                 method="Nelder-Mead",
    #                 lower = c(rep(0,M),0,0,0), upper = c(rep(1,M),1,10,1),
    #                 control = list(fnscale = -1,trace = 1, maxit = 5000))
    # parOpt<-objOpt$par
    # 
    # xbarOpt <- parOpt[1:M]
    # sigma_noiseOpt <- parOpt[M+1]
    # sigma_KOpt <- parOpt[M+2]
    # rhoOpt <- parOpt[M+3]
    # 
    # optimalValue <- objOpt$value
    # convergence <- objOpt$convergence

    lower <- c(rep(0,M),0,0,0)
    upper <- c(rep(1,M),sigma_noiseUpper,sigma_KUpper,1)
    out.GenSA <- GenSA(par = parStart,
                       lower = lower,
                       upper = upper,
                       fn = wrapperF_V_Max,
                       y = y, x = x, cond = cond,
                       control=list(verbose=T, max.time = maxTime))

    xbarOpt <- out.GenSA$par[1:M]
    sigma_noiseOpt <- out.GenSA$par[M+1]
    sigma_KOpt <- out.GenSA$par[M+2]
    rhoOpt <- out.GenSA$par[M+3]

    optimalValue <- out.GenSA$value
    convergence <- NA
     
  }
  
  #Rescale output:
  xbarOpt <- xbarOpt*max((xOrig-min(xOrig))) + min(xOrig)
  xbarStart <- xbarStart*max((xOrig-min(xOrig))) + min(xOrig)
  x <- x*max((xOrig-min(xOrig))) + min(xOrig)
  
  return(list(xbarOpt = xbarOpt,
              sigma_noiseOpt = sigma_noiseOpt,
              sigma_KOpt = sigma_KOpt,
              rhoOpt = rhoOpt,
              pseudoInputSelection = pseudoInputSelection,
              xbarStart = xbarStart,
              y = y, x = x,
              optimalValue = optimalValue,
              convergence = convergence))
}


# 
# #The following function allowing an adaptive choice of number of pseudo-inputs M.
# fitGpPseudoInputsAdapt<-function(y,x,M.range=seq(1,10)){
#   if(min(M.range)<1){stop("The number of pseudo inputs must be greater than 1!")}
#   result<-matrix(NA,ncol=2,nrow=length(M.range));
#   colnames(result)<-c("M","maximalLL");
#   iter<-1;
#   for(m in M.range){
#     value_m<-fitGpPseudoInputs(y,x,M=m)$likelihoodOpt;
#     value_m<-as.numeric(value_m)
#     result[iter,]<-c(m,value_m)
#     iter<-iter+1;
#   }
#   par(mfrow=c(1,1))
#   result<-as.data.frame(result)
#   plot(x=result$M,y=result$maximalLL,
#        xlab="M",ylab="maximal marginal likelihood",main="Adaptive Choice of number of pseudo inputs M.")
#   return(results)
# }

preferPseudoInputs <- function(x, y, M){
  
  N <- length(y)
  ySort <- y[order(x)]
  
  yDiff <- abs(ySort[2:N] - ySort[1:(N-1)])
  
  pTemp <- c(yDiff[1],  pmax(yDiff[1:(N-2)],yDiff[2:(N-1)]) , yDiff[N-1])
  p <- pTemp/max(pTemp)
  
  xbar <- x[sample.int(N, M, prob = p)]
  return(xbar)
}