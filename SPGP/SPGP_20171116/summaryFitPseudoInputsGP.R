#Different choices of number of pseudo-inputs

plotGpPseudoInputs<-function(fit, xPred = seq(0,1,length=100), cond = 0, ...) {
  
  #Load parameters
  xbarOpt <- fit$xbarOpt
  xbarStart <- fit$xbarStart
  sigma_noiseOpt <- fit$sigma_noiseOpt
  sigma_KOpt <- fit$sigma_KOpt
  rhoOpt <- fit$rhoOpt
  pseudoInputSelection <- fit$pseudoInputSelection
  y <- fit$y
  x <- fit$x
  
  M <- length(xbarOpt)
  N <- length(y)
  
  #Initialize plot
  plot(x, y, xlim = range(c(x, xbarOpt, xPred)),...)
  
  
  #Rescale input:
  xOrig <- x
  x <- (x-min(xOrig))/max((xOrig-min(xOrig)))
  xbarOpt <- (xbarOpt-min(xOrig))/max((xOrig-min(xOrig)))
  xPred <- (xPred-min(xOrig))/max((xOrig-min(xOrig)))
  
  if(pseudoInputSelection %in% c("none", "ml")) {
    
    pred <-    wrapperPredDistr(xNew = xPred, y = y, x = x, 
                                xbar = xbarOpt, 
                                sigma_noise = sigma_noiseOpt, 
                                sigma_K = sigma_KOpt, 
                                rho = rhoOpt, cond = cond)
  }
  
  if(pseudoInputSelection == "kl") {
    
    pred <-    wrapperApproxPosterior(xNew = xPred, y = y, x = x, 
                                xbar = xbarOpt, 
                                sigma_noise = sigma_noiseOpt, 
                                sigma_K = sigma_KOpt, 
                                rho = rhoOpt, cond = cond)
    
  }
  
  #Rescale outputs:
  xbarOpt <- xbarOpt*max((xOrig-min(xOrig))) + min(xOrig)
  xPred <-  xPred*max((xOrig-min(xOrig))) + min(xOrig)
  
  
  multiplier <- 1.96;
  
  lines(xPred, pred[1,], col = "red", lwd = 2)
  lines(xPred, as.numeric(pred[1,]) - sqrt( as.numeric( pred[2,]) ) * multiplier, col = "red", lty = 2)
  lines(xPred, as.numeric(pred[1,]) + sqrt( as.numeric( pred[2,]) ) * multiplier, col = "red", lty = 2)
  
  points(xbarOpt, rep(min(y), M), col="green3", pch=20)
  
  
  if(pseudoInputSelection %in% c("ml","kl")) {
    
    points(xbarStart, rep(max(y), M), col="blue", pch=20)
  
  }
  
}


mse <- function(fit, cond = 0, ...) {
  
  #Load parameters
  xbarOpt <- fit$xbarOpt
  xbarStart <- fit$xbarStart
  sigma_noiseOpt <- fit$sigma_noiseOpt
  sigma_KOpt <- fit$sigma_KOpt
  rhoOpt <- fit$rhoOpt
  pseudoInputSelection <- fit$pseudoInputSelection
  y <- fit$y
  x <- fit$x
  
  M <- length(xbarOpt)
  N <- length(y)

  if(pseudoInputSelection %in% c("none", "ml")) {
    
    predOnX <- wrapperPredDistr(xNew = x, y = y, x = x, 
                                xbar = xbarOpt, 
                                sigma_noise = sigma_noiseOpt, 
                                sigma_K = sigma_KOpt, 
                                rho = rhoOpt, cond = cond)
    
  }
  
  if(pseudoInputSelection == "kl") {
    
    predOnX <- wrapperApproxPosterior(xNew = x, y = y, x = x, 
                                      xbar = xbarOpt, 
                                      sigma_noise = sigma_noiseOpt, 
                                      sigma_K = sigma_KOpt, 
                                      rho = rhoOpt, cond = cond)
    
  }
  
  predError <- sum((y - as.numeric(predOnX[1,]))^2)
  
  return(predError)
  
}


# #This summary function allows you to analyze the fit of GP with pseudo inputs, closely.
# summaryPseudoInput <- function(fit,echo=T,return.data=F){
#   if(echo==T){
#     print("Gaussian Process with Sparse Pesudo-inputs ");
#     print(paste(fit$M," pseudo-inputs were chosen.",sep=""))
#     print(paste("The sigma_noise is ",fit$sigma_noiseOpt,sep=""))
#     print(paste("The marginal likelihood of pseudo-inputs is ",fit$likelihoodOpt,sep=""))
#     print("The optimal pseudo-inputs based on the observed data are")
#     print(as.data.frame(fit$xbarOpt))
#   }
#   if(return.data==T){
#     return(list(x=fit$x,y=fit$y)
#     )
#   }
# }