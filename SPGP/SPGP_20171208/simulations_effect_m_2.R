#Sparse Peudo-input Gaussian Process 
#Author: Hengrui Luo (luo.619@osu.edu) and Giovanni Nattino (nattino.1@osu.edu)

library(GenSA, lib.loc="C:/Users/imrp/Documents/R/win-library/3.4")
library(DiceKriging, lib.loc="C:/Users/imrp/Documents/R/win-library/3.4")
library(doParallel, lib.loc="C:/Users/imrp/Documents/R/win-library/3.4")
setwd("C:/Users/imrp/Desktop/GP")

library(GenSA)
library(DiceKriging)
library(doParallel)

###############################
# Source Files with Functions #
###############################

source("simulateGP.R")
source("fitPseudoInputsGP.R")
source("summaryFitPseudoInputsGP.R")


#########################################################
# Simulation 1: One simulation, investigate effect of m #
#########################################################

#Fit a realization of the GP with a GP with pseudo-inputs.
set.seed(12345)

# Number of observed points
N <- 100

#SD kernel
sigma_K <- 1

#Correlation
rho <- .1e-8

#Draw of a GP
simulation <- simulateGP_GaussianCorr(p = 1, 
                                      n = N, 
                                      theta = -log(rho), 
                                      sigma2 = sigma_K)
x <- as.vector(simulation$X)
yGp <- as.vector(simulation$Y)

y <- yGp + rnorm(length(yGp), sd = .05)

plot(x,y)

############################################################################
# Example (NO SIMULATION): Fit all the methods with fixed m and save plots #
############################################################################

#Number of Pseudo-inputs
M <- 10

start <- Sys.time()

#Provide random pseudo-inputs: grid 
#----------------------------------
fit <- fitGpPseudoInputs(y, x, M,
                         xbar = seq(0,1, length = M),
                         pseudoInputSelection = "none",
                         cond = 1e-14)
plotGpPseudoInputs(fit, cond = 1e-7)
mse(fit,cond = 1e-7)

save(fit, file = "fit_givenPI_grid.Rdata")

end <- Sys.time()
end-start


start <- Sys.time()

#Provide random pseudo-inputs: SRS from data points
#---------------------------------------------------
fit <- fitGpPseudoInputs(y, x, M,
                         xbar = sample(x = x, size = M),
                         pseudoInputSelection = "none",
                         cond = 1e-14)
plotGpPseudoInputs(fit,cond = 1e-7)
mse(fit, cond = 1e-7)

save(fit, file = "fit_givenPI_SRS.Rdata")

end <- Sys.time()
end-start


start <- Sys.time()

#Provide random pseudo-inputs: Preferential sampling from points
#---------------------------------------------------------------
xbar <- preferPseudoInputs(x = x, y = y, M =M)
fit <- fitGpPseudoInputs(y, x, M,
                         xbar = xbar,
                         pseudoInputSelection = "none",
                         cond = 1e-14)
plotGpPseudoInputs(fit,cond = 1e-7)
mse(fit,cond = 1e-7)

save(fit, file = "fit_givenPI_PS.Rdata")

end <- Sys.time()
end-start


start <- Sys.time()

#Maximize Marginal likelihood
#-----------------------------
fit <- fitGpPseudoInputs(y, x, M,
                         pseudoInputSelection = "ml",
                         cond = 1e-14)
plotGpPseudoInputs(fit,cond = 1e-7)
mse(fit,cond = 1e-7)

save(fit, file = "fit_maxML.Rdata")

end <- Sys.time()
end-start


start <- Sys.time()

#Minimize KL divergence
#-----------------------
fit <- fitGpPseudoInputs(y, x, M, 
                         pseudoInputSelection = "kl",
                         cond = 1e-14)
plotGpPseudoInputs(fit,cond = 1e-6)
mse(fit,cond = 1e-6)

save(fit, file = "fit_minKL.Rdata")

end <- Sys.time()
end-start


start <- Sys.time()

#Fit with full data
#--------------------
library(DiceKriging)
model <- km(~1, design = data.frame(x=x), covtype = "gauss",
            response = data.frame(y=y), nugget.estim = T)

xPred <- seq(0,1,length = 100)
p <- predict(model, data.frame(x=xPred), "UK")

par(mfrow=c(1,1))
plot(x, y)

lines(xPred, p$mean, col = "red")
lines(xPred, p$lower95, lty =2)
lines(xPred, p$upper95, lty =2)


#Save plots of results 
#---------------------
listFit <- c("Fixed pseudo inputs - grid" = "fit_givenPI_grid.Rdata",
             "Fixed pseudo inputs - SRS" =  "fit_givenPI_SRS.Rdata",
             "Fixed pseudo inputs - preferential sampling" =   "fit_givenPI_PS.Rdata",
             "Maximization marginal likelihood" = "fit_maxML.Rdata",
             "Minimization KL divergence" = "fit_minKL.Rdata")

library(DiceKriging)
model <- km(~1, design = data.frame(x=x), covtype = "gauss",
            response = data.frame(y=y), nugget.estim = T)

xPred <- seq(0,1,length = 100)
p <- predict(model, data.frame(x=xPred), "UK")


for(i in 1:length(listFit)) {
  load(file = listFit[[i]])
  
  plotGpPseudoInputs(fit, cond = 1e-7, 
                     main = names(listFit)[i])
  
  lines(xPred, p$mean, col = "darkviolet")
  lines(xPred, p$lower95, lty =2, col = "darkviolet")
  lines(xPred, p$upper95, lty =2, col = "darkviolet")
  
  #lines(x[order(x)], yGp[order(x)], col = "black")
  dev.copy(pdf, file = gsub(".Rdata",".pdf",listFit[[i]]), 
           height = 6, width =8)
  dev.off()
  
  
  cat("-------------------------------------------------")
  cat("\n\n")
  cat(names(listFit)[i])
  cat("\n")
  cat(paste0("sigma noise: ", fit$sigma_noiseOpt, "\n"))
  cat(paste0("sigma K: ", fit$sigma_noiseOpt, "\n"))
  cat(paste0("rho: ", fit$rhoOpt, "\n"))
  cat(paste0("MSE: ", mse(fit, cond = 1e-7), "\n"))
  cat("\n")
  cat("-------------------------------------------------")
  cat("\n")
  
}


####################################
# Then, simulations with several m #
####################################

#Number of Pseudo-inputs
vectorM <- c(5,7,10,12,15,20,30,40,50)

cl <- makeCluster(2)
registerDoParallel(cl)

start <- Sys.time()

results <- foreach(i = 1:length(vectorM), 
                   .packages = c("DiceKriging", "GenSA")) %dopar% {
  
  source("fitPseudoInputsGP.R")
  source("summaryFitPseudoInputsGP.R")
  
  M <- vectorM[i]
  
  fit_givenPI_grid <- fitGpPseudoInputs(y, x, M,
                           xbar = seq(0,1, length = M),
                           pseudoInputSelection = "none",
                           cond = 1e-14)
  
  fit_givenPI_SRS <- fitGpPseudoInputs(y, x, M,
                           xbar = sample(x = x, size = M),
                           pseudoInputSelection = "none",
                           cond = 1e-14)
  
  xbar <- preferPseudoInputs(x = x, y = y, M = M)
  fit_givenPI_PS <- fitGpPseudoInputs(y, x, M,
                           xbar = xbar,
                           pseudoInputSelection = "none",
                           cond = 1e-14)
  
  fit_maxML <- fitGpPseudoInputs(y, x, M,
                           pseudoInputSelection = "ml",
                           cond = 1e-14)
  
  fit_minKL <- fitGpPseudoInputs(y, x, M, 
                           pseudoInputSelection = "kl",
                           cond = 1e-14)
  
  resultTemp <- list(m = M,
                     fit_givenPI_grid = fit_givenPI_grid, 
                     fit_givenPI_SRS = fit_givenPI_SRS, 
                     fit_givenPI_PS = fit_givenPI_PS,
                     fit_maxML = fit_maxML, 
                     fit_minKL= fit_minKL)

  resultTemp
}

save(results, file = "simulation_effect_m_2.Rdata")

end <- Sys.time()

#####################################
# Extract information of simulation #
#####################################

load(file = "simulation_effect_m_2.Rdata")

mseResults <- data.frame(m = sapply(results, FUN = "[[", "m"),
                         fit_givenPI_grid = rep(NA, length(results)), 
                         fit_givenPI_SRS = rep(NA, length(results)), 
                         fit_givenPI_PS = rep(NA, length(results)),
                         fit_maxML = rep(NA, length(results)), 
                         fit_minKL= rep(NA, length(results)))


for(fitString in setdiff(names(results[[1]]),"m")) {
  
  pdf(paste0(fitString, "_2.pdf"), onefile = T)
  
  for(i in 1:length(results)) {
    resultsIter <- results[[i]]
    plotGpPseudoInputs(resultsIter[[fitString]], cond = 1e-5)
    mseResults[i,fitString] <- mse(resultsIter[[fitString]], cond = 1e-5)
  }
  
  dev.off()
  
}


colors = rainbow(5)

plot(mseResults$m, log(mseResults$fit_givenPI_grid, base = 10), 
     main ="Prediction error - Wiggly function", 
     type = "b",  lwd =2, 
     ylab = expression('log'[10]*'(Prediction error)'), xlab = "Number of pseudo-inputs", 
     col = colors[1], ylim = c(-3,-1))
points(mseResults$m, log(mseResults$fit_givenPI_SRS, base = 10), 
       type = "b", lwd =2,  col = colors[2])
points(mseResults$m, log(mseResults$fit_givenPI_PS, base = 10), 
       type = "b", lwd =2, col = colors[3])
points(mseResults$m, log(mseResults$fit_maxML, base = 10), 
       type = "b", lwd =2, col = colors[4])
points(mseResults$m, log(mseResults$fit_minKL, base = 10), 
       type = "b", lwd =2, col = colors[5])

legend("topright", col = colors, pch = 1, lwd =2, 
       legend = c("Grid","Random sampling", "Preferential sampling",
                  "Max Marginal Likelihood", "Min KL divergence"))

dev.copy(pdf, file = "summary_effect_m_2.pdf", height = 6, width = 8)
dev.off()
