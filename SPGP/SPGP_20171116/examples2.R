#Sparse Peudo-input Gaussian Process 
#Author: Hengrui Luo (luo.619@osu.edu) and Giovanni Nattino (nattino.1@osu.edu)

library(GenSA)

###############################
# Source Files with Functions #
###############################

source("simulateGP.R")
source("fitPseudoInputsGP.R")
source("summaryFitPseudoInputsGP.R")


#############
# Example 2 #
#############

#Fit a realization of the GP with a GP with pseudo-inputs.
set.seed(12)

# Number of observed points
N <- 100

#SD kernel
sigma_K <- 3

#Correlation
rho <- .001

#Draw of a GP
simulation <- simulateGP_GaussianCorr(p = 1, 
                                      n = N, 
                                      theta = -log(rho), 
                                      sigma2 = sigma_K)
x <- as.vector(simulation$X)
yGp <- as.vector(simulation$Y)

y <- yGp + rnorm(length(yGp), sd = .001)

#Number of Pseudo-inputs
M <- 5

start <- Sys.time()

#Provide random pseudo-inputs: grid 
#----------------------------------
fit <- fitGpPseudoInputs(y, x, M,
                         xbar = seq(0,1, length = M),
                         pseudoInputSelection = "none",
                         cond = 1e-14)
plotGpPseudoInputs(fit, cond = 1e-7)
mse(fit,cond = 1e-7)

save(fit, file = "fit_givenPI_grid_ex2.Rdata")

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

save(fit, file = "fit_givenPI_SRS_ex2.Rdata")

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
mse(fit,cond = 1e-10)

save(fit, file = "fit_givenPI_PS_ex2.Rdata")

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

save(fit, file = "fit_maxML_ex2.Rdata")

end <- Sys.time()
end-start


start <- Sys.time()

#Minimize KL divergence
#-----------------------
fit <- fitGpPseudoInputs(y, x, M, 
                         pseudoInputSelection = "kl",
                         cond = 1e-14)
plotGpPseudoInputs(fit,cond = 1e-7)
mse(fit,cond = 1e-7)

save(fit, file = "fit_minKL_ex2.Rdata")

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
listFit <- c("Fixed pseudo inputs - grid" = "fit_givenPI_grid_ex2.Rdata",
             "Fixed pseudo inputs - SRS" =  "fit_givenPI_SRS_ex2.Rdata",
             "Fixed pseudo inputs - preferential sampling" =   "fit_givenPI_PS_ex2.Rdata",
             "Maximization marginal likelihood" = "fit_maxML_ex2.Rdata",
             "Minimization KL divergence" = "fit_minKL_ex2.Rdata")

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
  
  lines(x, yGp, col = "black")
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

