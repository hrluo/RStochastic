#Sparse Peudo-input Gaussian Process 
#Author: Hengrui Luo (luo.619@osu.edu) and Giovanni Nattino (nattino.1@osu.edu)

library(GenSA)
library(DiceKriging)
library(doSNOW)
library(doParallel)
library(doRNG)

###############################
# Source Files with Functions #
###############################

source("simulateGP.R")
source("fitPseudoInputsGP.R")
source("summaryFitPseudoInputsGP.R")


########################################################
# Simulation 2: Several draws from GP, compare methods #
########################################################

# Number of observed points
N <- 100

# Number of pseudo-inputs
M <- 10

#SD kernel
sigma_K <- 1

#SD noise 
sigmaNoise <- 0.05

#Correlation
rho <-  .1e-8

#Create cluster to run in parallel
cl <- makeCluster(2)
registerDoSNOW(cl)

# Number of simulations
nDraws <- 50

#Reproducibility
set.seed(123456)

start <- Sys.time()

pb <- txtProgressBar(max = nDraws, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results <- foreach(i=1:nDraws, .packages = c("DiceKriging", "GenSA"),
                   .options.snow = opts) %dorng% {
  
  source("simulateGP.R")
  source("fitPseudoInputsGP.R")
  source("summaryFitPseudoInputsGP.R")
  
  #Draw a GP
  simulation <- simulateGP_GaussianCorr(p = 1, 
                                        n = N, 
                                        theta = -log(rho), 
                                        sigma2 = sigma_K)
  x <- as.vector(simulation$X)
  yGp <- as.vector(simulation$Y)
  #Add noise
  y <- yGp + rnorm(length(yGp), sd = sigmaNoise)
  
  fit_givenPI_grid <- fitGpPseudoInputs(y, x, M,
                           xbar = seq(0,1, length = M),
                           pseudoInputSelection = "none",
                           cond = 1e-14, maxTime = 60*60)
  
  fit_givenPI_SRS <- fitGpPseudoInputs(y, x, M,
                           xbar = sample(x = x, size = M),
                           pseudoInputSelection = "none",
                           cond = 1e-14, maxTime = 60*60)
  
  xbar <- preferPseudoInputs(x = x, y = y, M =M)
  fit_givenPI_PS <- fitGpPseudoInputs(y, x, M,
                           xbar = xbar,
                           pseudoInputSelection = "none",
                           cond = 1e-14, maxTime = 60*60)
  
  fit_maxML <- fitGpPseudoInputs(y, x, M,
                           pseudoInputSelection = "ml",
                           cond = 1e-14, maxTime = 60*60)
  
  fit_minKL <- fitGpPseudoInputs(y, x, M, 
                           pseudoInputSelection = "kl",
                           cond = 1e-14, maxTime = 60*60)
  
  resultTemp <- list(fit_givenPI_grid = fit_givenPI_grid, 
                     fit_givenPI_SRS = fit_givenPI_SRS, 
                     fit_givenPI_PS = fit_givenPI_PS,
                     fit_maxML = fit_maxML, 
                     fit_minKL= fit_minKL)
  
  resultTemp
}

save(results, file = "simulation_effect_draws_2.Rdata")

end <- Sys.time()

#####################################
# Extract information of simulation #
#####################################

load(file = "simulation_effect_draws_2.Rdata")
load(file = "C:/Users/Giovanni/Google Drive/PhD/Classes/STAT 8810/final project/results/Result effect draws - 2/simulation_effect_draws_2.Rdata")


mseResults <- data.frame(fit_givenPI_grid = rep(NA, length(results)), 
                         fit_givenPI_SRS = rep(NA, length(results)), 
                         fit_givenPI_PS = rep(NA, length(results)),
                         fit_maxML = rep(NA, length(results)), 
                         fit_minKL= rep(NA, length(results)))


for(fitString in names(results[[1]])) {
  
  pdf(paste0(fitString, "_2.pdf"), onefile = T)
  
  for(i in 1:length(results)) {
      resultsIter <- results[[i]]
      plotGpPseudoInputs(resultsIter[[fitString]], cond = 1e-6)
      mseResults[i,fitString] <- mse(resultsIter[[fitString]], cond = 1e-6)
  }
  
  dev.off()
 
}

colors = rainbow(5)
boxplot(log(mseResults$fit_givenPI_grid, base = 10), 
        log(mseResults$fit_givenPI_SRS, base = 10),
        log(mseResults$fit_givenPI_PS, base = 10),
        log(mseResults$fit_maxML, base = 10),
        log(mseResults$fit_minKL, base = 10), col = colors,
        names = c("Grid", "RS", "PS",
                  "Max ML", "Min KL"),
        ylab = expression('log'[10]*'(Prediction error)'),
        main = "Prediction error - Wiggly functions")

dev.copy(pdf, file = "summary_effect_draws_2.pdf", height = 6, width = 8)
dev.off()

