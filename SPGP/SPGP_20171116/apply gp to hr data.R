#Sparse Peudo-input Gaussian Process 
#Author: Hengrui Luo (luo.619@osu.edu) and Giovanni Nattino (nattino.1@osu.edu)

setwd("C:/Users/Giovanni/Google Drive/PhD/Classes/STAT 8810/final project/S1_Dataset")

library(GenSA)

###############################
# Source Files with Functions #
###############################

source("simulateGP.R")
source("fitPseudoInputsGP.R")
source("summaryFitPseudoInputsGP.R")


datA <- read.table(file = "Exercise A.txt")
datB <- read.table(file = "Exercise B.txt")
datC <- read.table(file = "Exercise C.txt")
datD <- read.table(file = "Exercise D.txt")
datE <- read.table(file = "Exercise E.txt")
datF <- read.table(file = "Exercise F.txt")
datG <- read.table(file = "Exercise G.txt")
datH <- read.table(file = "Exercise H.txt")

(names(datA) <- names(datB) <- names(datC) <- 
   names(datD) <- names(datE) <- names(datF) <- 
    names(datG) <- names(datH) <- c("time", "hr"))


par(mfrow = c(4,2), mar = c(3,3,2,0))
plot(datA$time, datA$hr, type = "l")
plot(datB$time, datB$hr, type = "l")
plot(datC$time, datC$hr, type = "l")
plot(datD$time, datD$hr, type = "l")
plot(datE$time, datE$hr, type = "l")
plot(datF$time, datF$hr, type = "l")
plot(datG$time, datG$hr, type = "l")
plot(datH$time, datH$hr, type = "l")
# dev.copy(pdf, file = "HR_all.pdf", 
#          height = 11, width = 8.5)
# dev.off()

# Example on exercise A 
#-----------------------
dat <- datA
selectionA <- dat$time>=200 & dat$time<=400 
dat <- dat[selectionA, ]
timePred <- seq(from = min(dat$time), to = max(dat$time), length = 1000)

#Fit full model
library(DiceKriging)
model <- km(~1, design = data.frame(x=dat$time), 
            response = data.frame(y=dat$hr), nugget.estim = T)
p <- predict(model, data.frame(x=timePred), "UK")

#Plot estimate full model
par(mfrow=c(1,1))
plot(dat$time, dat$hr, type = "l")
lines(timePred, p$mean, col = "red")
lines(timePred, p$lower95, lty =2)
lines(timePred, p$upper95, lty =2)

# Apply pseudo inputs methods
#----------------------------

selectionA <- "dat$time>=200 & dat$time<=400" 
selectionB <- "dat$time>=100 & dat$time<=500" 
selectionC <- "dat$time>=150 & dat$time<=350" 
selectionD <- "dat$time>=100 & dat$time<=500" 
selectionE <- "dat$time>=150 & dat$time<=350" 
selectionF <- "dat$time>=100 & dat$time<=500"
selectionG <- "dat$time>=110 & dat$time<=310" 
selectionH <- "dat$time>=100 & dat$time<=500" 

M <- 20

for(letter in LETTERS[1:8]) {
  
  cat(letter)
  cat("\n")
  
  dat <- eval(parse(text = paste0("dat", letter)))
  dat <- dat[eval(parse(text = eval(parse(text = paste0("selection", letter))))), ]

  #Provide random pseudo-inputs: grid 
  set.seed(123)
  cat("Grid: ")
  start <- Sys.time()
  fit <- fitGpPseudoInputs(dat$hr, dat$time,
                           sigma_noiseStart = 1, sigma_KStart = 1, 
                           xbar = seq(0,1, length = M),
                           pseudoInputSelection = "none",
                           cond = 1e-10, maxTime = 2*60*60)
  save(fit, file = paste0("fit_givenPI_grid_HR_",letter,".Rdata"))
  end <- Sys.time()
  
  cat(end-start)
  cat("\n")
  
  #Marginal likelihood
  set.seed(123)
  cat("ML: ")
  start <- Sys.time()
  fit <- fitGpPseudoInputs(dat$hr, dat$time, M = M,
                           sigma_noiseStart = 1, sigma_KStart = 1, 
                           pseudoInputSelection = "ml",
                           cond = 1e-10, maxTime = 2*60*60)
  save(fit, file = paste0("fit_maxML_HR_",letter,".Rdata"))
  end <- Sys.time()
  
  cat(end-start)
  cat("\n")
  
  #Minimize KL divergence
  set.seed(123)
  cat("KL: ")
  start <- Sys.time()
  fit <- fitGpPseudoInputs(dat$hr, dat$time, M = M, 
                           sigma_noiseStart = 1, sigma_KStart = 1,
                           pseudoInputSelection = "kl",
                           cond = 1e-10, maxTime = 2*60*60)
  save(fit, file = paste0("fit_minKL_HR_",letter,".Rdata"))
  end <- Sys.time()
  
  cat(end-start)
  cat("\n")
  
} 


#Save plots of results 
#---------------------
par()

listFit <- c("Fixed pseudo inputs - grid" = "fit_givenPI_grid_HR.Rdata",
             "Maximization marginal likelihood" = "fit_maxML_HR.Rdata",
             "Minimization KL divergence" = "fit_minKL_HR.Rdata")


library(DiceKriging)

for(letter in LETTERS[1:8]) {
  
  dat <- eval(parse(text = paste0("dat", letter)))
  dat <- dat[eval(parse(text = eval(parse(text = paste0("selection", letter))))), ]
  timePred <- seq(from = min(dat$time), to = max(dat$time), length = 1000)
  
  model <- km(~1, design = data.frame(x=dat$time), 
              response = data.frame(y=dat$hr), nugget.estim = T)
  p <- predict(model, data.frame(x=timePred), "UK")

  
  for(i in 1:length(listFit)) {
    
    load(file = gsub(".Rdata", 
                     paste0("_",letter,".Rdata"), 
                     listFit[[i]]))
    
    plotGpPseudoInputs(fit, cond = 1e-5, 
                       main = names(listFit)[i], 
                       xPred = timePred)
    
    lines(timePred, p$mean, col = "darkviolet")
    lines(timePred, p$lower95, lty =2, col = "darkviolet")
    lines(timePred, p$upper95, lty =2, col = "darkviolet")
    
    dev.copy(pdf, file = gsub(".Rdata", 
                              paste0("_",letter,".pdf"), 
                              listFit[[i]]), 
             height = 6, width =8)
    dev.off()
    
    
    cat("-------------------------------------------------")
    cat("\n\n")
    cat(names(listFit)[i])
    cat("\n")
    cat(paste0("sigma noise: ", fit$sigma_noiseOpt, "\n"))
    cat(paste0("sigma K: ", fit$sigma_noiseOpt, "\n"))
    cat(paste0("rho: ", fit$rhoOpt, "\n"))
    cat(paste0("MSE: ", mse(fit, cond = 1e-5), "\n"))
    cat("\n")
    cat("-------------------------------------------------")
    cat("\n")
    
  }
}


