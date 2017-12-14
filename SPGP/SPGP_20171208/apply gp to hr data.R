#Sparse Peudo-input Gaussian Process 
#Author: Hengrui Luo (luo.619@osu.edu) and Giovanni Nattino (nattino.1@osu.edu)

setwd("C:/Users/Giovanni/Google Drive/PhD/Classes/STAT 8810/final project/R code")

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

#Plot of the selected portion of HR records
par(mfrow = c(4,2), mar = c(5,5,3,0))
for(i in 1:4) {
  lettersSelected <- LETTERS[c(2*(i-1)+1,2*(i-1)+2)]
  
  dat <- eval(parse(text = paste0("dat", lettersSelected[1])))
  dat1 <- dat[eval(parse(text = eval(parse(text = paste0("selection", lettersSelected[1]))))), ]
  
  dat <- eval(parse(text = paste0("dat", lettersSelected[2])))
  dat2 <- dat[eval(parse(text = eval(parse(text = paste0("selection", lettersSelected[2]))))), ]

  cat("\n") 
  cat(nrow(dat1))  
  cat("\n") 
  cat(nrow(dat2))  
  cat("\n") 
  dat1$time <- dat1$time - min(dat1$time)
  dat2$time <- dat2$time - min(dat2$time)
  
  plot(dat1$time, dat1$hr, type = "l", 
       main = paste0("Exercise ",i),
       xlab = "time", ylab = "bpm",
       ylim = c(mean(dat1$hr)-20.2,mean(dat1$hr)+14.8))
  plot(dat2$time, dat2$hr, type = "l", 
       main = paste0("Recovery after exercise ",i),
       xlab = "time", ylab = "bpm",
       ylim = c(mean(dat2$hr)-20.2,mean(dat2$hr)+14.8))
}

dev.copy(pdf, file = "HR_all.pdf",
         height = 11, width = 8.5)
dev.off()


M <- 40

library(doSNOW)
library(parallel)
library(doRNG)

cl <- makeCluster(8)
registerDoSNOW(cl)

start <- Sys.time()

pb <- txtProgressBar(max = 8, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

resultTemp <- foreach(iter=1:8, 
                      .packages = c("GenSA"), .options.snow = opts,
                      .export = c(paste0("selection",LETTERS[1:8]),
                                  paste0("dat",LETTERS[1:8]))) %dorng% {
  
  source("simulateGP.R")
  source("fitPseudoInputsGP.R")
  source("summaryFitPseudoInputsGP.R")
                        
  letter <- LETTERS[iter]
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
                           cond = 1e-10, maxTime = 8*60*60)
  save(fit, file = paste0("fit_givenPI_grid_HR_",letter,".Rdata"))

  #Marginal likelihood
  set.seed(123)
  cat("ML: ")
  start <- Sys.time()
  fit <- fitGpPseudoInputs(dat$hr, dat$time, M = M,
                           sigma_noiseStart = 1, sigma_KStart = 1, 
                           pseudoInputSelection = "ml",
                           cond = 1e-10, maxTime = 8*60*60)
  save(fit, file = paste0("fit_maxML_HR_",letter,".Rdata"))

  #Minimize KL divergence
  set.seed(123)
  cat("KL: ")
  start <- Sys.time()
  fit <- fitGpPseudoInputs(dat$hr, dat$time, M = M, 
                           sigma_noiseStart = 1, sigma_KStart = 1,
                           pseudoInputSelection = "kl",
                           cond = 1e-10, maxTime = 8*60*60)
  save(fit, file = paste0("fit_minKL_HR_",letter,".Rdata"))

} 

stopCluster(cl)

#Save plots of results 
#---------------------
par()

listFit <- c("Fixed pseudo inputs - grid" = "fit_givenPI_grid_HR.Rdata",
             "Maximization marginal likelihood" = "fit_maxML_HR.Rdata",
             "Minimization KL divergence" = "fit_minKL_HR.Rdata")


mseResults <- data.frame(letter = LETTERS[1:8],
                         fit_givenPI_grid_HR = rep(NA, 8),
                         fit_maxML_HR = rep(NA, 8),
                         fit_minKL_HR = rep(NA, 8))

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
    
    if((which(LETTERS==letter) %% 2) ==1) {
      startTitle <- "Exercise "  
      xAxisOrigScale <- seq(min(dat$time),max(dat$time)+10, by = 50)
      xAxisNewScale <- seq(0, 200, by = 50)
      
    } else {
      startTitle <- "Recovery after exercise "  
      xAxisOrigScale <- seq(min(dat$time), max(dat$time)+10, by = 100)
      xAxisNewScale <- seq(0, 400, by = 100)
      
    }
    
    plotGpPseudoInputs(fit, cond = 1e-5, 
                       main = paste0(startTitle, 
                                     ceiling(which(LETTERS==letter) / 2),
                                     "\n",
                                     names(listFit)[i]), 
                       xPred = timePred,
                       xlab = "time",
                       ylab = "bpm", xaxt = "n", type = "l")
    axis(side = 1, 
         at = xAxisOrigScale,
         labels = xAxisNewScale)
    
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
    cat("\n")
    cat("-------------------------------------------------")
    cat("\n")
    
    mseResults[mseResults$letter == letter,
               gsub(".Rdata", "", listFit[[i]])] <- mse(fit, cond = 1e-5)
    
  }
}

save(mseResults, file = "mseResults.Rdata")


#Distributions of pseudo-inputs
#------------------------------

listFitMod <- c("Maximization marginal likelihood" = "fit_maxML_HR.Rdata",
                 "Minimization KL divergence" = "fit_minKL_HR.Rdata")


for(letter in LETTERS[1:8]) {
  
  dat <- eval(parse(text = paste0("dat", letter)))
  dat <- dat[eval(parse(text = eval(parse(text = paste0("selection", letter))))), ]
  timePred <- seq(from = min(dat$time), to = max(dat$time), length = 1000)
  
  
  for(i in 1:length(listFitMod)) {
    
    load(file = gsub(".Rdata", 
                     paste0("_",letter,".Rdata"), 
                     listFitMod[[i]]))
    
    if((which(LETTERS==letter) %% 2) ==1) {
      startTitle <- "Exercise "  
      xAxisOrigScale <- seq(min(dat$time),max(dat$time)+10, by = 50)
      xAxisNewScale <- seq(0, 200, by = 50)
      
    } else {
      startTitle <- "Recovery after exercise "  
      xAxisOrigScale <- seq(min(dat$time), max(dat$time)+10, by = 100)
      xAxisNewScale <- seq(0, 400, by = 100)
      
    }
    
    hist(fit$xbarOpt,
         main = paste0(startTitle, 
                            ceiling(which(LETTERS==letter) / 2),
                            "\n",
                            names(listFitMod)[i]),
                  xlab = "time",
                  ylab = "Number of pseudo-inputs", xaxt = "n")
    axis(side = 1, 
         at = xAxisOrigScale,
         labels = xAxisNewScale)
    

    dev.copy(pdf, file = gsub(".Rdata", 
                              paste0("_",letter,"_piLocation.pdf"), 
                              listFitMod[[i]]), 
             height = 6, width =8)
    dev.off()
    
  }
}
