# packages
library(MASS)
library(Matrix)
library(glasso)
library(matrixcalc)
library(cvTools)
library(igraph)

library(tidyverse)
library(reshape2)

# Sparse Scenario

n = 100
p = 20
set.seed(41118)

sparse = matrix(0, nrow=p, ncol=p)
sparse[1, 1] = 1
for (i in 2:nrow(sparse)){
  sparse[i,i] = 1
  sparse[i, i-1] = 0.5
  sparse[i-1, i] = 0.5
}

sparsecov = solve(sparse)
# cov2cor(sparsecov) 

# Try alternative data generation from https://www4.stat.ncsu.edu/~reich/BigData/code/glasso.html
sparsedata <- matrix(rnorm(n*p),n,p)%*%chol(sparsecov)

sparse_empcov = cov(sparsedata)

# Grouped Scenario

g1 = qr.Q(qr(matrix(runif((p/2)^2, -.5, .5), (p/2)))) # https://stats.stackexchange.com/questions/215497/how-to-create-an-arbitrary-covariance-matrix
g1cov <- crossprod(g1, g1*((p/2):1))
g1prec <- crossprod(g1, g1/((p/2):1))

g2 = qr.Q(qr(matrix(runif((p/2)^2, -.5, 5), (p/2))))
g2cov <- crossprod(g2, g2*((p/2):1))
g2prec <- crossprod(g2, g2/((p/2):1))

group = matrix(0, nrow=p, ncol=p)
group[1:(p/2),1:(p/2)] = g1prec
group[(p/2 + 1):p,(p/2 + 1):p] = g2prec

group = group*5 # Make values in precision matrix larger

groupcov = as.matrix(solve(group))

cov2cor(groupcov)

groupdata <- matrix(rnorm(n*p),n,p)%*%chol(groupcov)

group_empcov = cov(groupdata)

# Dense Scenario

pos = qr.Q(qr(matrix(runif(p^2, -1, 1), p)))
dense <- crossprod(pos, pos/(p:1))*10
densecov = solve(dense)
cov2cor(densecov)

densedata <- matrix(rnorm(n*p),n,p)%*%chol(densecov)

dense_empcov = cov(densedata)

# Functions adapated from https://github.com/czarrar/causality_primer/blob/master/journal/2014-07-12_cross_validation_graphical_lasso.Rmd

log_likelihood <- function(precision, emp_cov) {
  p      <- nrow(precision)
  logdet <- determinant(precision, logarithm=T)$modulus
  loglik <- logdet - sum(diag(emp_cov %*% precision))
  return(as.numeric(loglik))
}

glasso_cv <- function(ts, k=5, rholist=NULL, verbose=T) {
  library(glasso)
  library(cvTools)
  
  if (is.null(rholist)) {
    S       <- cov(ts)
    # Now make the list of rhos
    rholist <- seq(0, max(abs(S)), length=300)
    if (nrow(ts)<=nrow(S)) {
      rholist <- (rholist+0.1)
    }
  }
  n     <- nrow(ts) 
  folds <- cvFolds(n, k, type="consecutive")
  
  loglikes <- sapply(1:k, function(ki) {
    if (verbose) cat("Fold ", ki, "\n")
    S_train <- cov(ts[folds$which!=ki,])
    S_test  <- cov(ts[folds$which==ki,])
    GLP     <- glassopath(S_train, rholist, trace=0)
    loglike <- apply(GLP$wi, 3, function(P_train) log_likelihood(P_train, S_test))
    # loglike = numeric()
    # for (i in 1:length(rholist)){
    #   loglike[i] = log_likelihood(GLP$wi[,,i], S_test)
    # }
    loglike
  })
  
  ind     <- which.max(rowMeans(loglikes))
  rhomax  <- rholist[ind]
  S       <- cov(ts)
  a       <- glasso(S, rhomax)
  a$rhomax <- rhomax
  
  return(list(a, rowMeans(loglikes), rholist))
}

# Add BIC

log_likelihood_bic <- function(precision, emp_cov) {
  p      <- nrow(precision)
  logdet <- determinant(precision, logarithm=T)$modulus
  loglik <- logdet - sum(diag(emp_cov %*% precision))
  return(as.numeric(loglik))
}

glasso_bic <- function(ts, rholist=NULL, verbose=T) {
  library(glasso)
  library(cvTools)
  
  if (is.null(rholist)) {
    S       <- cov(ts)
    rholist <- seq(0, max(abs(S)), length=200)
    if (nrow(ts)<=nrow(S)) {
      rholist <- (rholist+0.1)
    }
  }
  n     <- nrow(ts)
  dim <- nrow(S)
  BICs <- numeric()
  loglike = numeric()
  for (i in 1:length(rholist)){
    GLP <- glasso(S, rholist[i], trace=FALSE)
    p_off_d <- sum(GLP$wi!=0 & col(S)<row(S))
    loglike[i] = log_likelihood_bic(GLP$wi, S)
    BICs[i] = -n*(loglike[i]) + p_off_d*log(n)
  }
  
  ind     <- which.min(BICs)
  rhomin  <- rholist[ind]
  S       <- cov(ts)
  a       <- glasso(S, rhomin)
  a$rhomin <- rhomin
  
  return(list(a, BICs, rholist))
}

# Extended BIC

glasso_bicext <- function(ts, rholist=NULL, verbose=T) {
  library(glasso)
  library(cvTools)
  
  if (is.null(rholist)) {
    S       <- cov(ts)
    rholist <- seq(0, max(abs(S)), length=200)
    if (nrow(ts)<=nrow(S)) {
      rholist <- (rholist+0.1)
    }
  }
  n     <- nrow(ts)
  dim <- nrow(S)
  BICs <- numeric()
  loglike = numeric()
  for (i in 1:length(rholist)){
    GLP <- glasso(S, rholist[i], trace=FALSE)
    p_off_d <- sum(GLP$wi!=0 & col(S)<row(S))
    loglike[i] = log_likelihood_bic(GLP$wi, S)
    BICs[i] = -n*(loglike[i]) + p_off_d*log(n) + 4*p_off_d*1*log(dim)
  }
  
  ind     <- which.min(BICs)
  rhomin  <- rholist[ind]
  S       <- cov(ts)
  a       <- glasso(S, rhomin)
  a$rhomin <- rhomin
  
  return(list(a, BICs, rholist))
}

#### Set Rounding Here

rounding = 1

# Generate graphs

sparse_round = round(sparse, rounding)
group_round = round(group, rounding)
dense_round = round(dense, rounding)

sparsetrue = ifelse(sparse!=0 & row(sparse)!=col(sparse),1,0)
grouptrue = ifelse(group!=0 & row(group)!=col(group),1,0)
densetrue = ifelse(dense!=0 & row(dense)!=col(dense),1,0)

sparsetrue_round = ifelse(sparse_round!=0 & row(sparse_round)!=col(sparse_round),1,0)
grouptrue_round = ifelse(group_round!=0 & row(group_round)!=col(group_round),1,0)
densetrue_round = ifelse(dense_round!=0 & row(dense_round)!=col(dense_round),1,0)

sparsegraph_t = graph_from_adjacency_matrix(sparsetrue, mode="undirected")
groupgraph_t = graph_from_adjacency_matrix(grouptrue, mode="undirected")
densegraph_t = graph_from_adjacency_matrix(densetrue, mode="undirected")

sparsegraph_t_r = graph_from_adjacency_matrix(sparsetrue_round, mode="undirected")
groupgraph_t_r = graph_from_adjacency_matrix(grouptrue_round, mode="undirected")
densegraph_t_r = graph_from_adjacency_matrix(densetrue_round, mode="undirected")

### Plot Sparse
#name nodes
V(sparsegraph_t)$name <- paste0(1:20)
#Get the coordinates of the Nodes
sparseCoords <- layout_with_fr(sparsegraph_t) %>% 
  as_tibble %>%
  bind_cols(data_frame(names = names(V(sparsegraph_t))))

sparsecoords_rescaled <- sapply(sparseCoords[-3], function(x) -1+((x-min(x))*2)/diff(range(x)))
rownames(sparsecoords_rescaled) <- sparseCoords$names

sparsegraph_t%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)

sparsegraph_t_r%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)

# sparsegraph_t%>%
#   plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F)


### Plot Group
#name nodes
V(groupgraph_t)$name <- paste0(1:20)
#Get the coordinates of the Nodes
groupCoords <- layout_with_fr(groupgraph_t) %>% 
  as_tibble %>%
  bind_cols(data_frame(names = names(V(groupgraph_t))))

groupcoords_rescaled <- sapply(groupCoords[-3], function(x) -1+((x-min(x))*2)/diff(range(x)))
rownames(groupcoords_rescaled) <- groupCoords$names

# par(mfrow=c(1,2))
groupgraph_t%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)

groupgraph_t_r%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)
# par(mfrow=c(1,1))

# groupgraph_t%>% 
#   plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F)


### Plot Dense
#name nodes
V(densegraph_t)$name <- paste0(1:20)
#Get the coordinates of the Nodes
denseCoords <- layout_with_fr(densegraph_t) %>% 
  as_tibble %>%
  bind_cols(data_frame(names = names(V(densegraph_t))))

densecoords_rescaled <- sapply(denseCoords[-3], function(x) -1+((x-min(x))*2)/diff(range(x)))
rownames(densecoords_rescaled) <- denseCoords$names

# par(mfrow=c(1,2))
densegraph_t%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)

densegraph_t_r%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)
# par(mfrow=c(1,1))
# densegraph_t%>% 
#   plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F)

### Set up simulation files

p=20
samps = c(15, 20, 35, 50, 75, 100, 150, 200, 250, 350, 500, 750, 1000, 1250)
dthetacv = list()
gthetacv = list()
sthetacv = list()
dthetabic = list()
gthetabic = list()
sthetabic = list()
dthetabicext = list()
gthetabicext = list()
sthetabicext = list()
cvmat = matrix(NA, nrow=3, ncol=length(samps))
bicmat = matrix(NA, nrow=3, ncol=length(samps))
bicextmat = matrix(NA, nrow=3, ncol=length(samps))
cvl1 = matrix(NA, nrow=3, ncol=length(samps))
bicl1 = matrix(NA, nrow=3, ncol=length(samps))
bicextl1 = matrix(NA, nrow=3, ncol=length(samps))
cvspar = matrix(NA, nrow=3, ncol=length(samps))
bicspar = matrix(NA, nrow=3, ncol=length(samps))
bicextspar = matrix(NA, nrow=3, ncol=length(samps))
cvspar_round = matrix(NA, nrow=3, ncol=length(samps))
bicspar_round = matrix(NA, nrow=3, ncol=length(samps))
bicextspar_round = matrix(NA, nrow=3, ncol=length(samps))
maxloglik = matrix(NA, nrow=6, ncol=length(samps))
minbic = matrix(NA, nrow=6, ncol=length(samps))
minbicext = matrix(NA, nrow=6, ncol=length(samps))

# Iterate for different sample sizes

for (i in 1:length(samps)){
  sdat = matrix(rnorm(samps[i]*p),samps[i],p)%*%chol(sparsecov)
  gdat = matrix(rnorm(samps[i]*p),samps[i],p)%*%chol(groupcov)
  ddat = matrix(rnorm(samps[i]*p),samps[i],p)%*%chol(densecov)
  scov = cov(sdat)
  gcov = cov(gdat)
  dcov = cov(ddat)
  dmodcv = glasso_cv(ddat, k=ifelse(samps[i]<=20, 5, 10))
  gmodcv = glasso_cv(gdat, k=ifelse(samps[i]<=20, 5, 10))
  smodcv = glasso_cv(sdat, k=ifelse(samps[i]<=20, 5, 10))
  dmodbic = glasso_bic(ddat)
  gmodbic = glasso_bic(gdat)
  smodbic = glasso_bic(sdat)
  dmodbicext = glasso_bicext(ddat)
  gmodbicext = glasso_bicext(gdat)
  smodbicext = glasso_bicext(sdat)
  dthetacv[[i]] = dmodcv[[1]]$wi
  cvmat[1, i] = dmodcv[[1]]$rhomax
  cvl1[1, i] = sum(abs(dmodcv[[1]]$wi))
  cvspar[1, i] = sum(dmodcv[[1]]$wi!=0 & col(dcov)<row(dcov))
  cvspar_round[1, i] = sum(round(dmodcv[[1]]$wi, rounding)!=0 & col(dcov)<row(dcov))
  gthetacv[[i]] = gmodcv[[1]]$wi
  cvmat[2, i] = gmodcv[[1]]$rhomax
  cvl1[2, i] = sum(abs(gmodcv[[1]]$wi))
  cvspar[2, i] = sum(gmodcv[[1]]$wi !=0 & col(gcov)<row(gcov))
  cvspar_round[2, i] = sum(round(gmodcv[[1]]$wi, rounding) !=0 & col(gcov)<row(gcov))
  sthetacv[[i]] = smodcv[[1]]$wi
  cvmat[3, i] = smodcv[[1]]$rhomax
  cvl1[3, i] = sum(abs(smodcv[[1]]$wi))
  cvspar[3, i] = sum(smodcv[[1]]$wi, 1!=0 & col(scov)<row(scov))
  cvspar_round[3, i] = sum(round(smodcv[[1]]$wi, rounding)!=0 & col(scov)<row(scov))
  dthetabic[[i]] = dmodbic[[1]]$wi
  dthetabicext[[i]] = dmodbicext[[1]]$wi
  bicmat[1, i] = dmodbic[[1]]$rhomin
  bicextmat[1, i] = dmodbicext[[1]]$rhomin
  bicl1[1, i] = sum(abs(dmodbic[[1]]$wi))
  bicextl1[1, i] = sum(abs(dmodbicext[[1]]$wi))
  bicspar[1, i] = sum(dmodbic[[1]]$wi!=0 & col(dcov)<row(dcov))
  bicextspar[1, i] = sum(dmodbicext[[1]]$wi!=0 & col(dcov)<row(dcov))
  bicspar_round[1, i] = sum(round(dmodbic[[1]]$wi, rounding) !=0 & col(dcov)<row(dcov))
  bicextspar_round[1, i] = sum(round(dmodbicext[[1]]$wi, rounding) !=0 & col(dcov)<row(dcov))
  gthetabic[[i]] = gmodbic[[1]]$wi
  gthetabicext[[i]] = gmodbicext[[1]]$wi
  bicmat[2, i] = gmodbic[[1]]$rhomin
  bicextmat[2, i] = gmodbicext[[1]]$rhomin
  bicl1[2, i] = sum(abs(gmodbic[[1]]$wi))
  bicextl1[2, i] = sum(abs(gmodbicext[[1]]$wi))
  bicspar[2, i] = sum(gmodbic[[1]]$wi!=0 & col(gcov)<row(gcov))
  bicextspar[2, i] = sum(gmodbicext[[1]]$wi!=0 & col(gcov)<row(gcov))
  bicspar_round[2, i] = sum(round(gmodbic[[1]]$wi, rounding) !=0 & col(gcov)<row(gcov))
  bicextspar_round[2, i] = sum(round(gmodbicext[[1]]$wi,rounding) !=0 & col(gcov)<row(gcov))
  sthetabic[[i]] = smodbic[[1]]$wi
  sthetabicext[[i]] = smodbicext[[1]]$wi
  bicmat[3, i] = smodbic[[1]]$rhomin
  bicextmat[3, i] = smodbicext[[1]]$rhomin
  bicl1[3, i] = sum(abs(smodbic[[1]]$wi))
  bicextl1[3, i] = sum(abs(smodbicext[[1]]$wi))
  bicspar[3, i] = sum(smodbic[[1]]$wi!=0 & col(scov)<row(scov))
  bicextspar[3, i] = sum(smodbicext[[1]]$wi!=0 & col(scov)<row(scov))
  bicspar_round[3, i] = sum(round(smodbic[[1]]$wi, rounding) !=0 & col(scov)<row(scov))
  bicextspar_round[3, i] = sum(round(smodbicext[[1]]$wi, rounding) !=0 & col(scov)<row(scov))
  maxloglik[1, i] = max(dmodcv[[2]])
  maxloglik[2, i] = max(gmodcv[[2]])
  maxloglik[3, i] = max(smodcv[[2]])
  maxloglik[4, i] = log_likelihood(dense, dcov)
  maxloglik[5, i] = log_likelihood(group, gcov)
  maxloglik[6, i] = log_likelihood(sparse, scov)
  minbic[1, i] = min(dmodbic[[2]])
  minbic[2, i] = min(gmodbic[[2]])
  minbic[3, i] = min(smodbic[[2]])
  minbic[4, i] = -samps[i]*log_likelihood(dense, dcov) + log(samps[i])*sum(dense !=0 & col(dcov)<row(dcov))
  minbic[5, i] = -samps[i]*log_likelihood(group, gcov) + log(samps[i])*sum(group !=0 & col(gcov)<row(gcov))
  minbic[6, i] = -samps[i]*log_likelihood(sparse, scov) + log(samps[i])*sum(sparse !=0 & col(scov)<row(scov))
  minbicext[1, i] = min(dmodbicext[[2]])
  minbicext[2, i] = min(gmodbicext[[2]])
  minbicext[3, i] = min(smodbicext[[2]])
  minbicext[4, i] = -samps[i]*log_likelihood(dense, dcov) + log(samps[i])*sum(dense !=0 & col(dcov)<row(dcov)) + 4*sum(dense !=0 & col(dcov)<row(dcov))*nrow(dense)
  minbicext[5, i] = -samps[i]*log_likelihood(group, gcov) + log(samps[i])*sum(group !=0 & col(gcov)<row(gcov)) + 4*sum(group !=0 & col(gcov)<row(gcov))*nrow(group)
  minbicext[6, i] = -samps[i]*log_likelihood(sparse, scov) + log(samps[i])*sum(sparse !=0 & col(scov)<row(scov)) + 4*sum(sparse !=0 & col(scov)<row(scov))*nrow(sparse)
}

# Transform output into dataframes and plot

likelihoods = as.data.frame(t(maxloglik))
likelihoods$n <- samps
colnames(likelihoods) <- c("Dense", "Grouped", "Sparse", "Dense True", "Grouped True", "Sparse True", "n")
likelihoods = melt(likelihoods, id.vars="n")
colnames(likelihoods)[colnames(likelihoods)=="variable"] <- "Model"

bics = as.data.frame(t(minbic))
bics$n <- samps
colnames(bics) <- c("Dense", "Grouped", "Sparse", "Dense True", "Grouped True", "Sparse True", "n")
bics = melt(bics, id.vars="n")



ggplot(likelihoods, aes(x=n, y=value, color=Model)) + geom_line() + ylab("Log Likelihood") + scale_color_manual(breaks = c("Dense", "Grouped", "Sparse", "Dense True", "Grouped True", "Sparse True"),
                                                                                                                values = c("#F8766D", "#00BA38", "#619CFF", "pink","green","blue"))
ggplot(filter(bics, n<500), aes(x=n, y=value, color=variable)) + geom_line() + ylab("BIC")

dsparse_r = sum(dense_round !=0 & col(dense_round)<row(dense_round))
gsparse_r = sum(group_round !=0 & col(group_round)<row(group_round))
ssparse_r = sum(sparse_round !=0 & col(sparse_round)<row(sparse_round))

dsparse = sum(dense !=0 & col(dense)<row(dense))
gsparse = sum(group !=0 & col(group)<row(group))
ssparse = sum(sparse !=0 & col(sparse)<row(sparse))

dl1 = sum(abs(dense))
gl1 = sum(abs(group))
sl1 = sum(abs(sparse))

cvrhos = as.data.frame(t(cvmat))
cvrhos$n = samps
cvrhos$method = "Cross-Validation"
colnames(cvrhos) = c("Dense", "Grouped", "Sparse", "n", "Method")
cvrhos = melt(cvrhos, id.vars=c("Method","n"))
bicrhos = as.data.frame(t(bicmat))
bicrhos$n = samps
bicrhos$method = "BIC"
colnames(bicrhos) = c("Dense", "Grouped", "Sparse", "n", "Method")
bicrhos = melt(bicrhos, id.vars=c("Method","n"))
bicextrhos = as.data.frame(t(bicextmat))
bicextrhos$n = samps
bicextrhos$method = "Extended BIC"
colnames(bicextrhos) = c("Dense", "Grouped", "Sparse", "n", "Method")
bicextrhos = melt(bicextrhos, id.vars=c("Method","n"))
rhosdf = rbind(cvrhos, bicrhos, bicextrhos)
colnames(rhosdf)[colnames(rhosdf)=="value"] <- "rho"
colnames(rhosdf)[colnames(rhosdf)=="variable"] <- "Model"
ggplot(rhosdf, aes(x=n, y=rho, color=Model)) + geom_point() + ylab("rho") + facet_grid(~Method)
ggplot(filter(rhosdf, rho<5), aes(x=n, y=rho, color=Model)) + geom_point() + ylab("lambda") + facet_grid(~Method)


cvspardf = as.data.frame(t(cvspar))
cvspardf$n = samps
cvspardf$method = "Cross-Validation"
colnames(cvspardf) = c("Dense", "Grouped", "Sparse", "n", "Method")
cvspardf = melt(cvspardf, id.vars=c("Method","n"))
bicspardf = as.data.frame(t(bicspar))
bicspardf$n = samps
bicspardf$method = "BIC"
colnames(bicspardf) = c("Dense", "Grouped", "Sparse", "n", "Method")
bicspardf = melt(bicspardf, id.vars=c("Method","n"))
bicextspardf = as.data.frame(t(bicextspar))
bicextspardf$n = samps
bicextspardf$method = "Extended BIC"
colnames(bicextspardf) = c("Dense", "Grouped", "Sparse", "n", "Method")
bicextspardf = melt(bicextspardf, id.vars=c("Method","n"))
spardf = rbind(cvspardf, bicspardf, bicextspardf)

colnames(spardf)[colnames(spardf)=="value"] <- "sparsity"

ggplot(spardf, aes(x=n, y=sparsity, color=variable)) + geom_point() + ylab("Non-sparse Entries") + geom_hline(yintercept=dsparse, linetype="dashed", color = "pink") + geom_hline(yintercept=gsparse, linetype="dashed", color = "green") + geom_hline(yintercept=ssparse, linetype="dashed", color = "blue") + facet_grid(~Method)

cvspardf_round = as.data.frame(t(cvspar_round))
cvspardf_round$n = samps
cvspardf_round$method = "Cross-Validation"
colnames(cvspardf_round) = c("Dense", "Grouped", "Sparse", "n", "Method")
cvspardf_round = melt(cvspardf_round, id.vars=c("Method","n"))
bicspardf_round = as.data.frame(t(bicspar_round))
bicspardf_round$n = samps
bicspardf_round$method = "BIC"
colnames(bicspardf_round) = c("Dense", "Grouped", "Sparse", "n", "Method")
bicspardf_round = melt(bicspardf_round, id.vars=c("Method","n"))
bicextspardf_round = as.data.frame(t(bicextspar_round))
bicextspardf_round$n = samps
bicextspardf_round$method = "Extended BIC"
colnames(bicextspardf_round) = c("Dense", "Grouped", "Sparse", "n", "Method")
bicextspardf_round = melt(bicextspardf_round, id.vars=c("Method","n"))
spardf_round = rbind(cvspardf_round, bicspardf_round, bicextspardf_round)

colnames(spardf_round)[colnames(spardf_round)=="value"] <- "sparsity"
colnames(spardf_round)[colnames(spardf_round)=="variable"] <- "Model"

ggplot(spardf_round, aes(x=n, y=sparsity/190, color=Model)) + scale_y_continuous(labels = scales::percent) + geom_point() + ylab("% Non-Sparse Entries") + geom_hline(yintercept=dsparse_r/190, linetype="dashed", color = "pink") + geom_hline(yintercept=gsparse_r/190, linetype="dashed", color = "green") + geom_hline(yintercept=ssparse_r/190, linetype="dashed", color = "blue") + facet_grid(~Method)
ggplot(spardf_round, aes(x=n, y=sparsity/190, color=Model)) + scale_y_continuous(labels = scales::percent) + geom_point() + ylab("% Non-Sparse Entries") + geom_hline(yintercept=dsparse/190, linetype="dashed", color = "pink") + geom_hline(yintercept=gsparse/190, linetype="dashed", color = "green") + geom_hline(yintercept=ssparse/190, linetype="dashed", color = "blue") + facet_grid(~Method)

cvl1df = as.data.frame(t(cvl1))
cvl1df$n = samps
cvl1df$method = "Cross-Validation"
colnames(cvl1df) = c("Dense", "Grouped", "Sparse", "n", "Method")
cvl1df = melt(cvl1df, id.vars=c("Method","n"))
bicl1df = as.data.frame(t(bicl1))
bicl1df$n = samps
bicl1df$method = "BIC"
colnames(bicl1df) = c("Dense", "Grouped", "Sparse", "n", "Method")
bicl1df = melt(bicl1df, id.vars=c("Method","n"))
bicextl1df = as.data.frame(t(bicextl1))
bicextl1df$n = samps
bicextl1df$method = "Extended BIC"
colnames(bicextl1df) = c("Dense", "Grouped", "Sparse", "n", "Method")
bicextl1df = melt(bicextl1df, id.vars=c("Method","n"))
l1df = rbind(cvl1df, bicl1df, bicextl1df)

colnames(l1df)[colnames(l1df)=="value"] <- "l1norm"
colnames(l1df)[colnames(l1df)=="variable"] <- "Model"

ggplot(l1df, aes(x=n, y=l1norm, color=Model)) + geom_point() + ylab("L1 Norm") + geom_hline(yintercept=dl1, linetype="dashed", color = "pink") + geom_hline(yintercept=gl1, linetype="dashed", color = "green") + geom_hline(yintercept=sl1, linetype="dashed", color = "blue") + facet_grid(~Method)


# Combine rhos and sparsity dfs to see effect of rounding

rhosparsitydf <- rhosdf %>% left_join(spardf)

ggplot(rhosparsitydf, aes(x=rho, y=sparsity, color=variable)) + geom_point() + facet_grid(~Method) + ylab("Non-sparse Entries") + geom_hline(yintercept=dsparse, linetype="dashed", color = "pink") + geom_hline(yintercept=gsparse, linetype="dashed", color = "green") + geom_hline(yintercept=ssparse, linetype="dashed", color = "blue")

ggplot(filter(rhosparsitydf, rho > 0.1), aes(x=rho, y=sparsity, color=variable)) + geom_point() + facet_grid(~Method) + ylab("Non-sparse Entries") + geom_hline(yintercept=dsparse, linetype="dashed", color = "pink") + geom_hline(yintercept=gsparse, linetype="dashed", color = "green") + geom_hline(yintercept=ssparse, linetype="dashed", color = "blue")


rhosparsity_rounddf <- rhosdf %>% left_join(spardf_round)

ggplot(rhosparsity_rounddf, aes(x=rho, y=sparsity, color=variable)) + geom_point() + facet_grid(~Method) + ylab("Non-sparse Entries") + geom_hline(yintercept=dsparse_r, linetype="dashed", color = "pink") + geom_hline(yintercept=gsparse_r, linetype="dashed", color = "green") + geom_hline(yintercept=ssparse_r, linetype="dashed", color = "blue")

comparesparsity = spardf %>% mutate(source = "sparse") %>% left_join(spardf_round, by = c("Method", "n", "variable"))

ggplot(comparesparsity, aes(x=sparsity.x, y=sparsity.y, color=variable)) + facet_grid(~Method) + geom_point()

# Generate graphs for n=15, n-=150, n=1250

s1250 = round(sthetacv[[14]], 1)
g1250 = round(gthetacv[[14]], 1)
d1250 = round(dthetacv[[14]], 1)
g1250_b = round(gthetabic[[14]], 1)

sum(d1250!=0 & col(d1250)<row(d1250))
sum(g1250!=0 & col(g1250)<row(g1250))

sadj1250 = ifelse(s1250!=0 & row(s1250)!=col(s1250),1,0)
gadj1250 = ifelse(g1250!=0 & row(g1250)!=col(g1250),1,0)
dadj1250 = ifelse(d1250!=0 & row(d1250)!=col(d1250),1,0)
gadj1250_b = ifelse(g1250_b!=0 & row(g1250_b)!=col(g1250_b),1,0)

sgraph1250 = graph_from_adjacency_matrix(sadj1250, mode="undirected")
ggraph1250 = graph_from_adjacency_matrix(gadj1250, mode="undirected")
dgraph1250 = graph_from_adjacency_matrix(dadj1250, mode="undirected")
ggraph1250_b = graph_from_adjacency_matrix(gadj1250_b, mode="undirected")

s150 = round(sthetacv[[7]], 1)
g150 = round(gthetacv[[7]], 1)
d150 = round(dthetacv[[7]], 1)

sadj150 = ifelse(s150!=0 & row(s150)!=col(s150),1,0)
gadj150 = ifelse(g150!=0 & row(g150)!=col(g150),1,0)
dadj150 = ifelse(d150!=0 & row(d150)!=col(d150),1,0)

sgraph150 = graph_from_adjacency_matrix(sadj150, mode="undirected")
ggraph150 = graph_from_adjacency_matrix(gadj150, mode="undirected")
dgraph150 = graph_from_adjacency_matrix(dadj150, mode="undirected")

s150_b = round(sthetabic[[7]], 1)
g150_b = round(gthetabic[[7]], 1)
d150_b = round(dthetabic[[7]], 1)

sadj150_b = ifelse(s150_b!=0 & row(s150_b)!=col(s150_b),1,0)
gadj150_b = ifelse(g150_b!=0 & row(g150_b)!=col(g150_b),1,0)
dadj150_b = ifelse(d150_b!=0 & row(d150_b)!=col(d150_b),1,0)

sgraph150_b = graph_from_adjacency_matrix(sadj150_b, mode="undirected")
ggraph150_b = graph_from_adjacency_matrix(gadj150_b, mode="undirected")
dgraph150_b = graph_from_adjacency_matrix(dadj150_b, mode="undirected")

plot(sgraph150_b)
plot(ggraph150_b)
plot(dgraph150_b)

s15 = round(sthetacv[[1]], 1)
g15 = round(gthetacv[[1]], 1)
d15 = round(dthetacv[[1]], 1)

sadj15 = ifelse(s15!=0 & row(s15)!=col(s15),1,0)
gadj15 = ifelse(g15!=0 & row(g15)!=col(g15),1,0)
dadj15 = ifelse(d15!=0 & row(d15)!=col(d15),1,0)

sgraph15 = graph_from_adjacency_matrix(sadj15, mode="undirected")
ggraph15 = graph_from_adjacency_matrix(gadj15, mode="undirected")
dgraph15 = graph_from_adjacency_matrix(dadj15, mode="undirected")

par(mfrow=c(1,3))
sparsegraph_t_r%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)
groupgraph_t_r%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)
densegraph_t_r%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)

par(mfrow=c(1,3))
sgraph15%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)
sgraph150%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)
sgraph1250%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)

par(mfrow=c(1,3))
ggraph15%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)
ggraph150%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)
ggraph1250_b%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)

par(mfrow=c(1,3))
dgraph15%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)
dgraph150%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)
dgraph1250%>% 
  plot(.,vertex.size=15, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)
par(mfrow=c(1,1))


sparsegraph_t_r%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)
sgraph15%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)
sgraph150%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)
sgraph150_b%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)
sgraph1250%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = sparsecoords_rescaled, rescale=F, asp=0)


groupgraph_t_r%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)
ggraph15%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)
ggraph150%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)
ggraph150_b%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)
ggraph1250%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)
ggraph1250_b%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = groupcoords_rescaled, rescale=F, asp=0)


densegraph_t_r%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)
dgraph15%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)
dgraph150%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)
dgraph150_b%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)
dgraph1250%>% 
  plot(.,vertex.size=8, label.color="black", edge.arrow.size=1, layout = densecoords_rescaled, rescale=F, asp=0)
