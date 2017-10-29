preferPseudoInputs<-function(x,y,M,beta=2){
  p0=exp(beta*x) / max(exp(beta*x))
  points0 = rep(0,length(x))
  points0[ sample(1:length(x), size=M, prob=p0, replace = FALSE) ] = 1
  pp0 = which(points0 == 1)
  return(list(xbar=x[pp0],y.xbar=y[pp0]))
}
set.seed(123)

xbarPrefer<-preferPseudoInputs(x,y,M=4)$xbar
par(mfrow=c(2,2))
plotGpPseudoInputs(M=4,useOpt = T,title="Naive")
plotGpPseudoInputs(M=4,useOpt = F,title="RANDOM")
plotGpPseudoInputs(M=4,useOpt = F,supplied.xbar = xbarPrefer,title="SPS")

xbarPrefer<-preferPseudoInputs(x,y,M=16)$xbar
par(mfrow=c(2,2))
plotGpPseudoInputs(M=16,useOpt = T,title="Naive")
plotGpPseudoInputs(M=16,useOpt = F,title="RANDOM")
plotGpPseudoInputs(M=16,useOpt = F,supplied.xbar = xbarPrefer,title="SPS")
