x1 <- arima.sim(model = list(order = c(1, 0, 0), ar = .9), n = 128) 
x2 <- arima.sim(model = list(order = c(1, 0, 0), ar = -.5), n = 128) 
x<-as.ts(c(x1,x2))
par(mfrow=c(1,1))
plot(x)
xfit<-arima(x, order = c(20,0,0))
tsdiag(xfit,main='AR(1)')

vsim<-x
library("costat")
vsim.TOS <- BootTOS(as.vector(vsim))
plot(vsim.TOS)
vsim.TOS
# Compute the autocovariance function when treated as a global stationary process.
vsim.acf <- acf(vsim, plot=TRUE)
plot(vsim.acf)
# Compute the localized autocovariance using a wavelet filter of 10 order by default.
# filter.number means how many wavelet basis we are going to use for approximating the local autocovariance.
vsim.lacv <- lacv(vsim, filter.number=6, lag.max=100)
plot(vsim.lacv)

max_lag=10
plot(vsim.lacv, lags=0:max_lag, lcol=1:(max_lag+1),ylim=c(-1,1))
for(i in 0:max_lag) abline(h=0.18^i, lty=2, col=i+1)

abline(v = 50,lty=3,col='black')

pre<-as.data.frame(vsim.lacv$lacv)
seg<-as.vector(pre[,50])
plot(0,0,xlim=c(0,101),ylim=c(0,1),xlab='Lag',ylab='ACF',col='white',main='ACF at time point 50')
for(s in 1:length(seg)){
  segments(s,0,s,seg[s],col='black')
}

plot(vsim.lacv, type="acf", the.time=50)

