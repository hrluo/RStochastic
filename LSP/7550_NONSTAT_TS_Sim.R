set.seed(0)
library("costat")
# This code is modefied from "Cardinali, Alessandro, and Guy P. Nason. "Costationarity of locally stationary time series using costat." Preprint (2012)."

# Simulate a stationary AR(1) model with parameter 0.8.
# Then compute and print the regular ACF computed on the realization
vsim <- arima.sim(model=list(ar=0.8), n=1024)
# Since we are using the Wavelet method to approximate the locally stationary component, n must be 2^k.
# Do the bootstrap test on stationarity on this data set simulated from AR(1)
vsim.TOS <- BootTOS(as.vector(vsim))
plot(vsim.TOS)
vsim.TOS
# Compute the autocovariance function when treated as a global stationary process.
vsim.acf <- acf(vsim, plot=TRUE)
plot(vsim.acf)
# Compute the localized autocovariance using a wavelet filter of 10 order by default.
# filter.number means how many wavelet basis we are going to use for approximating the local autocovariance.
vsim.lacv <- lacv(vsim, filter.number=10, lag.max=100)
plot(vsim.lacv)

# Plot the localized autocovariance curves at different lags. 
# Different lags collectively provides the estimated localized covariance function at different time points.
# dashed lines means the theoretic values of covariance functions from AR(1) for comparison purposes.
max_lag=10
plot(vsim.lacv, lags=0:max_lag, lcol=1:(max_lag+1),ylim=c(-1,1))
for(i in 0:max_lag) abline(h=0.8^i, lty=2, col=i+1)

# Plot the localized acf at time=512 in the style of the classicial acf plot, 
# and superimpose the theoretical values of the acf from the stationary approximation.
ACF_fixedtime<-function(t){
  plot(vsim.lacv,lags=0:max_lag, type="acf", the.time=t)
  for(i in 1:max_lag) segments(i+.5,0,i+.5, 0.8^i, col='red')
  legend("topright", legend=c("Wavelet estimated autocovariance function","Theoretic autocovariance function AR(1), rho=.8"),
         col=c("black", "red"), lty=1, lwd=2,bg="transparent", cex=.75)
}
ACF_fixedtime(512)
# A function that can simulate a time-varying AR process.
tvar1sim <- function (n=512,sd=1) 
{
    arvec <- seq(from = 0.9, to = -0.9, length.out = n)
    for (i in 1:n) arvec[i]<- 2*i+1 
    #Time varying AR coefficient is specified by arvec.
    x <- c(rnorm(1, mean = 0, sd = sd), rep(0, n - 1))
    for (i in 2:n) x[i] <- arvec[i] * x[i - 1] + rnorm(1, mean = 0, 
        sd = sd)
    return(x)
}
tvar1.sim <- tvar1sim()
ts.plot(tvar1.sim)
# Compute and plot the time-localized autocovariance
tvar1.lacv <- lacv(tvar1.sim, filter.number=4, lag.max=50)
max_lag=10
plot(tvar1.lacv, lags=0:max_lag, lcol=1:(max_lag+1),ylim=c(-1,1))
# It can be observed that as the lag increases, the covariance function decreased. 
# This is the typical situation when local components are actually correctly fitted.