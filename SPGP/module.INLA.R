library(INLA)
### simulate 2 realizations of the latent field
S0 <- GaussRF(x=x, y=y, model=model, grid=TRUE,
              param=c(0, variance, nugget, scale))
S0 = S0 - mean(S0) + mean0
S0.vec =  inla.matrix2vector(S0)

S1 <- GaussRF(x=x, y=y, model=model, grid=TRUE,
              param=c(0, variance, nugget, scale))
S1 = S1 - mean(S1) + mean1
S1.vec =  inla.matrix2vector(S1)

###map to probability
p0=exp(beta*S0.vec) / max(exp(beta*S0.vec))
p1=exp(beta*S1.vec) / max(exp(beta*S1.vec))

### Simulate the data on S0 using a preferential sampling design

### sample points
points0 = rep(0,n)
points0[ sample(1:n, size=n.data, prob=p0, replace = FALSE) ] = 1

### observe the marks
y0 = rep(NA,n)
pp0 = which(points0 == 1)
y0[pp0] = rnorm(n.data, mean = S0.vec[pp0], sd=sd.data)

loc0=rep(NA,n)
loc0[pp0]=1
#This piece of code is directly from
#http://www.r-inla.org/examples/case-studies/diggle09/simulation2