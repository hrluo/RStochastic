#Function that simulates a realization of a GP with Gaussian Correlation 
# over a 1:n grid in [0,1]^p.
#-----------------------------------------------------------------------
#Inputs: 
# - scalar n, size of grid
# - scalar p, dimension of domain
# - scalar theta, parameter of Gaussian Correlation (rho = exp(-theta))
# - scalar sigma2, variance of GP (i.e., Cov(y_i,y_i)) 
#Outputs: 
# - n^p array Y, realization of GP over in the n^p dimensional grid
# - vector Yvect, realization of GP in vector form
# - n*p matrix X, grid where the GP is evaluated (each column is one dimension)

simulateGP_GaussianCorr <- function(n, p, theta, sigma2) {
  
  X <- matrix(NA, nrow = n, ncol = p)

  for(i in 1:p) {
    
    x_i <- seq(0, 1, length = n)
    X[,i] <- x_i
    #Generate a uniform grid, when put into use, this can be substituted with real data inputs.
    
    diff_x_i <- abs(outer(x_i, x_i, "-"))
    #distance between any two of these two locations.
    
    R_i <- exp(- theta[i] * diff_x_i^2 )
    #Correlation function for direction i, i=1,...,p
    
    L_i <- t(chol(R_i + diag(n)*.Machine$double.eps*100))
    
    if(i == 1) {
      L <- L_i
    } else {
      L <- L %x% L_i
    }
    
  }
  L <- sqrt(sigma2) * L
  
  Z <- rnorm(n^p, mean = 0, sd = 1)
  #Normal random sample
  
  Yvect <- L %*% Z
  #Realization (in vector form)
  
  Y <- aperm(array(Yvect, dim = rep(n,p))) 
  #Grid with realizations of GP
  
  return(list(Y = Y, 
              Yvect = Yvect,
              X = X))  
}