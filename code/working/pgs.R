# This function will run a naive univariate particle Gibbs sampler.
# 
# Input:
#   y  : data
#   od : log observation density (mass)
#   es : evolution sampler
#   ps : parameter sampler
#   J  : number of particles
#   pr : prior sampler for x
#   iv : initial values
#   n.iter : number of MCMC iterations
#
# Output:
#   theta : an n.iter x p matrix of fixed parameter samples
#   x     : an n.iter x length(y) of latent states samples
library(smcUtils)

pgs <- function(y, od, er, ps, pr, iv, J=length(y)*5, n.iter=10) {
  n = length(y)
  x = iv$x; theta = iv$theta # Initial values
  w = rep(1,J)               # Particle weights
  xp = A = matrix(NA,J,n)    # For scalar x
  A[1,] = 1                  # Ancestor, current particle is included
  
  # Store samples
  theta.reps = matrix(NA,n.iter,length(theta)) 
  x.reps     = matrix(NA,n.iter,n)
  
  # MCMC
  cat("  |",rep(" ",20),"|100%\n  |",sep="")
  for (i in 1:n.iter) {
  	if (i%%(n.iter/20)==0) cat("*")
  	# Draw x	
    xp[1,]   = x             # Current trajectory is included
    xp[-1,1] = pr(J-1,theta) # Particle draws for x0
    w        = renormalize.weights(od(y[1], xp[,1]), log=T)
  	for (j in 2:n) {
  	  A[-1,j]  = resample(w,J-1,"stratified")$indices
  	  xp[-1,j] = er(xp[A[-1,j],j-1],theta)
  	  w        = renormalize.weights(od(y[j],xp[A[,j],j]), log=T)  
  	}
  	
  	# Sample path
  	k = sample(J,1,prob=w)
  	for (j in n:1) {
  	  x[j] = xp[k,j]
  	  k     = A[k,j]
  	}
  	
  	
  	theta = ps(x)
  	
  	theta.reps[i,] = theta
  	x.reps[i,] = x
  }
  cat("|\n")
  
  return(list(theta=theta.reps, x=x.reps))
}
