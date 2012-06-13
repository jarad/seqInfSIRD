require(rjags)

# Contains functions for SIR model
#
# SIRsim  : simulates SIR data
# SIRmcmc : uses JAGS to run mcmc <not implemented yet>

hazard = function(X,theta,N) {
  return(c(theta[1]*X[1]*X[2]/N, # Hazard for S -> I
           theta[2]*X[2]))       # Hazard for I -> R
}

fix.updates <- function(dx,x) {
  dx[1] = min(dx[1],x[1])       # max new infecteds is S
  dx[2] = min(dx[2],x[2]+dx[1]) # max new recovereds is S+(S->I)
  return(dx)
}

update = function(x,dx) {
  c(x[1]-dx[1],x[2]+dx[1]-dx[2])
}

SIRsim = function(X0,theta,p,n) {
  x = y = dx = matrix(NA,n,2)
  N = sum(X0)
  
  # First step
  dx[1,] = rpois(2,hazard(X0,theta,N))
  dx[1,] = fix.updates(dx[1,],X0)
   x[1,] = update(X0,dx[1,])
   y[1,] = rbinom(2,dx[1,],p)
  
  # Remaining steps
  for (i in 2:n) {
  	stopifnot(all(!is.na(x[i-1,])))
  	dx[i,] = rpois(2,hazard(x[i-1,],theta,N))
  	dx[i,] = fix.updates(dx[i,],x[i-1,])
  	 x[i,] = update(x[i-1,],dx[i,])
  	 y[i,] = rbinom(2,dx[i,],p)
  }

  return(list( x = data.frame(S    =  x[,1], I    =  x[,2]),
              dx = data.frame(StoI = dx[,1], ItoR = dx[,2]),
               y = data.frame(StoI =  y[,1], ItoR =  y[,2])))
}




