require(rjags)

# Contains functions for SIRD model
#
# SIRDsim  : simulates SIRD data
# SIRDmcmc : uses JAGS to run mcmc <not implemented yet>

hazard = function(X,gamma,N=sum(X)) {
  return(c(gamma[1]*X[1]*X[2]/N, # Hazard for S -> I
           gamma[2]*X[2],        # Hazard for I -> R
           gamma[3]*X[1],        # Hazard for S -> R
           gamma[4]*X[2]))       # Hazard for I -> D
}

fix.updates <- function(dx,x,gamma) {
  p = gamma[2]/sum(gamma[c(2,4)]) # Expected proportion to recover 
  
  dx[1] <- min(c(dx[1],x[1]))
  if ( dx[2]+dx[4] > x[2]+dx[1]) {
    dx[2] <- rbinom(1, x[2]+dx[1], p)
    dx[4] <- x[2]+dx[1]-dx[2]
  }
  dx[3] <- min(c(dx[3],x[1]-dx[1])) 
  return(dx)
}

update <- function(x,dx) {
	c(x[1]-dx[1]-dx[3],x[2]+dx[1]-dx[2],x[3]+dx[2]+dx[3],x[4]+dx[4])
}

SIRDsim = function(X0,gamma,p,n) {
  x = y = dx = matrix(NA,n,4)
  
  # First step
  dx[1,] = rpois(4,hazard(X0,gamma))
  dx[1,] = fix.updates(dx[1,],X0)
  
  x[1,] = update(X0,dx[1,])

  y[1,] = rbinom(4,dx[1,],p)
  
  # Remaining steps
  for (i in 2:n) {
  	stopifnot(all(!is.na(x[i-1,])))
  	dx[i,] = rpois(2,hazard(x[i-1,],gamma))
  	dx[i,] = fix.updates(dx[i,],x[i-1,])
  	
  	x[i,] = update(x[i-1,],dx[i,])
  	
  	y[i,] = rbinom(2,dx[i,],p)
  }

  return(list( x = data.frame(S    =  x[,1], I    =  x[,2], R = x[,3], D = x[,4]),
              dx = data.frame(StoI = dx[,1], ItoR = dx[,2], StoR = dx[,3], ItoD = dx[,4]),
               y = data.frame(StoI =  y[,1], ItoR =  y[,2], StoR =  y[,3], ItoD =  y[,4])))
}




