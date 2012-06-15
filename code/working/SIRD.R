require(rjags)

n.states      = 4
n.transitions = 4 # S->I, I->R, S->R, I->D


# Contains functions for SIRD model
#
# SIRDsim  : simulates SIRD data
# SIRDmcmc : uses JAGS to run mcmc <not implemented yet>

hazard = function(X,theta,N) {
  return(c(theta[1]*X[1]*X[2]/N, # Hazard for S -> I
           theta[2]*X[2],        # Hazard for I -> R
           theta[3]*X[1],        # Hazard for S -> R
           theta[4]*X[2]))       # Hazard for I -> D
}

fix.updates <- function(dx,x,theta) {
  p = theta[2]/sum(theta[c(2,4)]) # Expected proportion to recover 
  
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

SIRDsim = function(X0,theta,p,n) {
  x = y = dx = matrix(NA,n,4)
  N = sum(X0)
  
  # First step
  dx[1,] = rpois(n.transitions,hazard(X0,theta,N))
  dx[1,] = fix.updates(dx[1,],X0,theta)
   x[1,] = update(X0,dx[1,])
   y[1,] = rbinom(n.states,dx[1,],p)
  
  # Remaining steps
  for (i in 2:n) {
  	stopifnot(all(!is.na(x[i-1,])))
  	dx[i,] = rpois(n.transitions,hazard(x[i-1,],theta,N))
  	dx[i,] = fix.updates(dx[i,],x[i-1,],theta)
  	 x[i,] = update(x[i-1,],dx[i,])
  	 y[i,] = rbinom(n.states,dx[i,],p)
  }

  # Include initial X in x
  x = rbind(X0,x)

  return(list( x = data.frame(S    =  x[,1], I    =  x[,2], R    =  x[,3],    D =  x[,4]),
              dx = data.frame(StoI = dx[,1], ItoR = dx[,2], StoR = dx[,3], ItoD = dx[,4]),
               y = data.frame(StoI =  y[,1], ItoR =  y[,2], StoR =  y[,3], ItoD =  y[,4])))
}




