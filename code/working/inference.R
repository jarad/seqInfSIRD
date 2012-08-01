N.RXNS <<- 4
# Reactions
# [1] S -> I
# [2] I -> R
# [3] S -> R, vaccinations
# [4] I -> D

hazard.R <- function(X, N=sum(X)) {
  hazard    <- rep(NA, N.RXNS)
  hazard[1] <- X[1]*X[2]/N
  hazard[2] <- X[2]
  hazard[3] <- X[1]
  hazard[4] <- X[2]
  return(hazard)
}

simulate.one.step.R <- function(X, theta, N) {
  #if (X[2]==0) return(list(X=X, dX=rep(0, N.RXNS)))

  h  <- hazard.R(X, sum(X))
  dX <- rpois(N.RXNS, theta*h)

  # Restrict to maintain non-negative states
  dX[1] <- min(dX[1], X[1])
  if (dX[2]+dX[4] > X[2]+dX[1]) {
    dX[2] <- rbinom(1, X[2]+dX[1], h[2]/(h[2]+h[4]))
    dX[4] <- X[2]+dX[1]-dX[2]

  }
  dX[3] <- min(dX[3], X[1]-dX[1])            # If you become infected and vaccinated, infection occurs

  # Update states
  X[1] <- X[1]-dX[1]      -dX[3]
  X[2] <- X[2]+dX[1]-dX[2]      -dX[4]
  X[3] <- X[3]      +dX[2]+dX[3]
  X[4] <- X[4]                  +dX[4]
  stopifnot(all(X>=0))
  return(list(X=X, dX=dX))
}

simulate.R <- function(X, theta, N=colSums(X), n.steps=1) {
  n.reps    <- ncol(X)
  n.states  <- nrow(X)
  Xout      <- dX <- array(NA, dim=c(n.states, n.steps+1, n.reps))
  Xout[,1,] <- X

  for (j in 1:n.reps) {
    for (i in 1:n.steps) {
      tmp          <- simulate.one.step.R(Xout[,i,j], theta, N[j])
      Xout[,i+1,j] <- tmp$X
      dX[  ,i+1,j] <- tmp$dX
    }
  }

  return(list(X=Xout,dX=dX))
}

inference.one.step.R <- function(X0, dX, N, prop, hyper, sample=FALSE) {

  h  <- hazard.R(X0, N)
  for (i in 1:N.RXNS) {
  	# Should these have the prop[i] removed according to equations 
  	# 2e and 2f of fluControl/new/measlesPL.pdf?
    hyper[i,1] <- hyper[i,1] + ifelse(sample, rbinom(1,dX[i], prop[i]), prop[i]*dX[i])
    hyper[i,2] <- hyper[i,2] + prop[i]* h[i]
  }
  return(hyper)
}


one.step.R <- function(X, hyper, theta, prop, sample=FALSE) {
  # X    :   n.states x n.reps matrix of starting states
  # hyper:   N.RXNS x 2 x n.reps array of starting prior hyperparameters
  # theta:   N.RXNS x n.reps matrix of true theta values
  # prop :   N.RXNS x n.reps vector of sampling proportions

  N      = colSums(X)
  n.reps = ncol(X)
  dX     = matrix(NA,N.RXNS,n.reps) 

  for (i in 1:n.reps) {
    tmp        = simulate.one.step.R(X[,i], theta[,i], N[i])
    dX[,i]     = tmp$dX
    hyper[,,i] = inference.one.step.R(X[,i], dX[,i], N[i], prop[,i], hyper[,,i], sample)
    X[,i]      = tmp$X # Must be after inference step
  }

  return(list(X=X, hyper=hyper, newX=dX))
  # X    : n.states x n.reps matrix of ending states
  # newX : N.RXNS x n.reps matrix of transitions
  # hyper: N.RXNS x 2 x.nreps array of ending prior hyperparameters
}



inference.R <- function(X, dX, N, prop, hyper=NULL, sample=FALSE) {
  # X    :   n.states x n.steps x n.reps array of starting states
  # dX   :   N.RXNS x n.steps x n.reps array of transitions
  # hyper:   N.RXNS x 2 x n.reps array of starting prior hyperparameters
  # prop :   N.RXNS x n.reps vector of sampling proportions (constant for each rep)

  n.states <- dim(X)[1]
  n.steps  <- dim(X)[2]
  n.reps   <- dim(X)[3]

  if (is.null(hyper)) {
    hyper <- array(1, dim=c(N.RXNS,2,n.reps))
  }
  tmp <- hyper
  hyper <- array(NA, dim=c(N.RXNS,2,n.steps,n.reps))
  hyper[,,1,] <- tmp

  for (j in 1:n.reps) {
    for (s in 2:n.steps) {
      hyper[,,s,j] <- inference.one.step.R(X[,s-1,j], dX[,s,j], N[j], prop[,j], hyper[,,s-1,j], sample)
    }
  }

  return(hyper)
  # hyper:   N.RXNS x 2 x n.steps x n.reps array of starting prior hyperparameters
}

#################################################################################
# C wrappers
#################################################################################
dyn.load(paste("inference-",.Platform$r_arch,.Platform$dynlib.ext,sep=''))

hazard.C <- function(X, N=sum(X)) {
  stopifnot(length(X)>1)
  out <- .C("hazard", as.integer(X), as.integer(N), hazard=double(N.RXNS))
  return(out$hazard)
}

simulate.one.step.C <- function(X, theta, N) {
  n.states <- length(X)
  stopifnot(length(theta)==N.RXNS, n.states==4)
  out <- .C("simulate_one_step", X=as.integer(X), as.double(theta), as.integer(N), as.integer(N.RXNS),
            dX=integer(n.states))
  return(list(X=out$X, dX=out$dX))
}

simulate.C <- function(X, theta, N=colSums(X), n.steps=1) {
  n.reps    <- ncol(X)
  n.states  <- nrow(X)
  dims      <- c(n.states, n.steps+1, n.reps)
  Xout      <- dX <- array(NA, dim=dims)
  Xout[,1,] <- X

  stopifnot(length(N)==n.reps)
  out <- .C("simulate", X=as.integer(Xout), as.double(theta), as.integer(N), as.integer(n.steps),
            as.integer(n.states), as.integer(n.reps), as.integer(N.RXNS), dX=as.integer(dX),
            NAOK=TRUE)

  return(list(X=array(out$X, dim=dims),dX=array(out$dX, dim=dims)))
}

inference.one.step.C <- function(X0, dX, N, prop, hyper, sample=FALSE) {
  n.states <- length(X0)
  stopifnot(n.states==4, length(dX)==N.RXNS, length(prop)==N.RXNS, length(hyper)==2*N.RXNS)
  out <- .C("inference_one_step", as.integer(X0), as.integer(dX), as.integer(N),
            as.double(prop), hyper=as.double(hyper), as.integer(N.RXNS), as.integer(sample))

  return(matrix(out$hyper, N.RXNS, 2))
}

inference.C <- function(X, dX, N, prop, hyper=NULL, sample=FALSE) {
  # prop represents the sampling proportion for each transition: 1 is perfect sampling, 0 is none
  n.reps   <- ncol(X)
  n.states <- nrow(X)

  if (is.null(hyper)) {
    hyper <- array(1, dim=c(N.RXNS,2,n.reps))
  }

  stopifnot(dim(dX)==c(N.RXNS, n.reps), length(N)==n.reps, dim(prop)==c(N.RXNS, n.reps),
            dim(hyper)==c(N.RXNS,2,n.reps))
  #out <- .C("inference",

  return(hyper)
}

one.step.C <- function(X, hyper, theta, prop, sample=FALSE) {
  # X     : n.states x n.reps matrix of starting states
  # hyper : N.RXNS x 2 x n.reps array of starting prior hyperparameters
  # theta : N.RXNS x n.reps matrix of true theta values
  # prop  : N.RXNS x n.reps vector of sampling proportions
  # sample: whether to sample the number of reactions

  N        <- colSums(X)
  n.reps   <- ncol(X)
  n.states <- nrow(X)
  hyper.old <- hyper
  dX        = matrix(0,N.RXNS,n.reps) # No need to initialize as a matrix

  stopifnot(dim(hyper)==c(N.RXNS,2,n.reps), dim(theta)==c(N.RXNS, n.reps), dim(prop)==c(N.RXNS, n.reps))
  out <- .C("one_step", X=as.integer(X), hyper=hyper, as.double(theta),
            as.double(prop), as.integer(N), as.integer(n.reps), as.integer(n.states), 
            as.integer(N.RXNS), as.integer(sample),
            dX=as.integer(dX))

  hyper = array(out$hyper, dim=c(N.RXNS,2,n.reps))  
  dX    = matrix(out$dX, N.RXNS,n.reps)

  return(list(X=matrix(out$X, n.states, n.reps), hyper=hyper, newX=dX))
  # X    : n.states x n.reps matrix of ending states
  # hyper: N.RXNS x 2 x n.reps array of ending prior hyperparameters
  # newX : N.RXNS x n.reps matrix of transitions
}
