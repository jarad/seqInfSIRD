# Functions for discrete-time compartment models

dyn.load("inference.so")

hazard.part = function(sys, engine="R") {
    stopifnot(all(sys$Pre%in%c(0,1))) # Dimerizations or above not supported 
 
    if (engine=="R") 
    {
      h = numeric(sys$r)
      for (i in 1:sys$r) h[i] = prod(sys$X^sys$Pre[i,])
      return(h)
    } else if (engine=="C") 
    {
      out = .C("hazard_part", 
               as.integer(sys$s), as.integer(sys$r), as.integer(sys$X),
               as.integer(t(sys$Pre)), h=integer(sys$r))
      return(out$h)
    } else 
    {
        stop("`engine' must be `C' or `R'")
    }
}

hazard = function(sys, engine="R") {
    stopifnot(all(sys$Pre%in%c(0,1))) # Dimerizations or above not supported 
 
    if (engine=="R") 
    {
      h = numeric(sys$r)
      for (i in 1:sys$r) h[i] = sys$theta[i]*prod(sys$X^sys$Pre[i,])
      return(h)
    } else if (engine=="C") 
    {
      out = .C("hazard", 
               as.integer(sys$s), as.integer(sys$r), as.integer(sys$X),
               as.integer(t(sys$Pre)), h=double(sys$r), as.double(sys$theta))
      return(out$h)
    } else 
    {
        stop("`engine' must be `C' or `R'")
    }
}

