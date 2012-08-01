# Functions for discrete-time compartment models

dyn.load("inference.so")

hazard.part = function(sys, engine="R") {
    stopifnot(all(sys$Pre%in%c(0,1,2))) # Dimerizations or above not supported 
 
    if (engine=="R") 
    {
      h = rep(1,sys$r)
      for (i in 1:sys$r) 
      {
        for (j in 1:sys$s) 
        {
          switch(sys$Pre[i,j]+1,
                 {}, 
                 { h[i] = h[i]*sys$X[j]},
                 { h[i] = h[i]*sys$X[j]*(sys$X[j]-1)/2},
                 { stop("Pre must have elements 0, 1, and 2.") })
        }
      }          
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

