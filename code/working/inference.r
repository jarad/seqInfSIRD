# Functions for discrete-time compartment models

dyn.load("inference.so")

hazard.part = function(sys, engine="R") { 
    if (engine=="R") 
    {
      h = rep(1,sys$r)
      for (i in 1:sys$r) 
      {
        for (j in 1:sys$s) 
        {
          h[i] = h[i]*choose(sys$X[j], sys$Pre[i,j]) 
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
    if (engine=="R") 
    {
      h = hazard.part(sys)
      for (i in 1:sys$r) h[i] = h[i]*sys$theta[i]
      return(h)
    } else if (engine=="C") 
    {
      out = .C("hazard", 
               as.integer(sys$s), as.integer(sys$r), as.integer(sys$X),
               as.integer(t(sys$Pre)), as.double(sys$theta), h=double(sys$r))
      return(out$h)
    } else 
    {
        stop("`engine' must be `C' or `R'")
    }
}

