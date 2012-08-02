# Functions for discrete-time compartment models

dyn.load("inference.so")

hazard.part = function(sys, engine="R") 
{ 
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

hazard = function(sys, engine="R") 
{ 
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


sim.poisson = function(sys, engine="R")
{
    h = hazard(sys)
    if (engine=="R")
    {
        return(rpois(sys$r,h))
    } else if (engine=="C")
    {
        out = .C("sim_poisson",
                 as.integer(sys$s), as.integer(sys$r), as.integer(sys$X),
                 as.integer(t(sys$Pre)), as.double(h), r=integer(sys$r))
        return(out$r)
    } else 
    {
        stop("'engine' must be 'C' or 'R'")
    }
}

update.species = function(sys, nRxns, engine="R") 
{
    if (engine=="R")
    {
        return(as.numeric(sys$X+sys$stoich%*%nRxns))
    } else if (engine=="C")
    {
        out = .C("update_species",
                 as.integer(sys$s), as.integer(sys$r), X=as.integer(sys$X),
                 as.integer(sys$stoich), as.integer(nRxns))
        return(out$X)
    } else 
    {
        stop("'engine' must be 'C' or 'R'")
    }
}


