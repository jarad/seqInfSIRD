# Functions for discrete-time compartment models

hazard.part = function(sys, engine="R") {
    if (engine!="R") stop("Only R is implemented")
    stopifnot(all(sys$Pre%in%c(0,1))) # Dimerizations or above not supported 
  
    hp = numeric(sys$r)
    for (i in 1:sys$r) hp[i] = prod(sys$X^sys$Pre[i,])
    return(hp)
}

hazard = function(sys, theta,...) {
    return(theta*hazard.part(sys,...))
}
