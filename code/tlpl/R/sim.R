sim = function(sys, n, tau, while.max=1000) 
{
    out = .C("sim",
             as.integer(sys$s), as.integer(sys$r), as.integer(sys$stoich), as.integer(t(sys$Pre)), 
             as.double(sys$theta), as.double(tau), as.integer(n), as.integer(while.max),
             X=as.integer(rep(sys$X,n)))
    return(matrix(out$X, n, sys$s, byrow=T))
}

