if (!is.loaded("sim")) dyn.load("../src/tlpl.so")

sim = function(sys, n, tau, while.max=1000) 
{
    check.system(sys)

    out = .C("sim",
             as.integer(sys$s), as.integer(sys$r), as.integer(sys$stoich), as.integer(t(sys$Pre)), 
             as.double(sys$theta), as.double(tau), as.integer(n), as.integer(while.max),
             X=as.integer(rep(sys$X,n+1)))
    return(matrix(out$X, n+1, sys$s, byrow=T))
}

