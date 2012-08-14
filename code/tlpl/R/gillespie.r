if (!is.loaded("sim")) dyn.load("../src/tlpl.so")

is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


gillespie = function(sys, n, tau) 
{
    # Error checking
    check.model(sys)
    stopifnot(tau>0, n>0)

    # if tau is constant
    if (length(tau)==1) tau=rep(tau,n)

    # 
    if (!is.wholenumber(n)) 
    {
        warning("Since n is not a whole number, using round(n) instead.")
        n = round(n)
    }

    out = .C("gillespie_wrap",
             as.integer(sys$s), as.integer(sys$r), as.integer(sys$stoich), as.integer(t(sys$Pre)), 
             as.double(sys$theta), as.double(tau), as.integer(n),
             X=as.integer(rep(sys$X,n+1)))
    return(matrix(out$X, n+1, sys$s, byrow=T))
}

