if (!is.loaded("sim")) dyn.load("../src/tlpl.so")

is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


sim = function(sys, n, tau, while.max=1000) 
{
    # Error checking
    check.system(sys)
    stopifnot(tau>0, n>0, while.max>0)

    # if tau is constant
    if (length(tau)==1) tau=rep(tau,n)

    # 
    if (!is.wholenumber(n)) 
    {
        warning("Since n is not a whole number, using round(n) instead.")
        n = round(n)
    }
    if (!is.wholenumber(while.max)) 
    {
        warning("Since while.max is not an integer, using round(while.max) instead.")
        while.max=round(while.max)
    }


    out = .C("sim",
             as.integer(sys$s), as.integer(sys$r), as.integer(sys$stoich), as.integer(t(sys$Pre)), 
             as.double(sys$theta), as.double(tau), as.integer(n), as.integer(while.max),
             X=as.integer(rep(sys$X,n+1)))
    return(matrix(out$X, n+1, sys$s, byrow=T))
}

