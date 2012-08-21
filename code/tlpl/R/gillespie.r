
is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


random.system = function(r=1,s=1)
{
    Pre  = matrix(rpois(s*r, 1),r,s)
    Post = matrix(rpois(s*r, 1),r,s)
    stoich = t(sys$Post-sys$Pre)
    theta = rgamma(r,100,100)
    X = rpois(s,10)+1    

    return(list(r=r,s=s,Pre=Pre,Post=Post,stoich=stoich,theta=theta,X=X))
}


hazard.part = function(sys, engine="R")
{
    engine = pmatch(engine, c("R","C"))

    check.system(sys)

    switch(engine,
    {
        # R implementation
        hp = numeric(sys$r)
        for (i in 1:sys$r) hp[i] = sum(lchoose(sys$X, sys$Pre[i,]))
        return(exp(hp))
    },
    {
        # C implementation
        out = .C("hazard_part_wrap",
                 as.integer(sys$s), as.integer(sys$r), as.integer(t(sys$Pre)), as.integer(sys$X), hp=integer(sys$r))
        return(out$hp)
    })
}


gillespie = function(sys, n, tau) 
{
    # Error checking
    check.system(sys)
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

