#
is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


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
                 as.integer(sys$s), as.integer(sys$r), as.integer(t(sys$Pre)), 
                 as.integer(sys$X), hp=integer(sys$r))
        return(out$hp)
    })
}

hazard = function(sys, tau=1, engine="R")
{
    engine = pmatch(engine, c("R","C"))
    
    check.system(sys)

    switch(engine,
    {
        # R implementation
        hp = hazard.part(sys)
        return(list(h=hp*sys$theta,hp=hp))
    },
    {
        # C implementation
        out = .C("hazard_wrap",
                 as.integer(sys$s), as.integer(sys$r), as.integer(t(sys$Pre)),
                 as.double(sys$theta), as.integer(sys$X), as.double(tau),
                 hp=integer(sys$r), h=double(sys$r))
        return(list(h=out$h,hp=out$hp))
    })
}

update.species = function(sys, nr, engine="R")
{
    engine = pmatch(engine, c("R","C"))
    
    check.system(sys)

    switch(engine,
    {
        # R implementation
        return(as.numeric(sys$X + sys$stoich %*% nr))
    },
    {
        # C implementation
        out = .C("update_species_wrap",
                 as.integer(sys$s), as.integer(sys$r), as.integer(sys$stoich),
                 as.integer(nr), X=as.integer(sys$X))
        return(out$X)
    })

}


tau.leap.one.step = function(sys, tau=1, while.max=1000, engine="R")
{
    engine = pmatch(engine, c("R","C"))
    
    check.system(sys)
    stopifnot(tau>0, while.max>0)

    h = hazard(sys,tau,engine="R")$h*tau

    switch(engine,
    {
        # R implementation
        count = 0
        X = rep(-1,sys$s)
        while (any(X<0))
        {
            nr = rpois(sys$r,h)
            X  = update.species(sys,nr,engine="R")
            count = count + 1
            if (count > while.max) 
                stop("R:tau_leap_one_step: Too many unsuccessful simulation iterations.")
        }
        return(list(X=X,nr=nr))
    },
    {
        # C implementation
        out = .C("tau_leap_one_step_wrap",
                 as.integer(sys$s), as.integer(sys$r), as.integer(sys$stoich),
                 as.double(h), as.integer(while.max), 
                 nr=integer(sys$r), X=as.integer(sys$X))
        return(list(X=out$X, nr=out$nr))
    },
    {
        stop(paste("No 'engine' matching: ",engine,".\n")) 
    })

}

tau.leap = function(sys, n=1, tau=1, while.max=1000, engine="R")
{
    engine = pmatch(engine, c("R","C"))
    
    check.system(sys)
    stopifnot(tau>0, while.max>0, n>0)
    if (length(tau)==1) tau=rep(tau,n)
    stopifnot(length(tau)==n)

    switch(engine,
    {
        # R implementation
        X = matrix(sys$X, n+1, sys$s, byrow=T)
        for (i in 1:n) { 
            sys$X = X[i,]
            X[i+1,] = tau.leap.one.step(sys,tau[i],while.max,engine="R")$X
        }
        return(X)
    },
    {
        # C implementation
        out = .C("tau_leap_wrap",
                 as.integer(sys$s), as.integer(sys$r), as.integer(sys$stoich), 
                 as.integer(t(sys$Pre)), as.double(sys$theta), as.double(tau), 
                 as.integer(n), as.integer(while.max),
                 X=as.integer(rep(sys$X,n+1)))
        return(matrix(out$X, n+1, sys$s, byrow=T))    
    },
    {
        stop(paste("No 'engine' matching: ",engine,".\n")) 
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

