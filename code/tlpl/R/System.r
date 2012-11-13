# Contains functions for creating stochastic chemical kinetic models or systems


# Creates a test system for use in testing functions
test.system = function(i=1)
{
    foundi = FALSE
    switch(i,
    {            # Constitutive production
        X = c(0)
        Pre  = matrix(0,1,1)
        Post = matrix(1,1,1)
        s = length(X)
        r = nrow(Pre)
        theta =  rgamma(r, 100, 100)
        stoich = t(Post-Pre)
        foundi = TRUE
    }, 
    {
        X = c(0,1,2)
        Pre = rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1),c(1,1,0),c(1,0,1),c(0,1,1),c(1,1,1))
        s = length(X)
        r = nrow(Pre)
        theta =  rgamma(r, 100, 100)
        Post = matrix(1,nrow=r,ncol=s)
        stoich = t(Post-Pre)
        foundi = TRUE
    })

    # Used in place of a default switch case
    if (!foundi) stop("Test system doesn't exist.")

    return(list(r=r,s=s,Pre=Pre,Post=Post,stoich=stoich,theta=theta,X=X))
}


null.system = function(r=2,s=3)
{
    Pre  = matrix(0,r,s)
    Post = matrix(0,r,s)
    stoich = t(Post-Pre)
    theta = rep(0,r)
    X = rep(0,s)    

    return(list(r=r,s=s,Pre=Pre,Post=Post,stoich=stoich,theta=theta,X=X))
}

# 
random.system = function(seed=1,r=NULL,s=NULL)
{
    set.seed(seed)
    if (is.null(r)) r = rpois(1,2)+1
    if (is.null(s)) s = rpois(1,3)+1
    Pre  = matrix(rpois(r*s, 1),r,s)
    Post = matrix(rpois(r*s, 1),r,s)
    stoich = t(Post-Pre)
    theta = rgamma(r,100,100)
    X = rpois(s,10)+1    

    return(list(r=r,s=s,Pre=Pre,Post=Post,stoich=stoich,theta=theta,X=X))
}



check.system = function(sys) {
    stopifnot(sys$s == ncol(  sys$Pre),
              sys$s == ncol(  sys$Post),
              sys$s == nrow(  sys$stoich),
              sys$s == length(sys$X),
              sys$r == nrow(  sys$Pre),
              sys$r == nrow(  sys$Post),
              sys$r == ncol(  sys$stoich),
              sys$r == length(sys$theta),
              sys$r == length(sys$lmult),
              all(sys$theta>=0),
              all(sys$X    >=0),
              all(sys$Pre  >=0),
              all(sys$Post >=0),
              all(sys$mult >=0))

    if (!all.equal(sys$stoich,t(sys$Post-sys$Pre))) warning("sys$stoich!=t(sys$Post-sys$Pre))")
}


