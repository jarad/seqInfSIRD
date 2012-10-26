# Calculates the quantiles for a mixture of beta distributions
mix.beta.quantiles = function(a,b,probs,w=1)
{
    n = length(a)
    p = length(probs)
    q = rep(NA,p)
    if (w==1) w = rep(1/n,n)

    for (i in 1:length(probs)) 
    {
        tmp = uniroot(function(x) probs[i]-sum(w*pbeta(x,a,b)), c(0,1))
        q[i] = tmp$root
    }

    return(q)
}

# Calculates the quantiles for a mixture of beta distributions
mix.gamma.quantiles = function(a,b,probs,w=1)
{
    n = length(a)
    p = length(probs)
    q = rep(NA,p)
    if (w==1) w = rep(1/n,n)

    for (i in 1:length(probs)) 
    {
        tmp = uniroot(function(x) probs[i]-sum(w*pgamma(x,a,b)), c(0,10))
        q[i] = tmp$root
    }

    return(q)
}


# This function will calculate the desired quantiles for tlpl output.
tlpl_quantile = function(tlpl, probs=c(.025,.5,.975), which="xpr")
{
    which = tolower(which)
    n = dim(tlpl$X)[1]
    s = dim(tlpl$X)[2]
    r = dim(tlpl$hyper$prob$a)[2]
    p = length(probs)

    do =list()
    do$x = grepl("x", which)
    do$p = grepl("p", which)
    do$r = grepl("r", which)
 
    if (do$x) { X.quantiles = array(NA, dim=c(n,s,p)) } else { X.quantiles=NULL }
    if (do$p) { p.quantiles = array(NA, dim=c(n,r,p)) } else { p.quantiles=NULL }
    if (do$r) { r.quantiles = array(NA, dim=c(n,r,p)) } else { r.quantiles=NULL }

    for (i in 1:n) 
    {
        for (j in 1:s) 
        {
            if (do$x) X.quantiles[i,j,] = quantile(tlpl$X[i,j,], probs=probs)
        }
        for (j in 1:r)
        {
            if (do$p) p.quantiles[i,j,] = mix.beta.quantiles( tlpl$hyper$prob$a[i,j,], tlpl$hyper$prob$b[i,j,], probs)
            if (do$r) r.quantiles[i,j,] = mix.gamma.quantiles(tlpl$hyper$rate$a[i,j,], tlpl$hyper$rate$b[i,j,], probs)
        }
    }

    return(list(X.quantiles=X.quantiles, p.quantiles=p.quantiles, r.quantiles=r.quantiles))
}

