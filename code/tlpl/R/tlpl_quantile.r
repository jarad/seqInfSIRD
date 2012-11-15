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
    n = dim(tlpl$X)[3]
    s = dim(tlpl$X)[1]
    r = dim(tlpl$hyper$prob$a)[1]
    p = length(probs)

    do =list()
    do$x = grepl("x", which)
    do$p = grepl("p", which)
    do$r = grepl("r", which)
 
    if (do$x) { X.quantiles = array(NA, dim=c(s,p,n)) } else { X.quantiles=NULL }
    if (do$p) { p.quantiles = array(NA, dim=c(r,p,n)) } else { p.quantiles=NULL }
    if (do$r) { r.quantiles = array(NA, dim=c(r,p,n)) } else { r.quantiles=NULL }

    for (i in 1:n) 
    {
        for (j in 1:s) 
        {
            if (do$x) X.quantiles[j,,i] = quantile(tlpl$X[j,,i], probs=probs)
        }
        for (j in 1:r)
        {
            if (do$p) p.quantiles[j,,i] = mix.beta.quantiles( tlpl$hyper$prob$a[j,,i], tlpl$hyper$prob$b[j,,i], probs)
            if (do$r) r.quantiles[j,,i] = mix.gamma.quantiles(tlpl$hyper$rate$a[j,,i], tlpl$hyper$rate$b[j,,i], probs)
        }
    }

    return(list(X.quantiles=X.quantiles, p.quantiles=p.quantiles, r.quantiles=r.quantiles))
}

