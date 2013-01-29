library(Hmisc)

load(paste("LW", round(delta*100), ".RData", sep=""))

source("quantile-settings.R")

lwq = list()

n = dim(lw[[1]]$X)[3]
p = length(probs)
s = dim(lw[[1]]$X)[1]
r = dim(lw[[1]]$p)[1]

Xq = array(NA, dim=c(s,p,n))
pq = rq = array(NA, dim=c(r,p,n))


my.wtd.quantile = function(x,...)
{
  if (all(is.na(x)) | length(x)==0) return(NA)

  wtd.quantile(x,normwt=T,...)
}

for (i in 1:length(lw))
{
  for (j in 1:n)
  {
    Xq[,,j] = t(apply(lw[[i]]$X[,,j], 1, my.wtd.quantile, probs=probs, weights=lw[[i]]$weights[,j]))
    pq[,,j] = t(apply(lw[[i]]$p[,,j], 1, my.wtd.quantile, probs=probs, weights=lw[[i]]$weights[,j]))
    rq[,,j] = t(apply(lw[[i]]$r[,,j], 1, my.wtd.quantile, probs=probs, weights=lw[[i]]$weights[,j]))
  }    

  lwq[[i]] = list(X.quantiles = Xq, 
                  p.quantiles = pq, 
                  r.quantiles = rq)
}

save(lwq, file=paste("LW", round(delta*100), "q.RData", sep=""))

q()

