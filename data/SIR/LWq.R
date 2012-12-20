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

for (i in 1:length(lw))
{
  for (j in 1:n)
  {
    Xq[,,j] = apply(lw[[i]]$X[,,i], 1, wtd.quantile, probs=probs, weights=lw[[i]]$weights[,j])
    pq[,,j] = apply(lw[[i]]$p[,,i], 1, wtd.quantile, probs=probs, weights=lw[[i]]$weights[,j])
    rq[,,j] = apply(lw[[i]]$r[,,i], 1, wtd.quantile, probs=probs, weights=lw[[i]]$weights[,j])
    
    lwq[[i]] = list(X.quantiles = Xq, 
                    p.quantiles = pq, 
                    r.quantiles = rq)
  }
}

save(lwq, file=paste("LWq", round(delta*100), ".RData", sep=""))

#q("no")

