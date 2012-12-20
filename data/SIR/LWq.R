library(Hmisc)

load(paste("LW", round(delta*100), ".RData", sep=""))

source("quantile-settings.R")

lwq = list()

for (i in 1:length(lw))
{
  lwq[[i]] = list(X.quantiles = apply(lw[[i]]$X, c(1,3), quantile, probs=probs, na.rm=T), 
                  p.quantiles = apply(lw[[i]]$p, c(1,3), quantile, probs=probs, na.rm=T), 
                  r.quantiles = apply(lw[[i]]$r, c(1,3), quantile, probs=probs, na.rm=T))
}

save(lwq, file=paste("LWq", round(delta*100), ".RData", sep=""))

q("no")

