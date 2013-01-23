library(tlpl)
library(plyr)

load("sims.RData")
source("filter-settings.R")

sckm$theta = rep(0,sckm$r)
prior$X=sckm$X

pl = llply(lapply(sims, function(x) return(list(y=t(x$y), tau=1))), 
           tlpl, sckm=sckm, prior=prior, n.particles=n.particles, verbose=0,
           .progress = ifelse(interactive(), "text", "none"))
save(pl, file="PL.RData")

q()

