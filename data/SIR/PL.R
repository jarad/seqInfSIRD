library(tlpl)
library(plyr)

load("sims.RData")
source("settings.R")

sys$theta = rep(0,sys$r)
prior$X=rmultinom(n.particles, N, sys$X)

pl = llply(lapply(sims, function(x) return(list(y=x$y, tau=1))), 
           tlpl, sckm=sys, prior=prior, n.particles=n.particles, 
           verbose=0, engine="C",
           .progress = progress_text(style=ifelse(interactive(), 3, 1)), 
           .inform=TRUE)

save(pl, file="PL.RData")

q(ifelse(interactive(),"ask","no"))

