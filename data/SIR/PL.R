library(tlpl)
library(plyr)

load("sims.RData")
source("settings.R")

sckm$theta = rep(0,sckm$r)
prior$X=rmultinom(n.particles, N, sckm$X)

pl = llply(lapply(sims, function(x) return(list(y=t(x$y), tau=1))), 
           tlpl, sckm=sckm, prior=prior, n.particles=n.particles, verbose=0,
           .progress = progress_text(style=ifelse(interactive(), 3, 1)), .inform=TRUE)
save(pl, file="PL.RData")

q(ifelse(interactive(),"ask","no"))

