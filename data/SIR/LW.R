# Performs LW for preset tuning parameter
library(plyr)

load("sims.RData")
source("settings.R")
source("liu_west.r")

sckm$lmult = log(sckm$mult) # Temporary

prior$X = rmultinom(n.particles, N, sckm$X)

lw = llply(lapply(sims, function(x) return(x$y)), liu_west,
           sckm=sckm, 
           n.particles=n.particles,
           delta=delta, 
           prior=prior,
           method=resampling.function,
           .progress = progress_text(style=ifelse(interactive(), 3, 1)))

save(lw, file=paste("LW", round(delta*100), ".RData", sep=""))

q(ifelse(interactive(),"ask","no"))

