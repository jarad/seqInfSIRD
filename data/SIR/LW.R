# Performs LW for preset tuning parameter

load("sims.RData")
source("filter-settings.R")
source("liu_west.r")

lw = list()

for (i in 1:n.sims)
{
  d = data.matrix(data[data$sim==i & data$time != 0,c("StoI","ItoR")])
  lw[[i]] = liu_west(d, sckm, n.particles=n.particles, delta=delta, prior=prior,
                     method=resampling.function)
}

save(lw, file=paste("LW", round(delta*100), ".RData", sep=""))

q("no")

