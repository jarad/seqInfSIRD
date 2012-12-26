library(tlpl)

load("PL.RData")
source("quantile-settings.R")

n.sims = length(pl)
plq = list()

for (i in 1:n.sims)
{
  cat("Simulation",i,"(",round(i/n.sims*100),"%)\n")
  plq[[i]] = tlpl_quantile(pl[[i]], probs=probs, verbose=0)
}

save(plq, file="PLq.RData")

q("no")

