library(tlpl)

load("PL.RData")

n.sims = length(pl)
plq = list()

for (i in 1:n.sims)
{
  cat("Simulation",i,"(",round(i/n.sims*100),"%)\n")
  plq[[i]] = tlpl_quantile(pl[[i]], verbose=0)
}

save(plq, file="PLq.RData")

q("no")

