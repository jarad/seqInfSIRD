library(plyr)

load("sims.RData")

sir_plot = function(sim) 
{
  par(mfrow=c(1,2))
  plot(0,0, type="n", ylim=range(sim$X), xlim=c(0,n), main="Truth")
  ns = ncol(sim$X)
  for (i in 1:ns) lines(0:n, sim$X[,i], col=i, lty=i, lwd=2)
  legend("topright", c("S","I","R"), col=1:ns, lty=1:ns, lwd=2)

  no = ncol(sim$y)
  plot(0,0, type="n", ylim=range(sim$y), xlim=c(0,n), main="Observations")
  for (i in 1:no) lines(1:n, sim$y[,i], col=i, lty=i, lwd=2)
  legend("topright", c("S->I","I->R"), col=1:no, lty=1:no, lwd=2)

}


pdf("sims.pdf")
l_ply(sims, sir_plot)
dev.off()


q()

