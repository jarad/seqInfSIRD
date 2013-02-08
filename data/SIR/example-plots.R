# Creates plots to see how the methods compare to the data on a single run

load("sims.RData")
n.sims = length(sims)

source("quantile-settings.r")
ll = which(probs==.025); ll=1 # current the former doesn't match
ul = which(probs==.975)


load("PLq.RData")
load("LW90q.RData"); lw90q = lwq
load("LW95q.RData"); lw95q = lwq
load("LW99q.RData"); lw99q = lwq

methods=c("Truth","PL","LW90","LW95","LW99")

pdf("example-plots.pdf")
for (i in 1:n.sims)
{
  sim = sims[[i]]
  plq0 = plq[[i]]
  lw90q0 = lw90q[[i]]
  lw95q0 = lw95q[[i]]
  lw99q0 = lw99q[[i]]

  par(mfrow=c(3,3))

  # States 
  for (j in 1:sckm$s)
  {
    plot(0:n,sim$X[,j], type="l", ylim=c(0,N), xlim=c(0,n), main=sckm$states[j])
    lines(0:n, plq0  $X.quantiles[j,ll,], col=2)
    lines(0:n, plq0  $X.quantiles[j,ul,], col=2)

    lines(0:n, lw90q0$X.quantiles[j,ll,], col=3)
    lines(0:n, lw90q0$X.quantiles[j,ul,], col=3)

    lines(0:n, lw95q0$X.quantiles[j,ll,], col=4)
    lines(0:n, lw95q0$X.quantiles[j,ul,], col=4)

    lines(0:n, lw99q0$X.quantiles[j,ll,], col=5)
    lines(0:n, lw99q0$X.quantiles[j,ul,], col=5)

    legend("left", methods, col=1:5, lty=1)
  }

  # Probabilities
  for (j in 1:sckm$r)
  {
    plot(0,0, type="n", ylim=range(plq0$p.quantiles[j,,]), xlim=c(0,n), main=paste("p:",sckm$rxns[j]))
    abline(h=sim$probs[j])
    lines(0:n, plq0  $p.quantiles[j,ll,], col=2)
    lines(0:n, plq0  $p.quantiles[j,ul,], col=2)

    lines(0:n, lw90q0$p.quantiles[j,ll,], col=3)
    lines(0:n, lw90q0$p.quantiles[j,ul,], col=3)

    lines(0:n, lw95q0$p.quantiles[j,ll,], col=4)
    lines(0:n, lw95q0$p.quantiles[j,ul,], col=4)

    lines(0:n, lw99q0$p.quantiles[j,ll,], col=5)
    lines(0:n, lw99q0$p.quantiles[j,ul,], col=5)

    legend("left", methods, col=1:5, lty=1)
  }

  plot(0,0,type="n", axes=F, xlab="", ylab="")

  # Rates
  for (j in 1:sckm$r)
  {
    plot(0,0, type="n", ylim=range(plq0$r.quantiles[j,,]), xlim=c(0,n), main=paste("r:",sckm$rxns[j]))
    abline(h=sim$rates[j])
    lines(0:n, plq0  $r.quantiles[j,ll,], col=2)
    lines(0:n, plq0  $r.quantiles[j,ul,], col=2)

    lines(0:n, lw90q0$r.quantiles[j,ll,], col=3)
    lines(0:n, lw90q0$r.quantiles[j,ul,], col=3)

    lines(0:n, lw95q0$r.quantiles[j,ll,], col=4)
    lines(0:n, lw95q0$r.quantiles[j,ul,], col=4)

    lines(0:n, lw99q0$r.quantiles[j,ll,], col=5)
    lines(0:n, lw99q0$r.quantiles[j,ul,], col=5)

    legend("left", methods, col=1:5, lty=1)
  }

}
dev.off()





