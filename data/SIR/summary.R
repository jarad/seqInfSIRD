# Read in truth, means, and quantiles and 
# calculate MSE, MAD, MAPE, coverage
#
load("sims.RData")

source("quantile-settings.R")

load("PLq.RData")
load("LW90q.RData"); lw90q = lwq
load("LW95q.RData"); lw95q = lwq
load("LW99q.RData"); lw99q = lwq


# Does PL work?
plot.empty = function() plot(0,0,type="n", frame=F, axes=F, xlab="", ylab="")


plotq = function(q)
{
}


  states = c("S","I","R")
  rxns = c("S->I","I->R")

# PL
  for (j in 1:n.sims) { 
    n.last = which(rowSums(sims[[j]]$nr)==0)[1]
    if (is.na(n.last)) n.last=n
 
    par(mfrow=c(sckm$s,3))
    for (i in 1:sckm$s)
    {
      plot(0,0, type='n', main=states[i], ylim=c(0,N), xlim=c(0,n.last),
           xlab="Time", ylab="Count")
      lines(0:n, sims[[j]]$X[,i], col="red")
      lines(0:n, plq[[j]]$X.quantiles[i,1,]) 
      lines(0:n, plq[[j]]$X.quantiles[i,3,], col="gray") 
      lines(0:n, plq[[j]]$X.quantiles[i,5,]) 
    }

    for (i in 1:sckm$r)
    {
      plot(0,0,, type='n', main=rxns[i], ylim=range(plq[[j]]$r.quantiles[i,,]), xlim=c(0,n.last),
           xlab="Time", ylab="")
      abline(h=sims[[j]]$rates[i], col="red")
      lines(0:n, plq[[j]]$r.quantiles[i,1,]) 
      lines(0:n, plq[[j]]$r.quantiles[i,3,], col="gray") 
      lines(0:n, plq[[j]]$r.quantiles[i,5,]) 
    }
    plot.empty()

    for (i in 1:sckm$r)
    {
      plot(0,0,, type='n', main=rxns[i], ylim=range(plq[[j]]$p.quantiles[i,,]), xlim=c(0,n.last),
           xlab="Time", ylab="")
      abline(h=sims[[j]]$probs[i], col="red")
      lines(0:n, plq[[j]]$p.quantiles[i,1,]) 
      lines(0:n, plq[[j]]$p.quantiles[i,3,], col="gray") 
      lines(0:n, plq[[j]]$p.quantiles[i,5,]) 
    }
    plot.empty()

    readline("")
  }

# LW99
load("LW99q.RData")
  for (j in 1:n.sims) { 
    n.last = which(rowSums(sims[[j]]$nr)==0)[1]
    if (is.na(n.last)) n.last=n
 
    par(mfrow=c(3,sckm$s))
    for (i in 1:sckm$s)
    {
      plot(0,0, type='n', main=states[i], ylim=c(0,N), xlim=c(0,n.last),
           xlab="Time", ylab="Count")
      lines(0:n, sims[[j]]$X[,i], col="red")
      lines(0:n, lwq[[j]]$X.quantiles[i,1,]) 
      lines(0:n, lwq[[j]]$X.quantiles[i,3,], col="gray") 
      lines(0:n, lwq[[j]]$X.quantiles[i,5,]) 
    }

    for (i in 1:sckm$r)
    {
      plot(0,0,, type='n', main=rxns[i], ylim=range(lwq[[j]]$r.quantiles[i,,]), xlim=c(0,n.last),
           xlab="Time", ylab="")
      abline(h=sims[[j]]$rates[i], col="red")
      lines(0:n, lwq[[j]]$r.quantiles[i,1,]) 
      lines(0:n, lwq[[j]]$r.quantiles[i,3,], col="gray") 
      lines(0:n, lwq[[j]]$r.quantiles[i,5,]) 
    }
    plot.empty()

    for (i in 1:sckm$r)
    {
      plot(0,0,, type='n', main=rxns[i], ylim=range(lwq[[j]]$p.quantiles[i,,]), xlim=c(0,n.last),
           xlab="Time", ylab="")
      abline(h=sims[[j]]$probs[i], col="red")
      lines(0:n, lwq[[j]]$p.quantiles[i,1,]) 
      lines(0:n, lwq[[j]]$p.quantiles[i,3,], col="gray") 
      lines(0:n, lwq[[j]]$p.quantiles[i,5,]) 
    }
    plot.empty()
  }


