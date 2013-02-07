# Read in truth, means, and quantiles and 
# calculate MSE, MAD, MAPE, coverage
#
load("sims.RData")

source("quantile-settings.R")
N.METHODS = 4
load("PLq.RData")
load("LW90q.RData"); lw90q = lwq
load("LW95q.RData"); lw95q = lwq
load("LW99q.RData"); lw99q = lwq
method = c("PL","LW90","LW95","LW99")

# 
if (is.null(sckm$states)) sckm$states = c("S","I","R")
if (is.null(sckm$rxns))   sckm$rxns   = c("S->I","I->R")

states = sckm$states
rxns   = sckm$rxns
ns     = sckm$s
nr     = sckm$r
comp = c(states, paste("p:",rxns), paste("r:",rxns))


nn = c(length(plq), ns+2*nr, n+1, N.METHODS)

## Convert to data frames for saving
library(reshape2)
array.to.df = function(arr) 
{
  df = melt(arr, varnames=c("parameter","time","method"))
  df$parameter = as.factor(comp[df$parameter])
  df$time      = df$time-1
  df$method    = method[df$method]
  return(df)
}


#########################################################



# MSE
which.md = which(probs==0.5)
SE = array(NA, dim=nn)
for (i in 1:nn[1])
{
  for (j in 1:sckm$s)
  {
    SE[i,j,,1] = (sims[[i]]$X[,j] - plq  [[i]]$X.quantiles[j,which.md,])^2
    SE[i,j,,2] = (sims[[i]]$X[,j] - lw90q[[i]]$X.quantiles[j,which.md,])^2
    SE[i,j,,3] = (sims[[i]]$X[,j] - lw95q[[i]]$X.quantiles[j,which.md,])^2
    SE[i,j,,4] = (sims[[i]]$X[,j] - lw99q[[i]]$X.quantiles[j,which.md,])^2
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s + j
    SE[i,k,,1] = (sims[[i]]$probs[j] - plq  [[i]]$p.quantiles[j,which.md,])^2
    SE[i,k,,2] = (sims[[i]]$probs[j] - lw90q[[i]]$p.quantiles[j,which.md,])^2
    SE[i,k,,3] = (sims[[i]]$probs[j] - lw95q[[i]]$p.quantiles[j,which.md,])^2
    SE[i,k,,4] = (sims[[i]]$probs[j] - lw99q[[i]]$p.quantiles[j,which.md,])^2
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s + sckm$r + j
    SE[i,k,,1] = (sims[[i]]$rates[j] - plq  [[i]]$r.quantiles[j,which.md,])^2
    SE[i,k,,2] = (sims[[i]]$rates[j] - lw90q[[i]]$r.quantiles[j,which.md,])^2
    SE[i,k,,3] = (sims[[i]]$rates[j] - lw95q[[i]]$r.quantiles[j,which.md,])^2
    SE[i,k,,4] = (sims[[i]]$rates[j] - lw99q[[i]]$r.quantiles[j,which.md,])^2
  }
}

MSE = apply(SE, 2:4, mean, na.rm=T)


# MAD
which.md = which(probs==0.5)
AD = array(NA, dim=nn)
for (i in 1:nn[1])
{
  for (j in 1:sckm$s)
  {
    AD[i,j,,1] = abs(sims[[i]]$X[,j] - plq  [[i]]$X.quantiles[j,which.md,])
    AD[i,j,,2] = abs(sims[[i]]$X[,j] - lw90q[[i]]$X.quantiles[j,which.md,])
    AD[i,j,,3] = abs(sims[[i]]$X[,j] - lw95q[[i]]$X.quantiles[j,which.md,])
    AD[i,j,,4] = abs(sims[[i]]$X[,j] - lw99q[[i]]$X.quantiles[j,which.md,])
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s + j
    AD[i,k,,1] = abs(sims[[i]]$probs[j] - plq  [[i]]$p.quantiles[j,which.md,])
    AD[i,k,,2] = abs(sims[[i]]$probs[j] - lw90q[[i]]$p.quantiles[j,which.md,])
    AD[i,k,,3] = abs(sims[[i]]$probs[j] - lw95q[[i]]$p.quantiles[j,which.md,])
    AD[i,k,,4] = abs(sims[[i]]$probs[j] - lw99q[[i]]$p.quantiles[j,which.md,])
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s + sckm$r + j
    AD[i,k,,1] = abs(sims[[i]]$rates[j] - plq  [[i]]$r.quantiles[j,which.md,])
    AD[i,k,,2] = abs(sims[[i]]$rates[j] - lw90q[[i]]$r.quantiles[j,which.md,])
    AD[i,k,,3] = abs(sims[[i]]$rates[j] - lw95q[[i]]$r.quantiles[j,which.md,])
    AD[i,k,,4] = abs(sims[[i]]$rates[j] - lw99q[[i]]$r.quantiles[j,which.md,])
  }
}

MAD = apply(AD, 2:4, mean, na.rm=T)



# MAPE
which.md = which(probs==0.5)
APE = array(NA, dim=nn)
for (i in 1:nn[1])
{
  for (j in 1:sckm$s)
  {
    APE[i,j,,] = AD[i,j,,]/sims[[i]]$X[,j]
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s + j
    APE[i,k,,] = AD[i,k,,]/sims[[i]]$probs[j]    
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s + sckm$r + j
    APE[i,k,,] = AD[i,k,,]/sims[[i]]$rates[j]
  }
}

MAPE = apply(APE, 2:4, mean, na.rm=T)





# Coverage
which.ci = c(1,5)
CV = array(NA, dim=nn)

vIsTRUE = Vectorize(isTRUE)
inside = function(x,a,b) return(ifelse(vIsTRUE(x>=a & x<=b), TRUE, FALSE))
for (i in 1:nn[1])
{
  for (j in 1:sckm$s)
  {
    CV[i,j,,1] = inside(sims[[i]]$X[,j], plq  [[i]]$X.quantiles[j,which.ci[1],], plq[[i]]$X.quantiles[j,which.ci[2],])
    CV[i,j,,2] = inside(sims[[i]]$X[,j], lw90q[[i]]$X.quantiles[j,which.ci[1],], plq[[i]]$X.quantiles[j,which.ci[2],])
    CV[i,j,,3] = inside(sims[[i]]$X[,j], lw95q[[i]]$X.quantiles[j,which.ci[1],], plq[[i]]$X.quantiles[j,which.ci[2],])
    CV[i,j,,4] = inside(sims[[i]]$X[,j], lw99q[[i]]$X.quantiles[j,which.ci[1],], plq[[i]]$X.quantiles[j,which.ci[2],])
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s+j
    CV[i,k,,1] = inside(sims[[i]]$probs[j], plq  [[i]]$p.quantiles[j,which.ci[1],], plq[[i]]$p.quantiles[j,which.ci[2],])
    CV[i,k,,2] = inside(sims[[i]]$probs[j], lw90q[[i]]$p.quantiles[j,which.ci[1],], plq[[i]]$p.quantiles[j,which.ci[2],])
    CV[i,k,,3] = inside(sims[[i]]$probs[j], lw95q[[i]]$p.quantiles[j,which.ci[1],], plq[[i]]$p.quantiles[j,which.ci[2],])
    CV[i,k,,4] = inside(sims[[i]]$probs[j], lw99q[[i]]$p.quantiles[j,which.ci[1],], plq[[i]]$p.quantiles[j,which.ci[2],])
  }
  for (j in 1:sckm$r)
  {
    k = sckm$s + sckm$r + j
    CV[i,k,,1] = inside(sims[[i]]$rates[j], plq  [[i]]$r.quantiles[j,which.ci[1],], plq[[i]]$r.quantiles[j,which.ci[2],])
    CV[i,k,,2] = inside(sims[[i]]$rates[j], lw90q[[i]]$r.quantiles[j,which.ci[1],], plq[[i]]$r.quantiles[j,which.ci[2],])
    CV[i,k,,3] = inside(sims[[i]]$rates[j], lw95q[[i]]$r.quantiles[j,which.ci[1],], plq[[i]]$r.quantiles[j,which.ci[2],])
    CV[i,k,,4] = inside(sims[[i]]$rates[j], lw99q[[i]]$r.quantiles[j,which.ci[1],], plq[[i]]$r.quantiles[j,which.ci[2],])
  }
}

MCV = apply(CV, 2:4, mean, na.rm=T)





MSE  = array.to.df(MSE ); MSE $statistic = "MSE"
MCV  = array.to.df(MCV ); MCV $statistic = "MCV"
MAD  = array.to.df(MAD ); MAD $statistic = "MAD"
MAPE = array.to.df(MAPE); MAPE$statistic = "MAPE"

write.csv(rbind(MSE,MCV,MAD,MAPE), file="summary-statistics.csv", row.names=F)


q(ifelse(interactive(),"ask","no"))



