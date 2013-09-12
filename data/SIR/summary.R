# Read in truth, means, and quantiles and 
# calculate MSE, MAD, MAPE, coverage
#
load("sims.RData")
source("settings.R")

states = sys$states
rxns   = sys$reactions
ns     = sys$s
nr     = sys$r
comp = c(states, paste("p:",rxns), paste("r:",rxns))


nn = c(length(probs), ns+2*nr, n+1)

## Convert to data frames for saving
library(reshape2)
array_to_df = function(arr) 
{
  df = melt(arr, varnames=c("parameter","time"))
  df$parameter = as.factor(comp[df$parameter])
  df$time      = df$time-1
  return(df)
}


#########################################################



# MSE
which.md = which(probs==0.5)
SE = array(NA, dim=nn)
for (i in 1:nn[1])
{
  for (j in 1:sys$s) 
  {
    SE[i,j,] = (sims[[i]]$X[   ,j] - q[[i]]$X.quantiles[j,which.md,])^2
  }  

  for (j in 1:sys$r)
  {
    k = sys$s + j
    SE[i,k,] = (sims[[i]]$probs[j] - q[[i]]$p.quantiles[j,which.md,])^2
  }
  for (j in 1:sys$r)
  {
    k = sys$s + sys$r + j
    SE[i,k,] = (sims[[i]]$rates[j] - q[[i]]$r.quantiles[j,which.md,])^2
  }
}

MSE = apply(SE, 2:3, mean, na.rm=T)
seMSE = sqrt(apply(SE, 2:3, sd, na.rm=T)/apply(SE, 2:3, function(x) length(na.omit(x))))


# MAD
which.md = which(probs==0.5)
AD = array(NA, dim=nn)
for (i in 1:nn[1])
{
  for (j in 1:sys$s) 
  {
    AD[i,j,] = abs(sims[[i]]$X[   ,j] - q[[i]]$X.quantiles[j,which.md,])
  }  

  for (j in 1:sys$r)
  {
    k = sys$s + j
    AD[i,k,] = abs(sims[[i]]$probs[j] - q[[i]]$p.quantiles[j,which.md,])
  }
  for (j in 1:sys$r)
  {
    k = sys$s + sys$r + j
    AD[i,k,] = abs(sims[[i]]$rates[j] - q[[i]]$r.quantiles[j,which.md,])
  }
}


MAD = apply(AD, 2:3, mean, na.rm=T)
seMAD = sqrt(apply(AD, 2:3, sd, na.rm=T)/apply(AD, 2:3, function(x) length(na.omit(x))))



# MAPE
which.md = which(probs==0.5)
APE = array(NA, dim=nn)
for (i in 1:nn[1])
{
  for (j in 1:sys$s)
  {
    APE[i,j,] = AD[i,j,]/(sims[[i]]$X[,j]+1) # To ensure finite results
  }
  for (j in 1:sys$r)
  {
    k = sys$s + j
    APE[i,k,] = AD[i,k,]/sims[[i]]$probs[j]    
  }
  for (j in 1:sys$r)
  {
    k = sys$s + sys$r + j
    APE[i,k,] = AD[i,k,]/sims[[i]]$rates[j]
  }
}

MAPE = apply(APE, 2:3, mean, na.rm=T)
seMAPE = sqrt(apply(APE, 2:3, sd, na.rm=T)/apply(APE, 2:3, function(x) length(na.omit(x))))






# Coverage
which.ci = c(1,5)
CV = array(NA, dim=nn)

vIsTRUE = Vectorize(isTRUE)
inside = function(x,a,b) return(ifelse(vIsTRUE(x>=a & x<=b), TRUE, FALSE))
for (i in 1:nn[1])
{
  for (j in 1:sys$s)
  {
    CV[i,j,] = inside(sims[[i]]$X[     ,j], q[[i]]$X.quantiles[j,which.ci[1],],
                                            q[[i]]$X.quantiles[j,which.ci[2],])
  }
  for (j in 1:sys$r)
  {
    k = sys$s+j
    CV[i,k,] = inside(sims[[i]]$probs[j], q[[i]]$p.quantiles[j,which.ci[1],],
                                          q[[i]]$p.quantiles[j,which.ci[2],])
  }
  for (j in 1:sys$r)
  {
    k = sys$s + sys$r + j
    CV[i,k,] = inside(sims[[i]]$rates[j], q[[i]]$r.quantiles[j,which.ci[1],],
                                          q[[i]]$r.quantiles[j,which.ci[2],])
  }
}

MCV = apply(CV, 2:3, mean, na.rm=T)
seMCV = sqrt(MCV*(1-MCV)/apply(CV, 2:3, function(x) length(na.omit(x))))

  MSE  = array_to_df(  MSE );   MSE $statistic =   "MSE"
seMSE  = array_to_df(seMSE ); seMSE $statistic = "seMSE"
  MCV  = array_to_df(  MCV );   MCV $statistic =   "MCV"
seMCV  = array_to_df(seMCV ); seMCV $statistic = "seMCV"
  MAD  = array_to_df(  MAD );   MAD $statistic =   "MAD"
seMAD  = array_to_df(seMAD ); seMAD $statistic = "seMAD"
  MAPE = array_to_df(  MAPE);   MAPE$statistic =   "MAPE"
seMAPE = array_to_df(seMAPE); seMAPE$statistic = "seMAPE"

df = rbind(MSE,seMSE,MCV,seMCV,MAD,seMAD,MAPE,seMAPE)
df$method = as.factor(method)
df$statistic = as.factor(df$statistic)


write.csv(df, file=paste(method,"sum.csv", sep=""), row.names=F)


q(ifelse(interactive(),"ask","no"))



