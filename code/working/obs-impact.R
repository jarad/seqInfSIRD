#############################################################
# Compare impact of observations on the inference
#############################################################

# Generate a scenario
N.week <- 60
scen.params <- base.params
scen.params$initP <- c(1,1,1,1) # use full observations first
set.seed(5)
sampScen <- generate.scenario(scen.params,tauLeap,N.week)
sampScen$theta
# [1] 0.7159145 0.5000000 0.0000000 0.0020000
plot(sampScen$X[,2]/2,type="l")
lines(sampScen$Y[,1],col="green")
scen.params$trueTheta <- sampScen$theta
.UNKNOWNP <- FALSE
N.parts <- 16000

#######################
# now generate various observations by sub-sampling
#  Full obs case
Yfull <- sampScen$Y
scen.params$initP <- c(1,1,1,1)
scenFull <- plSIR(N.parts,N.week,LOOPN=1,Y=Yfull,trueX=sampScen$X,verbose="CI",model.params=scen.params)

##########################
# Full obs of S->I only
Ybeta <- array(0,dim=c(N.week+1,4))
Ybeta[,1] <-sampScen$Y[,1]
scen.params$initP <- c(1,0,0,0)
scenBeta <- plSIR(N.parts,N.week,LOOPN=1,Y=Ybeta,trueX=sampScen$X,verbose="CI",model.params=scen.params)

##########################
# Typical observations of S->I
Ymed <- array(0,dim=c(N.week+1,4))
Ymed[,1] <-rbinom(N.week+1,sampScen$Y[,1],0.2)
scen.params$initP <- c(0.2,0,0,0)
scenMed <- plSIR(N.parts,N.week,LOOPN=1,Y=Ymed,trueX=sampScen$X,verbose="CI",model.params=scen.params)

##########################
# Typical low observations (like in Zimbabwe)
Ylow <- array(0,dim=c(N.week+1,4))
Ylow[,1] <-rbinom(N.week+1,sampScen$Y[,1],0.05)
scen.params$initP <- c(0.05,0,0,0)
scenLow <- plSIR(N.parts,N.week,LOOPN=1,Y=Ylow,trueX=sampScen$X,verbose="CI",model.params=scen.params)


##########################
# 100%-observations of I->R
YgamFull <- array(0,dim=c(N.week+1,4))
YgamFull[,2] <-sampScen$Y[,2]
scen.params$initP <- c(0,1,0,0)
scenGammaFull <- plSIR(N.parts,N.week,LOOPN=1,Y=YgamFull,trueX=sampScen$X,verbose="CI",model.params=scen.params)


##########################
# 20%-observations of each S->I and I->R
Ygam <- Ymed
Ygam[,2] <- rbinom(N.week+1,sampScen$Y[,2],0.2)
scen.params$initP <- c(0.2,0.2,0,0)
scenGamma <- plSIR(N.parts,N.week,LOOPN=1,Y=Ygam,trueX=sampScen$X,verbose="CI",model.params=scen.params)

##########################
# 5%-observations of S->I and I->R
YgamLow <- Ylow
YgamLow[,2] <- rbinom(N.week+1,sampScen$Y[,2],0.05)
scen.params$initP <- c(0.05,0.05,0,0)
scenGammaLow <- plSIR(N.parts,N.week,LOOPN=1,Y=YgamLow,trueX=sampScen$X,verbose="CI",model.params=scen.params)

##########################
# observe 20% of infecteds as well as ALL deaths (very implicit info about current infecteds)
Ydeaths <- Ymed
Ydeaths[,4] <- sampScen$Y[,4]
scen.params$initP <- c(0.2,0,0,1)
scenDeaths <- plSIR(N.parts,N.week,LOOPN=1,Y=Ydeaths,trueX=sampScen$X,verbose="CI",model.params=scen.params)



#########################################
# Impact of unknown P
.UNKNOWNP <- TRUE

# tight prior for P around true p=0.2
scen.params$hyperPrior <- c(22.5,30,15,30,0,1000,0.6,80,   10,40,10,40,0.01,10,0.01,10)
scen.params$initP <- c(0.2,0.2,0,0)
scenUncertainPBase <- plSIR(N.parts,N.week,LOOPN=1,Y=Ygam,trueX=sampScen$X,verbose="CI",model.params=scen.params)

# Uniform (very vague) prior for P
scen.params$hyperPrior <- c(22.5,30,15,30,0,1000,0.6,80,   1,1,1,1,0.01,10,0.01,10)
scen.params$initP <- c(0.2,0.2,0,0)
scenUncertainPUnif <- plSIR(N.parts,N.week,LOOPN=1,Y=Ygam,trueX=sampScen$X,verbose="CI",model.params=scen.params)

# vague prior for P
scen.params$hyperPrior <- c(22.5,30,15,30,0,1000,0.6,80,   1,5,1,5,0.01,10,0.01,10)
scen.params$initP <- c(0.2,0.2,0,0)
scenUncertainPVague <- plSIR(N.parts,N.week,LOOPN=1,Y=Ygam,trueX=sampScen$X,verbose="CI",model.params=scen.params)

save(sampScen,scenFull,scenBeta, scenMed,scenLow,scenGammaFull,scenGamma,scenGammaLow,scenDeaths,scenUncertainPBase,scenUncertainPUnif, file="SIRDsims-obs-impact.RData")
readline("hit enter")

###########################################
## Plot the resulting densities against each other to see the effect
stat.ndx <- list();
stat.ndx[[1]] <- 1:3;  stat.ndx[[2]] <- 7:9;
stat.ndx[[3]] <- 4:6;  stat.ndx[[4]] <- 10:12;


# Plot 1: Effect of changing S->I sampling
par(mfrow=c(2,2),mar=c(3,3,2,1)+0.1, oma=c(0,0,0,1),mgp=c(1,1,0))
XLab <- c("S", "I", "S->I", "I->R")
truth <- c( sampScen$X[N.week+1,1], sampScen$X[N.week+1,2], sampScen$theta[1], sampScen$theta[2])

for (k in 1:4) {
  plot( scenMed$density[[k]], col=1, lty=3, main=XLab[k], xlab="", yaxt="n", ylab="Posterior  Density", lwd=2, do.p="F")
  lines( scenBeta$density[[k]], col=2, lty=2, lwd=2, do.p="F")
  lines( scenFull$density[[k]], col=3, lty=1, lwd=2, do.p="F")
  lines( scenLow$density[[k]], col=4, lty=4, lwd=2, do.p="F")
  #points(scenMed$stat[1,N.week+1,stat.ndx[[k]] ],rep(0,3),pch=19,col=1, cex=1.5)
  #points(scenLow$stat[1,N.week+1,stat.ndx[[k]] ],rep(0,3),pch=19,col=2, cex=1.5)
  #points(scenGamma$stat[1,N.week+1,stat.ndx[[k]] ],rep(0,3),pch=19,col=3, cex=1.5)
  #points(scenDeaths$stat[1,N.week+1,stat.ndx[[k]] ],rep(0,3),pch=19,col=4, cex=1.5)

  abline(v=truth[k], col="red")
  if (k >2)
     legend("topright", c("Full Obs", "100% S->I ", "20% S->I", "5% S->I"), col=c(3,2,1,4),lty=c(1,2,3,4))
}


dev.new()
par(mfrow=c(2,2),mar=c(3,3,2,1)+0.1, oma=c(0,0,0,1),mgp=c(1,1,0))
for (k in 1:4) {
  plot( scenFull$density[[k]], col=1, main=XLab[k], xlab="", yaxt="n", ylab="Posterior  Density")
  lines( scenGamma$density[[k]], col=2)
  lines( scenGammaLow$density[[k]], col=3)
  lines( scenDeaths$density[[k]], col=4)
  #points(scenMed$stat[1,N.week+1,stat.ndx[[k]] ],rep(0,3),pch=19,col=1, cex=1.5)
  #points(scenLow$stat[1,N.week+1,stat.ndx[[k]] ],rep(0,3),pch=19,col=2, cex=1.5)
  #points(scenGamma$stat[1,N.week+1,stat.ndx[[k]] ],rep(0,3),pch=19,col=3, cex=1.5)
  #points(scenDeaths$stat[1,N.week+1,stat.ndx[[k]] ],rep(0,3),pch=19,col=4, cex=1.5)

  abline(v=truth[k])
  if (k >2)
     legend("topright", c("Full Obs", "20% S->I/I->R ", "5% S->I/I->R", "20% S->I, 100% I->D"), lty=c(1,1,1,1),col=1:4)
}


uncertain.P <- array(0,dim=c(4,61,24))
uncertain.P[1,,] <- scenGamma$stat[1,,]
uncertain.P[2,,] <- scenUncertainPBase$stat[1,,]
uncertain.P[3,,] <- scenUncertainPVague$stat[1,,]
uncertain.P[4,,] <- scenUncertainPUnif$stat[1,,]

plot.ci(uncertain.P,sampScen$X,scen.params$trueTheta,in.legend=c("Known P","Tight Prior", "Vague Prior", "Uniform Prior"),ltype=1,col=8)



par(mfrow=c(2,2),mar=c(3,3,2,1)+0.1, oma=c(0,0,0,1),mgp=c(1,1,0))
for (k in 1:4) {
  plot( scenGamma$density[[k]], col=1, main=XLab[k], xlab="", yaxt="n", ylab="Posterior  Density")
  lines( scenUncertainPBase$density[[k]], col=2)
  lines( scenUncertainPUnif$density[[k]], col=3)
  lines( scenUncertainPVague$density[[k]], col=4)

  abline(v=truth[k])
  if (k >2)
     legend("topright", c("Known p", "Tight Prior", "Uniform Prior", "Vague Prior"), lty=c(1,1,1,1),col=1:4)
}

#savePlot(filename="20120723-density-comp-p", type="pdf")
