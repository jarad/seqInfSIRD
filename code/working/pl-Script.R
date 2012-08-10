source("inference.R")
source("SMLib.R");
require(smcUtils)
require(Hmisc)
source("pfSIR.R")
source("plSIR-plot.R")

options(error=recover)
require(plotrix)

dyn.load(paste("gillespieExactStep-",.Platform$r_arch,.Platform$dynlib.ext,sep=''))
dyn.load(paste("inference-",.Platform$r_arch,.Platform$dynlib.ext,sep=''))

#base.params <- list(    initP=c(0.25, 0.25, 0.5, 0.5),
#    initX=c(19980,20, 0,0),    hyperPrior = c(22.5,30,15,30,0,1000,1.6,80),
#    trueTheta = array(c(0.75, 0.5, 0, 0.02),dim=c(4,1)) )

N.RXNS <- 4

# initP are the sampling proportions: S->I, I->R, S->R, I->D (in that order)
# initX is the initial state of S/I/R/D
# hyperPrior is the Gamma prior on the 4 thetas: (alpha_i,beta_i)_{i=1}^4
# trueTheta are the actual thetas for S->I, I->R, S->R, I->D
base.params <- list(    initP=c(0.05, 0, 0, 0),
    initX=c(15980,20, 0,0),    hyperPrior = c(22.5,30,15,30,0,1000,0.6,80),
    trueTheta = array(c(0.8, 0.5, 0, 0.002),dim=c(4,1)) )



# RUN particleSampledSIR
#
tt <- plSIR(5000, 50, LOOPN=5,verbose="CI")
lw <- particleSampledSIR(5000, 50, LOOPN=5,aLW=0.99,verbose="f",trueX=tt$trueX,Y=tt$Y)  # Liu-West

# compare different methods
tt <- plSIR(5000, 52, LOOPN=10,verbose="f")
lw <- particleSampledSIR(5000, 50, LOOPN=10,aLW=0.99,verbose="f",trueX=tt$trueX,Y=tt$Y) # liu-west
str <- particleSampledSIR(5000, 50, LOOPN=10,aLW=1.99,verbose="f",trueX=tt$trueX,Y=tt$Y) # storvik


plot(c(0,50), c(0.75,0.75), type="l",col="red",xlab="Periods t",ylab="Beta Posterior",lwd=2)
plot.ci.mcError(str$stat,1,col1="#ffeedd99", col2="black")
plot.ci.mcError(tt$stat,1,col1="#eeddff99", col2="blue")
plot.ci.mcError(lw$stat,1,col1="#ddffee99", col2="green")

#####################
# To see a more realistic example
temp <- plSIR(5000, 70, LOOPN=10,verbose="HIST")
plot.ci(temp$stat[1,,],temp$trueX,temp$theta)
plot.ci.mcError(temp$stat,1)


#################################
# Try to handle the actual weekly data set from Harare: 59 weeks total
#################################
Nweek <- 59
harare <- read.csv("../../data/harare/harareClean.csv")
harY <- array(0,dim=c(Nweek,4))
harY[,1] <- diff(harare$total)
harare.params <- list()

Nweek <- 30
harare.den.I <- list(); harare.den.SI <- list(); harare.den.IR <- list(); harare.den.SIprop <- list()
harare.quants <- array(0, dim=c(6,Nweek+1,24))


.UNKNOWNP <- F
harare.params$initP <- c(0.1,0,0,0)
harare.params$initX <- c(2000,5,0,0)

.ndx <-1
harare.params$hyperPrior <- c(120,100,100,100)
hararePL <- plSIR(10000,Nweek,Y=harY,verbose="CI",model.params=harare.params)
harare.quants[.ndx,,1:(3*(N.STATES+N.RXNS))] <- hararePL$stat
harare.den.I[[.ndx]] <- hararePL$density[[2]]
harare.den.SI[[.ndx]] <- hararePL$density[[3]]
harare.den.IR[[.ndx]] <- hararePL$density[[4]]

.ndx <-2
harare.params$hyperPrior <- c(3,2,100,100) # tight on I->R
hararePL <- plSIR(10000,Nweek,Y=harY,verbose="CI",model.params=harare.params)
harare.quants[.ndx,,1:(3*(N.STATES+N.RXNS))] <- hararePL$stat
harare.den.I[[.ndx]] <- hararePL$density[[2]]
harare.den.SI[[.ndx]] <- hararePL$density[[3]]
harare.den.IR[[.ndx]] <- hararePL$density[[4]]

.ndx <- 3
harare.params$hyperPrior <- c(3,2,4,4)
hararePL <- plSIR(10000,Nweek,LOOPN=1,Y=harY,verbose="CI",model.params=harare.params)
harare.quants[.ndx,,1:(3*(N.STATES+N.RXNS))] <- hararePL$stat
harare.den.I[[.ndx]] <- hararePL$density[[2]]
harare.den.SI[[.ndx]] <- hararePL$density[[3]]
harare.den.IR[[.ndx]] <- hararePL$density[[4]]


####
.UNKNOWNP <- T

.ndx <- 4
harare.params$hyperPrior <- c(120,120,100,100,1,1,1,1)   # tight theta, vague prop
hararePL <- plSIR(10000,Nweek,LOOPN=1,Y=harY,verbose="CI",model.params=harare.params)
harare.quants[.ndx,,] <- hararePL$stat
harare.den.I[[.ndx]] <- hararePL$density[[2]]
harare.den.SI[[.ndx]] <- hararePL$density[[3]]
harare.den.IR[[.ndx]] <- hararePL$density[[4]]
harare.den.SIprop[[.ndx]] <- hararePL$density[[5]]

.ndx <-5
harare.params$hyperPrior <- c(3,2,2,2,1,1,1,1)  # very vague priors all around
hararePL <- plSIR(10000,Nweek,LOOPN=1,Y=harY,verbose="CI",model.params=harare.params)
harare.quants[.ndx,,] <- hararePL$stat
harare.den.I[[.ndx]] <- hararePL$density[[2]]
harare.den.SI[[.ndx]] <- hararePL$density[[3]]
harare.den.IR[[.ndx]] <- hararePL$density[[4]]
harare.den.SIprop[[.ndx]] <- hararePL$density[[5]]

.ndx<-6
harare.params$hyperPrior <- c(3,2,100,100,1,5,0.1,5)  #tighter prior on prop
hararePL <- plSIR(10000,Nweek,LOOPN=1,Y=harY,verbose="CI",model.params=harare.params)
harare.quants[.ndx,,] <- hararePL$stat
harare.den.I[[.ndx]] <- hararePL$density[[2]]
harare.den.SI[[.ndx]] <- hararePL$density[[3]]
harare.den.IR[[.ndx]] <- hararePL$density[[4]]
harare.den.SIprop[[.ndx]] <- hararePL$density[[5]]

cols <-rainbow(7,alpha=0.4)
leg <- c("Tight","Vague S->I", "Both Vague","Tight + p","Vague +p","Tight on p")
par(mfrow=c(2,4))
plot.ci(harare.quants,which.dim=list(x=2,theta=c(1,2),prop=1),col=cols,in.legend=leg,shade=T)

cols <-rainbow(7)
plot( harare.den.I[[1]],do.p="F",col=cols[1],main="I",ylab="",xlab="",xlim=c(0,50))
legend("topright",leg,col=cols,lwd=2,lty=1)
for (kk in 2:6) 
   lines(harare.den.I[[kk]],col=cols[kk],do.p="F")

plot( harare.den.SI[[1]],type="l",col=cols[1],xlim=c(0.8,2),main="S->I",ylab="",xlab="")
for (kk in 2:6) 
   lines(harare.den.SI[[kk]],col=cols[kk])

plot( harare.den.IR[[1]],type="l",col=cols[1],xlim=c(0.8,1.9),main="I->R",ylab="",xlab="")
for (kk in 2:6) 
   lines(harare.den.IR[[kk]],col=cols[kk])
legend("topright",leg,col=cols,lwd=2,lty=1)     

plot( harare.den.SIprop[[4]],type="l",col=cols[4],xlim=c(0.04,0.36),main="S->I Prop",ylab="",xlab="")
for (kk in 5:6) 
   lines(harare.den.SIprop[[kk]],col=cols[kk])
# 



plot.ci(hararePL$stat[1,,],which.dim=2)
hararePL <- plSIR(5000,Nweek,LOOPN=5,Y=harY,verbose="CI")
plot(c(0,60),c(0.6,0.6),type="l",col="red",ylab="Beta posterior",xlab="Periods")
plot.ci.mcError(hararePL$stat,1,col1="#ddddddaa", col2="blue")

# The prior on GAMMA has a strong impact on the POSTERIOR of BETA:
base.params$hyperPrior <- c(1,1.2,20,40,0,1000,0,1000)
hararePL3 <- plSIR(5000,Nweek-1,LOOPN=10,Y=harY,verbose="CI")
plot(c(0,60),c(0.6,0.6),type="l",col="red",ylab="Beta posterior",xlab="Periods")
plot.ci.mcError(hararePL$stat,1,col1="#ddddddaa", col2="blue") # close to 0.5

base.params$hyperPrior <- c(2.25,3,1,2,0,1000,0,1000)  # close to 1 due to the wide prior of Gamma
hararePL <- plSIR(5000,Nweek-1,LOOPN=10,Y=harY,verbose="CI")


################################################
# To run SIRDsims.RData examples
################################################

load("../../data/SIRDsims/SIRDsims.RData")
n.sims <- length(sims)
N.week <- max(n) +1
STOICH_MATRIX <- array(0,dim=c(4,2))
STOICH_MATRIX[1:3,] <- stoich


sim.PL <- array(0, dim=c(2*n.sims,N.week,24))
sim.LW <- array(0, dim=c(2*n.sims,N.week,24))
sim.LW2 <- array(0, dim=c(2*n.sims,N.week,24))
sim.SV <- array(0, dim=c(2*n.sims,N.week,24))

smc.PL <- array(0, dim=c(2*n.sims,18))
smc.LW <- array(0, dim=c(2*n.sims,18))
smc.LW2 <- array(0, dim=c(2*n.sims,18))
smc.SV <- array(0, dim=c(2*n.sims,18))

kd.PL <- list(); kd.LW <- list(); kd.LW2 <- list(); kd.SV <- list();
i.pl <- 1; i.lw <- 1; i.sv <- 1; i.lw2 <- 1;
ln <- 3*(N.RXNS + N.STATES)
toSave <- c(1:6,13:18)
hp <- c(as.vector(cbind(prior$theta$a,prior$theta$b)),as.vector(cbind(prior$p$a,prior$p$b)))

sims.params <- list(    initP=rep(0,N.RXNS),
    initX=X0, trueTheta = rep(0,N.RXNS)
    )

offset <- n.sims
i.pl = i.lw = i.sv = i.lw2 = n.sims*4 + 1;
ln <- 3*(2*N.RXNS + N.STATES)
toSave <- c(1:6,13:24)
.UNKNOWNP <- T

for (j in 3:n.sims)
{
  sims.params$initP <- probs[j,]
  sims.params$trueTheta <- thetas[j,]
  simY <- as.matrix(sims[[j]]$y)
  simY[,3:4] <- 0
  simX <- as.matrix(sims[[j]]$x)
  sims.params$hyperPrior <- hp[1:(4*N.RXNS)]
  
  simPL <- plSIR(10000,n[j],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  simLW <- particleSampledSIR(aLW=0.99,10000,n[j],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  simLW2 <- particleSampledSIR(aLW=0.96,10000,n[j],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)

  simSV <- particleSampledSIR(aLW=2,10000,n[j],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  
  # store all the quantiles; the final quantiles of interest
  sim.PL[j+offset,1:(n[j]+1),1:ln] <- simPL$stat[1,,]
  sim.LW[j+offset,1:(n[j]+1),1:ln] <- simLW$stat[1,,]
  sim.LW2[j+offset,1:(n[j]+1),1:ln] <- simLW2$stat[1,,]
  sim.SV[j+offset,1:(n[j]+1),1:ln] <- simSV$stat[1,,]
  
  smc.PL[j+offset,1:length(toSave)] <- as.vector(simPL$stat[1,n[j]+1,toSave])
  smc.LW[j+offset,1:length(toSave)] <- as.vector(simLW$stat[1,n[j]+1,toSave])
  smc.LW2[j+offset,1:length(toSave)] <- as.vector(simLW2$stat[1,n[j]+1,toSave])
  smc.SV[j+offset,1:length(toSave)] <- as.vector(simSV$stat[1,n[j]+1,toSave])
  
  # and the final density estimates
  for (k in 1:6) {
     kd.PL[[i.pl]] <- simPL$density[[k]]
     i.pl <- i.pl + 1
     kd.LW[[i.lw]] <- simLW$density[[k]]
     i.lw <- i.lw + 1
     kd.LW2[[i.lw2]] <- simLW2$density[[k]]
     i.lw2 <- i.lw2 + 1

     kd.SV[[i.sv]] <- simSV$density[[k]]
     i.sv <- i.sv + 1
  }
}

key = paste("simulation ",rep(1:n.sims,each=4),", ",c("S","I","S->I","I->R"),sep="")
save(kd.PL,key,file="../../data/SIRDsims/SIRDsims-smcPL-density.RData")
save(kd.LW,key,file="../../data/SIRDsims/SIRDsims-smcLW99-density.RData")
save(kd.LW2,key,file="../../data/SIRDsims/SIRDsims-smcLW96-density.RData")
save(kd.SV,key,file="../../data/SIRDsims/SIRDsims-smcSV-density.RData")

colnames(smc.PL) = paste(rep(c("S","SI","I","IR"),each=3),c("50","2.5","97.5"))
write.csv(smc.PL,"../../data/SIRDsims/SIRDsims-smcPL-quantiles.csv",row.names=F)
colnames(smc.LW) = paste(rep(c("S","SI","I","IR"),each=3),c("50","2.5","97.5"))
write.csv(smc.LW,"../../data/SIRDsims/SIRDsims-smcLW99-quantiles.csv",row.names=F)
colnames(smc.LW2) = paste(rep(c("S","SI","I","IR"),each=3),c("50","2.5","97.5"))
write.csv(smc.LW2,"../../data/SIRDsims/SIRDsims-smcLW96-quantiles.csv",row.names=F)

colnames(smc.SV) = paste(rep(c("S","SI","I","IR"),each=3),c("50","2.5","97.5"))
write.csv(smc.SV,"../../data/SIRDsims/SIRDsims-smcSV-quantiles.csv",row.names=F)

#####################################
# Make a plot to compare MCMC vs SMC
######################################
# Compare the sequential run for scenario 1
seq.mcmc <- read.csv("SIRDsims-mcmc-seq-quants.csv")
load("../../data/SIRDsims/SIRDsims-mcmc-density.RData")

comp.quants <- array(0, dim=c(2,n[1],12))

comp.quants[2,,]<-as.matrix(seq.mcmc)
comp.quants[1,,]<-simPL$stat[1,2:(n[1]+1),1:12]
 plot.ci(comp.quants,simX[2:(n[1]+1),],sims.params$trueTheta,col=8,in.legend=c("PL","MCMC"),ltype=c(1,1))

#########################################
# Create a "hash" plot comparing CIs across methods/scenarios

which.col <- c("IR 2.5", "IR 97.5")
plot(smc.PL[1,which.col], c(0,0), type="l", lwd=2.5, ylim=c(-0.5,25),xlim=c(0.75,1.25),col="#ff0000",yaxt="n",
  xlab="S->I Posterior", ylab="Scenarios", main="I->R")
cl = c("#ff0000","#000000","#777777","blue")

for (j in 1:n.sims) {
  lines(smc.PL[j,which.col], c(j-1+0,j-1+0), lwd=2.5, col=cl[1])
  lines(smc.LW[j,which.col], c(j-1-0.2,j-1-0.2), lwd=2, col=cl[2])
  lines(smc.LW2[j,which.col], c(j-1-0.1,j-1-0.1), lwd=2, col=cl[3])
  lines(smc.SV[j,which.col], c(j-1+0.1,j-1+0.1), lwd=2, col=cl[4])
}
legend("topleft", c("PL", "Liu-West", "LW96", "Storvik"), lty=rep(1,4),col=cl,lwd=rep(2.5,4))



######################################
# Compare terminal densities

mcmc.quantile <- read.csv("../../data/SIRDsims/SIRDsims-mcmc-quantiles.csv")
pl.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcPL-quantiles.csv")
lw.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcLW99-quantiles.csv")
lw2.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcLW96-quantiles.csv")
sv.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcSV-quantiles.csv")

plot.ndx <- 1 # which simulation to look at 

par(mfrow=c(1,2),mar=c(4,3,2,1), oma=c(0,0,0,1), mgp=c(1,1,0))
# Make 2 panels, each panel shows the 5 different algorithm outputs (in terms of terminal densities)

cols <- terrain.colors(6)
XLab <- c("S", "I", "S->I", "I->R")
stat.ndx <- c(1,7,4,10);
for (k in 3:4) {
  if (k==2) {
      plot( kd[[4*plot.ndx+k-4]], col=cols[1], main=XLab[k], xlab="", yaxt="n", ylab="Posterior Density",ylim=c(0,20))
      lines( kd.LW2[[4*plot.ndx+k-4]], col=cols[2] , do.p="F", vert="T")
      lines( kd.LW[[4*plot.ndx+k-4]], col=cols[3] , do.p="F", vert="T")
      lines( kd.SV[[4*plot.ndx+k-4]], col=cols[4] , do.p="F", vert="T")
       lines( kd.PL[[4*plot.ndx+k-4]], col=cols[5] , do.p="F", vert="T")
  
  }
  else {
        plot( kd[[4*plot.ndx+k-4]], col=cols[1], main=XLab[k], xlab="", yaxt="n", ylab="Posterior Density",lwd=2)
         lines( kd.LW2[[4*plot.ndx+k-4]], col=cols[2], lwd=2)
        lines( kd.LW[[4*plot.ndx+k-4]], col=cols[3], lwd=2)
        lines( kd.SV[[4*plot.ndx+k-4]], col=cols[4], lwd=2)
        lines( kd.PL[[4*plot.ndx+k-4]], col=cols[5], lwd=2)
  }
  abline(v=mcmc.quantile[plot.ndx,stat.ndx[k]],pch=3,col=cols[1], cex=1.5)
  abline(v=pl.quantile[plot.ndx,stat.ndx[k]],pch=4,col=cols[5], cex=1.5)
  abline(v=lw2.quantile[plot.ndx,stat.ndx[k]],pch=5,col=cols[2], cex=1.5)
  abline(v=lw.quantile[plot.ndx,stat.ndx[k]],pch=6,col=cols[3], cex=1.5)
  abline(v=sv.quantile[plot.ndx,stat.ndx[k]],pch=7,col=cols[4], cex=1.5)
  if (k > 2 )
     legend("topright", c("MCMC", "LW96", "Liu-West", "Storvik", "PL"), lty=rep(1,5),col=cols[1:5])
}
savePlot(filename=paste("../../data/SIRDsims/density-comp",plot.ndx,".pdf",sep=""), type="pdf")


diff.pl <- array(0,dim=c(16,8))
diff.lw <- array(0,dim=c(16,8))
diff.sv <- array(0,dim=c(16,8))
for (j in 1:n.sims) {
  for (k in 1:4) {
   diff.pl[j,2*k-1] <- mcmc.quantile[j,3*k-2]-pl.quantile[j,3*k-2]
   diff.pl[j,2*k] <- (pl.quantile[j,3*k]-pl.quantile[j,3*k-1])/(mcmc.quantile[j,3*k]-mcmc.quantile[j,3*k-1])
   diff.lw[j,2*k-1] <- mcmc.quantile[j,3*k-2]-lw.quantile[j,3*k-2]
   diff.lw[j,2*k] <- (lw.quantile[j,3*k]-lw.quantile[j,3*k-1])/(mcmc.quantile[j,3*k]-mcmc.quantile[j,3*k-1])
   diff.sv[j,2*k-1] <- mcmc.quantile[j,3*k-2]-sv.quantile[j,3*k-2]
   diff.sv[j,2*k] <- (sv.quantile[j,3*k]-sv.quantile[j,3*k-1])/(mcmc.quantile[j,3*k]-mcmc.quantile[j,3*k-1])
  }
}
apply(diff.pl, 2, mean)
apply(diff.lw, 2, mean)
apply(diff.sv, 2, mean)


###################################
# Save all the plots
for (j in 1:1) {
  nn <- n[j]+1
  stt <- array(0, dim=c(3,nn,24))
  stt[1,,] <- sim.PL[j,1:nn,]; stt[2,1:nn,] <- sim.LW[j,1:nn,]; stt[3,,] <-  sim.SV[j,1:nn,]
  plot.ci(stt,as.matrix(sims[[j]]$x[1:nn,]),thetas[j,],plot.diff=1,col=2,in.legend=c("PL","LW","Storvik"))
  #savePlot(filename=paste("../SIRDsims",j,".eps",sep=""), type="eps")
  #savePlot(filename=paste("../SIRDsims",j,".pdf",sep=""), type="pdf")
  plot.ci(stt,as.matrix(sims[[j]]$x[1:nn,]),thetas[j,],plot.diff=0,col=5,in.legend=c("PL","LW","Storvik"))
  #savePlot(filename=paste("../SIRDsims",j,"Abs.eps",sep=""), type="eps")
  #savePlot(filename=paste("../SIRDsims",j,"Abs.pdf",sep=""), type="pdf")
  
}

###########################################
# Compare the number of particles impact
###########################################
quants.PL <- array(0, dim=c(4,12))
quants.LW <- array(0, dim=c(4,12))
quants.SV <- array(0, dim=c(4,12))
kdp.PL <- list(); kdp.LW <- list(); kdp.SV <- list();
i.pl <- 1; i.lw <- 1; i.sv <- 1;
.ndx <- 12


for (j in 1:4)
{
  sims.params$initP <- probs[.ndx,]
  sims.params$trueTheta <- thetas[.ndx,]
  simY <- as.matrix(sims[[.ndx]]$y)
  simY[,3:4] <- 0
  simX <- as.matrix(sims[[.ndx]]$x)
  
  # keep doubling number of particles each time
  simPL <- plSIR(2000*2^j,n[.ndx],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  simLW <- particleSampledSIR(2000*2^j,n[.ndx],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  simSV <- particleSampledSIR(aLW=2,2000*2^j,n[.ndx],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  
  quants.PL[j,] <- as.vector(simPL$stat[1,n[.ndx]+1,1:12])
  quants.LW[j,] <- as.vector(simLW$stat[1,n[.ndx]+1,1:12])
  quants.SV[j,] <- as.vector(simSV$stat[1,n[.ndx]+1,1:12])
  
  for (k in 1:4) {
     kdp.PL[[i.pl]] <- simPL$density[[k]]
     i.pl <- i.pl + 1
     kdp.LW[[i.lw]] <- simLW$density[[k]]
     i.lw <- i.lw + 1
     kdp.SV[[i.sv]] <- simSV$density[[k]]
     i.sv <- i.sv + 1
  }
}

#key = paste("number of particles 2000, 4000, 8000, 16000",c("S","I","S->I","I->R"),sep="")
#save(kdp.PL,kdp.LW, kdp.SV, key,quants.PL, quants.LW, quants.SV, file="../../data/SIRDsims/SIRDsims-smc-noParts-density.RData")

# plot the densities for PL output as an example
par(mfcol=c(2,2),mar=c(4,4,2,1), oma=c(0,0,0,1))
XLab <- c("S", "I", "S->I", "I->R")
stat.ndx <- c(1,7,4,10);
for (k in 1:4) {
  plot( kd[[4*.ndx+k-4]], col=1, main=XLab[k], xlab="")
  lines( kdp.LW[[4+k-4]], col=2)
  lines( kdp.LW[[8+k-4]], col=3)
  lines( kdp.LW[[12+k-4]], col=4)
  lines( kdp.LW[[16+k-4]], col=5)
  if (k==2)
     legend("topright", c("MCMC", "2000", "4000", "8000","16000"), lty=c(1,1,1,1,1),col=1:5)
}
savePlot(filename=paste("../../data/SIRDsims/density-comp-parts-PL",.ndx,".pdf",sep=""), type="pdf")


###################################
# Compare the impact of the prior

sims.params <- list(    initP=rep(0,N.RXNS),
    initX=X0,    hyperPrior = as.vector(rbind(prior$theta$a,prior$theta$b))/2,
    trueTheta = rep(0,N.RXNS) )

.ndx <- 1
simY <- as.matrix(sims[[.ndx]]$y)
simX <- as.matrix(sims[[.ndx]]$x)
sims.params$trueTheta <- thetas[.ndx,]
sims.params$initP <- probs[.ndx,]

pl.prior <- array(0, dim=c(4,n[.ndx]+1,24))
for (j in 1:3)
{
  
  simPL <- plSIR(10000,n[.ndx],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  
  # store all the quantiles; the final quantiles of interest
  pl.prior[j,1:(n[.ndx]+1),] <- simPL$stat[1,,]
  
  sims.params$hyperPrior <- sims.params$hyperPrior*2

}

# mis-specified S-> I
sims.params$hyperPrior <- as.vector(rbind(prior$theta$a,prior$theta$b))
sims.params$hyperPrior[1] <- 100  # rather than 150
simPL <- plSIR(10000,n[.ndx],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  
pl.prior[2,1:(n[.ndx]+1),] <- simPL$stat[1,,]

plot.ci(pl.prior,simX[1:(n[.ndx]+1),],thetas[.ndx,], col=c("#000000","#7777ff","#ff3333"), in.legend=c("Wide", "Low", "Narrow"),ltype=1:3)
savePlot(filename="20120723-prior-sensitivity", type="pdf")
