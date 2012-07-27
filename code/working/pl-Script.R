source("inference.R")
source("SMLib.R");
require(smcUtils)
require(Hmisc)
source("pfSIR.R")

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
plot.ci(temp$stat[1,,],temp$trueX,temp$theta,1)
plot.ci.mcError(temp$stat,1)


#################################
# Try to handle the actual weekly data set from Harare: 59 weeks total
#################################
Nweek <- 59
harare <- read.csv("data/harare/harareClean.csv")
harY <- array(0,dim=c(Nweek,4))
harY[,1] <- diff(harare$total)
hararePL <- plSIR(5000,Nweek,LOOPN=1,Y=harY,verbose="HIST")
plot.ci(hararePL$stat[1,,],NULL,hararePL$theta,1)
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
N.week <- 52
X0 = c(16000,2,0,0)

load("../../data/SIRDsims/SIRDsims.RData")
load("../../data/SIRDsims/SIRDsims-mcmc-density.RData")
n.sims <- length(sims)

sim.PL <- array(0, dim=c(n.sims,N.week,24))
sim.LW <- array(0, dim=c(n.sims,N.week,24))
sim.SV <- array(0, dim=c(n.sims,N.week,24))
smc.PL <- array(0, dim=c(n.sims,12))
smc.LW <- array(0, dim=c(n.sims,12))
smc.SV <- array(0, dim=c(n.sims,12))
kd.PL <- list(); kd.LW <- list(); kd.SV <- list();
i.pl <- 1; i.lw <- 1; i.sv <- 1;

sims.params <- list(    initP=rep(0,N.RXNS),
    initX=X0,    hyperPrior = as.vector(rbind(prior$theta$a,prior$theta$b)),
    trueTheta = rep(0,N.RXNS) )


for (j in 1:1)
{
  sims.params$initP <- probs[j,]
  sims.params$trueTheta <- thetas[j,]
  simY <- as.matrix(sims[[j]]$y)
  simY[,3:4] <- 0
  simX <- as.matrix(sims[[j]]$x)
  
  simPL <- plSIR(10000,n[j],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  simLW <- particleSampledSIR(aLW=0.98,12000,n[j],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  simSV <- particleSampledSIR(aLW=2,12000,n[j],LOOPN=1,Y=simY,trueX=simX,verbose="C",model.params=sims.params)
  
  # store all the quantiles; the final quantiles of interest
  sim.PL[j,1:(n[j]+1),] <- simPL$stat[1,,]
  sim.LW[j,1:(n[j]+1),] <- simLW$stat[1,,]
  sim.SV[j,1:(n[j]+1),] <- simSV$stat[1,,]
  smc.PL[j,] <- as.vector(simPL$stat[1,n[j]+1,1:12])
  smc.LW[j,] <- as.vector(simLW$stat[1,n[j]+1,1:12])
  smc.SV[j,] <- as.vector(simSV$stat[1,n[j]+1,1:12])
  
  # and the final density estimates
  for (k in 1:4) {
     kd.PL[[i.pl]] <- simPL$density[[k]]
     i.pl <- i.pl + 1
     kd.LW[[i.lw]] <- simLW$density[[k]]
     i.lw <- i.lw + 1
     kd.SV[[i.sv]] <- simSV$density[[k]]
     i.sv <- i.sv + 1
  }
}

key = paste("simulation ",rep(1:n.sims,each=4),", ",c("S","I","S->I","I->R"),sep="")
save(kd.PL,key,file="../../data/SIRDsims/SIRDsims-smcPL-density.RData")
save(kd.LW,key,file="../../data/SIRDsims/SIRDsims-smcLW-density.RData")
save(kd.SV,key,file="../../data/SIRDsims/SIRDsims-smcSV-density.RData")

colnames(smc.PL) = paste(rep(c("S","SI","I","IR"),each=3),c("50","2.5","97.5"))
write.csv(smc.PL,"../../data/SIRDsims/SIRDsims-smcPL-quantiles.csv",row.names=F)
colnames(smc.LW) = paste(rep(c("S","SI","I","IR"),each=3),c("50","2.5","97.5"))
write.csv(smc.LW,"../../data/SIRDsims/SIRDsims-smcLW-quantiles.csv",row.names=F)
colnames(smc.SV) = paste(rep(c("S","SI","I","IR"),each=3),c("50","2.5","97.5"))
write.csv(smc.SV,"../../data/SIRDsims/SIRDsims-smcSV-quantiles.csv",row.names=F)

#####################################
# Make a plot to compare MCMC vs SMC
######################################
# Compare the sequential run for scenario 1
seq.mcmc <- read.csv("SIRDsims-mcmc-seq-quants.csv")
comp.quants <- array(0, dim=c(2,n[1],12))

comp.quants[2,,]<-as.matrix(seq.mcmc)
comp.quants[1,,]<-simPL$stat[1,2:(n[1]+1),1:12]
 plot.ci(comp.quants,simX[2:(n[1]+1),],sims.params$trueTheta,col=8,in.legend=c("PL","MCMC"),sir.plotCI=1,ltype=c(1,1))




######################################
# Compare terminal densities

mcmc.quantile <- read.csv("../../data/SIRDsims/SIRDsims-mcmc-quantiles.csv")
pl.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcPL-quantiles.csv")
lw.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcLW-quantiles.csv")
sv.quantile <- read.csv("../../data/SIRDsims/SIRDsims-smcSV-quantiles.csv")

plot.ndx <- 1 # which simulation to look at 

par(mfrow=c(2,2),mar=c(4,3,2,1), oma=c(0,0,0,1), mgp=c(1,1,0))
# Make 4 panels, each panel shows the 4 different algorithm outputs (in terms of terminal densities)

XLab <- c("S", "I", "S->I", "I->R")
stat.ndx <- c(1,7,4,10);
for (k in 1:4) {
  if (k==2) {
      plot( kd[[4*plot.ndx+k-4]], col=1, main=XLab[k], xlab="", yaxt="n", ylab="Posterior Density",ylim=c(0,20))
      lines( kd.PL[[4*plot.ndx+k-4]], col=2 , do.p="F", vert="T")
      lines( kd.LW[[4*plot.ndx+k-4]], col=3 , do.p="F", vert="T")
      lines( kd.SV[[4*plot.ndx+k-4]], col=4 , do.p="F", vert="T")
  }
  else {
        plot( kd[[4*plot.ndx+k-4]], col=1, main=XLab[k], xlab="", yaxt="n", ylab="Posterior Density")
        lines( kd.LW[[4*plot.ndx+k-4]], col=3)
        lines( kd.SV[[4*plot.ndx+k-4]], col=4)
        lines( kd.PL[[4*plot.ndx+k-4]], col=2)
  }
  abline(v=mcmc.quantile[plot.ndx,stat.ndx[k]],pch=4,col=1, cex=1.5)
  abline(v=pl.quantile[plot.ndx,stat.ndx[k]],pch=5,col=2, cex=1.5)
  abline(v=lw.quantile[plot.ndx,stat.ndx[k]],pch=6,col=3, cex=1.5)
  abline(v=sv.quantile[plot.ndx,stat.ndx[k]],pch=7,col=4, cex=1.5)
  if (k > 2 )
     legend("topright", c("MCMC", "PL", "Liu-West", "Storvik"), lty=c(1,1,1,1),col=1:4)
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
  plot.ci(stt,as.matrix(sims[[j]]$x[1:nn,]),thetas[j,],plot.diff=1,col=5,in.legend=c("PL","LW","Storvik"),sir.plotCI=1)
  #savePlot(filename=paste("../SIRDsims",j,".eps",sep=""), type="eps")
  #savePlot(filename=paste("../SIRDsims",j,".pdf",sep=""), type="pdf")
  plot.ci(stt,as.matrix(sims[[j]]$x[1:nn,]),thetas[j,],plot.diff=0,col=5,in.legend=c("PL","LW","Storvik"),sir.plotCI=1)
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

plot.ci(pl.prior,simX[1:(n[.ndx]+1),],thetas[.ndx,],1, col=c("#000000","#7777ff","#ff3333"), in.legend=c("Wide", "Low", "Narrow"),ltype=1:3)
savePlot(filename="20120723-prior-sensitivity", type="pdf")
