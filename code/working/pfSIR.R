dyn.load(paste("gillespieExactStep-",.Platform$r_arch,.Platform$dynlib.ext,sep=''))
dyn.load(paste("inference-",.Platform$r_arch,.Platform$dynlib.ext,sep=''))

#base.params <- list(    initP=c(0.25, 0.25, 0.5, 0.5),
#    initX=c(19980,20, 0,0),    hyperPrior = c(22.5,30,15,30,0,1000,1.6,80),
#    trueTheta = array(c(0.75, 0.5, 0, 0.02),dim=c(4,1)) )

N.RXNS <- 4

# initP are the sampling proportions: S->I, I->R, S->R, I->D (in that order
# initX is the initial state of S/I/R/D
# hyperPrior is the Gamma prior on the 4 thetas: (alpha_i,beta_i)_{i=1}^4
# trueTheta are the actual thetas for S->I, I->R, S->R, I->D
base.params <- list(    initP=c(0.2, 0.2, 0.5, 0.5),
    initX=c(24980,20, 0,0),    hyperPrior = c(22.5,30,15,30,0,1000,0.6,80),
    trueTheta = array(c(0.75, 0.5, 0, 0.002),dim=c(4,1)) )

##########################################################################
#Use particles to filter for SIR model with discrete-time sampled 
#transitions 
#
#  model.params are the initialization parameters -- see the default case base.params above
#  N is the number of particles
#  T is the horizon in periods of length dt (eg T=50,dt=1)
#  LOOPN is the number of repetitions to do
#  aLW is the Liu-West weight alpha; if aLW > 1 then Storvik algorithm is done instead
#  verbose=CI/HIST will output some summary plots using the given color
#  if trueX/Y is provided, that would be used as the true scenario
#  model.propagate.func is the method to use for particle mutation/true system evolution: tauLeap/gillespie
#  resample.func is the resampling method: multinomial.resample/branchMinVar
#
particleSampledSIR <- function(N, T, dt=1, model.params=base.params, LOOPN=1,aLW=0.99,
                verbose="CI",col="blue",trueX=NULL,Y=NULL,
                model.propagate.func=tauLeap,resample.func=multinomial.resample)
{
    initP <- model.params$initP
    initX <- model.params$initX
    hyperPrior <- model.params$hyperPrior
    trueTheta <- model.params$trueTheta


    saved.stats <- array(0, dim=c(LOOPN,ceil(T/dt)+1,24))  # all the saved summary statistics

    if (is.null(trueX)) # generate a scenario on the fly
    {
       trueX <- array(0, dim=c(ceil(T/dt)+1,4))
       Y <- array(0,dim=c(ceil(T/dt)+1,4))
       
       trueX[1,] <- initX

       # Fix the true state trajectory for replicability
       #set.seed(21)
       # Construct the true process and the observations
       for (i in 1:(T/dt)) 
       {
           out <- model.propagate.func(t(t(trueX[i,])),t(t(trueTheta)),t(initP))
           
           Y[i+1,] <- out$Y
           trueX[i+1,1:dim(out$X)[1]] <- out$X
       }
       #browser()
    }
    
    for (loop in 1:LOOPN) {    # run the filter LOOPN times to understand MC variance if needed
    
       # Initialize particles
       X <- t(array(rep(initX, N),dim=c(N.RXNS,N)))
       lambda <- t(array(rep(initP, N), dim=c(N.RXNS,N)))
       Suff <- t(array(rep(hyperPrior,N), dim=c(2*N.RXNS,N))) # sufficient conjugate statistics
       fixTheta <- array(0,dim=c(N,N.RXNS))

       # fixed theta particles if implementing Liu-West
       for (jj in 1:N.RXNS)
          fixTheta[,jj] <- rgamma(N, hyperPrior[2*jj-1],hyperPrior[2*jj])
       
       theta <- fixTheta
       p.weights <- array(1, N)/N

       # start at t=0
       curt <- 0 
       i <- 1
       totalWeight <- 1
       saved.stats[loop,i,] <- saveStats(X, fixTheta)

       # continue until there are some infecteds or up to T
       while (curt < T-dt & trueX[i,2] >0) {
         if (aLW > 1 & is.null(Suff) == 0 && i %% 3 == 1)  # Storvik filter: sample from the posterior Gamma mixture
           for (jj in 1:N.RXNS)
               theta[,jj] <- rgamma(N,Suff[,2*jj-1],Suff[,2*jj])
           
         else if (aLW < 1) # use the Liu-West thetas
            theta <- fixTheta
        
        # propagate particles 
        out <- model.propagate.func(t(X), t(theta),lambda,curY=Y[i+1,],hyper=Suff)
        X <- t(out$X); dX <- t(out$dX); Suff <- out$hyper
        
        # update weights
        p.weights <- updateWeights(dX,Y[i+1,],lambda,p.weights)
        curt <- (i-1)*dt

         if ( i %%15  == 6 & verbose =="HIST") # plot posterior of infectiousness parameter every 15 steps
         {
            
            if (aLW < 1)  # for Liu and West just histogram the particles
              hist(theta[,1],35,freq=F,main=sprintf("t=%d",curt),xlab='S->I Rate') 
            else { 
               # for Storvik construct exact Gamma pdf on a grid using S-particles
               gridx <- seq(0.4,1.1,by=0.0025)
               gridy <- array(0, length(gridx))
               for (jj in 1:length(gridx))
                  gridy[jj] <- sum(dgamma(gridx[jj],Suff[,1],Suff[,2]))
               gridy <- gridy/sum(gridy)
               print(c("i=",i-1,gridx[min(which(cumsum(gridy)>0.025))], gridx[max(which(cumsum(gridy)<0.975))]))

               plot(gridx,gridy,type="l",col=col,main=sprintf("t=%d",curt),xlab='S->I Rate')
            }
             
            abline(v=trueTheta[1], col="red")
            #browser()
         }
         totalWeight <- totalWeight*sum(p.weights)

         ##### resample and update all the particles
         ESS <- 1/sum(N^2*p.weights^2) 
       #if (i %% 3 == 0 | ESS < essThreshold) {
         newIndex <- resample.func(p.weights) 
         
         X <- X[newIndex,] 
         if (is.null(Suff) == 0)
             Suff <- Suff[newIndex,]
         lambda <- lambda[newIndex,]
         fixTheta <- fixTheta[newIndex,]
         theta <- theta[newIndex,]
         
         #### Liu and West move
         meanLam <- varLam <- array(0,N.RXNS)
         if (aLW < 1) # & curt > 1
         {
            for (jj in 1:4) {
              meanLam[jj] <- mean(fixTheta[,jj])
              varLam[jj] <- var(fixTheta[,jj])
              fixTheta[,jj] <- pmax(0, aLW*fixTheta[,jj] + (1-aLW)*meanLam[jj] +
                  sqrt( (1-aLW^2)*varLam[jj])*rnorm(N))
            }
          }
          p.weights <- rep(1/N, N)
       #} 

        i <- i+1
        saved.stats[loop,i,] <- saveStats(X,theta)
        # end of main loop
      }
   }
   
   if (verbose=='CI')  # plot some CI over time 
      plot.ci(saved.stats[1,1:i,],trueX[1:i,],trueTheta,1,col)
 
   return( list(stat=saved.stats,trueX=trueX,Y=Y))
}


#####################################################
# move one step of the SIR as a continuous-time Markov chain
#
# In this case there are no sufficient hyper-parameters, so last parameter is never used
# X is a matrix: each column has 4 rows for SIRD states 
# theta is a matrix, each column has 4 rows for SIRD rates 
# 
gillespieStep <- function(X, theta, prop,curY = NULL,hyper=NULL)
{
    hyper <- array(1, dim=c(4,dim(X)[2]))
        
    out <- gillespieExactStep.C(sir0=matrix(X[1:3,],nrow=3),prior=hyper,th=theta[1:3,])
    dX <- array(0, dim=c(4,dim(X)[2]))
    dX[1,] <- X[1,]-out$newSim[1,]
    dX[2,] <- dX[1,]
    dX[3:4,] <- 0
       
    if (is.null(curY)) 
        Y <- c(rbinom(2, dX[1:2,], prop[,1:2]),0,0)
    else
        Y <- curY

    return(list(X=out$newSim,dX=dX,Y=Y, hyper=NULL))
}



#######################################################
# move one step of the SIR as a tau-leaping Poisson approx
#
# X is a matrix: each column has 4 rows for SIRD states 
# theta is a matrix, each column has 4 rows for SIRD rates 
tauLeap<- function(X, theta, prop, curY=NULL, hyper=NULL)
{
    h <- t(X)
    dX <- X
    hyper2 <- array(1, dim=c(N.RXNS,2,dim(X)[2]))
    
    out <- one.step.C(X, hyper2, theta, t(prop), sample=T)
    
    for (j in 1:dim(X)[2])
    {    
       #out <- simulate.one.step.C(X[,j],theta[,j],sum(X[,j]) )
       #dX[,j] <- out$dX
       #X[,j] <- out$X
       h[j,]  <- hazard.R(X[,j], sum(X[,j]))   
    }
    newX <- out$newX[,1]
    newX[is.nan(newX)] <- 0  # to take care of the case when some proportions are zero
       
    if (is.null(curY))
        Y <- rbinom(N.RXNS, newX, prop)
    else
        Y <- curY
    #Y[is.nan(Y)] <- 0

    # update the hyperparameters
    if (!is.null(hyper)) 
      for (i in 1:N.RXNS) {
        hyper[,2*i-1] <- hyper[,2*i-1] + Y[i]
        hyper[,2*i] <- hyper[,2*i] + prop[,i]* h[,i]
      }
    
    return(list(X=out$X,dX=out$newX,hyper=hyper,Y=Y))
}

#######################################################
# Update the weights of the particles using the binomial sampling
updateWeights <- function(dX, curY, prop, weights)
{
    new.weights <- weights
   
    for (jj in 1:N.RXNS)
      new.weights <- new.weights*dbinom(curY[jj], dX[,jj], prop[,jj])
      
    return(new.weights)

}

#######################################################
# Predictive likelihood of the next observation using Poisson approximations
predictiveLikelihood <- function(X, nextY, theta, prop, weights)
{
    h <- X  # just to set the size of h correctly
    
    new.weights <- weights
    for (j in 1:dim(X)[1])
       h[j,]  <- hazard.R(X[j,], sum(X[j,]))   
       
    for (jj in 1:N.RXNS) {
        new.weights <- new.weights*dpois(nextY[jj],prop[,jj]*theta[,jj]*h[,jj])      
    }
    #browser()
    return(new.weights)

}

#######################################################
# Summary statistics: for each X-coordinate and theta parameter
# save the 95% CI and the mean generated by the particle cloud
saveStats <- function(X, theta)
{
    summ.stat <- array(0,dim=c(24,1))
    for (jj in 1:dim(X)[2])
    {
      summ.stat[(jj-1)*6+1] <- mean(X[,jj])
      summ.stat[((jj-1)*6+2):((jj-1)*6+3)] <- quantile(X[,jj],c(0.025,0.975))

      summ.stat[(jj-1)*6+4] <- mean(theta[,jj])
      summ.stat[((jj-1)*6+5):((jj-1)*6+6)] <- quantile(theta[,jj],c(0.025,0.975))
    }
    return(summ.stat)
}

################################################# 
# Resample using Crisan min-variance method
# branching from Crisan (2006) p .10
# similar to residual sampling but even less variance
branchMinVar <- function(p.weights)
{
        N <- length(p.weights)
        branch <- p.weights/mean(p.weights)
        fracBranch <- branch - floor(branch)
        ub <- runif(N-1)
        newNdx <- array(N, dim=c(N,1))
        curNdx <- 1
        gb <- N
        hb <- N
        for (j in 1:(N-1))
        {
          if (fracBranch[j] + gb - branch[j] - floor(gb - branch[j]) < 1)
            ob <- floor(branch[j]) + (ub[j] > 1-fracBranch[j]/(gb-floor(gb)+1e-8))*(hb - floor(gb))
          else
            ob <- floor(branch[j]) + (hb - floor(gb)) + (1-hb+floor(gb))*(ub[j] < 1 - (1-fracBranch[j])/(1-gb+floor(gb)))

          gb <- gb - branch[j]
          if (ob == 0)
            next
          
          hb <- hb - ob
          newNdx[curNdx:(curNdx+ob-1)] <- j
          curNdx <- curNdx + ob
        }
        #browser()
        # last index saved already as N

        return (newNdx)
}

####################################################
# Plotting Commands
# Plot the 95% CI over time: right now hard-coded to show beta/gamma only
# trueX
# shows
plot.ci <- function(saved.stats,trueX,trueTheta,sir.plotCI=0,col="blue")
{
   par(mfcol=c(2,2),mar=c(4,4,2,1), oma=c(0,0,0,1))
   XLab <- c("S", "I", "R", "D")
   TLab <- c("S->I", "I->R", "S->R","I->D")
   for (jj in c(1,2))
   {
       len <- dim(saved.stats)[1]
       plot(1:len,saved.stats[,(6*jj-5)],col=col,lwd=2,type="l",ylab=XLab[jj],
        xlab='',ylim=c(min(saved.stats[,(6*jj-4)]),max(saved.stats[,(6*jj-3)])))
       lines(1:len,saved.stats[,(6*jj-4)],col=col,lty=3)
       lines(1:len,saved.stats[,(6*jj-3)],col=col,lty=3)

       lines(1:len,trueX[,jj],col ='red',xlim=c(0,len),lwd=1.5)
       plot(1:len,saved.stats[,(6*jj-2)],col=col,xlab='time',ylab=TLab[jj],type="l",
       ylim=c(trueTheta[jj]*0.75,trueTheta[jj]*1.25), lwd=2 )
       abline(h=trueTheta[jj],col='red',lwd=1.5)
       if (sir.plotCI == 1) {
         lines(1:len,saved.stats[,(6*jj-1)],col=col,lty=3)
         lines(1:len,saved.stats[,6*jj],col=col,lty=3)
         
       }
   }
}

    
plot.ci.mcError <- function(in.stats, which.var,col1="#eeeeeeaa", col2="blue")
#  plots the Monte Carlo variability around the mean/95-quantiles of the posterior for the 
#  specified parameter which.var. The shaded region is the min/max across the MC runs in in.stats
#  col1: color of the shaded region; col2: color of the mean line
#  
#  First, need to initialize via: plot(c(1,len),y.range,type="n",xaxs="i")
{
  ii <-1
  #for (ii in 1:3) {
    ndx <- which.var*6-3+ii  
    
    len <- dim(in.stats)[2]
    plm <- rbind(apply(in.stats[,,ndx],2, min),apply(in.stats[,,ndx],2,max),apply(in.stats[,,ndx],2,mean))

  #}
    ndx1 <- which.var*6-3+2
    ndx2 <- which.var*6-3+3  
    
    plm2 <- rbind(apply(in.stats[,,ndx1],2, min),apply(in.stats[,,ndx2],2,max),apply(in.stats[,,ndx1],2,mean),apply(in.stats[,,ndx2],2,mean))

    #polygon(c(1:len,len:1),c(plm2[1,1:len],plm2[2,len:1]),col=col1,border=NA)
    lines(1:len,plm2[1,],col=col2,lwd=0.5,lty=2)
    lines(1:len,plm2[2,],col=col2,lwd=0.5,lty=2)
    
    lines(1:len,plm2[3,],col=col2,lwd=2,lty=2)
    lines(1:len,plm2[4,],col=col2,lwd=2,lty=2)
    polygon(c(1:len,len:1),c(plm[1,1:len],plm[2,len:1]),col=col1,border=NA)
    lines(1:len,plm[3,],col=col2,lwd=3)
  
}  


##########################################################################
#Use particles to filter for SIR model with discrete-time sampled 
#transitions 
#
#  model.params are the initialization parameters -- see the default case base.params above
#  LOOPN is the number of repetitions to do
#  aLW is the Liu-West weight alpha; if aLW > 1 then particle learning is done
#  verbose=CI/HIST will output some summary plots using the given color
#  if trueX/Y is provided, that would be used as the true scenario
#  model.propagate.func is the method to use for particle mutation/true system evolution: tauLeap/gillespie
#  resample.func is the resampling method: multinomial.resample/branchMinVar
#
plSIR <- function(N, T, dt=1, model.params=base.params, LOOPN=1,verbose="CI",
                col="blue",trueX=NULL,Y=NULL,
                model.propagate.func=tauLeap,resample.func=multinomial.resample)
{
    initP <- model.params$initP
    initX <- model.params$initX
    hyperPrior <- model.params$hyperPrior
    trueTheta <- model.params$trueTheta


    saved.stats <- array(0, dim=c(LOOPN,ceil(T/dt)+1,24))

    if (is.null(trueX)) # generate a scenario
    {
       trueX <- array(0, dim=c(ceil(T/dt)+1,4))
       Y <- array(0,dim=c(ceil(T/dt)+1,4))
       
       trueX[1,] <- initX

       # Fix the true state trajectory for replicability
       #set.seed(21)
       # Construct the true process and the observations
       for (i in 1:(T/dt)) 
       {
           out <- model.propagate.func(t(t(trueX[i,])),t(t(trueTheta)),t(initP))
           
           Y[i+1,] <- out$Y
           trueX[i+1,1:dim(out$X)[1]] <- out$X
       }
       #browser()
    }
    
    for (loop in 1:LOOPN) {
    
       # Initialize particles
       X <- t(array(rep(initX, N),dim=c(N.RXNS,N)))
       lambda <- t(array(rep(initP, N), dim=c(N.RXNS,N)))
       Suff <- t(array(rep(hyperPrior,N), dim=c(2*N.RXNS,N))) # sufficient conjugate statistics
       theta <- array(0,dim=c(N,N.RXNS))
       for (jj in 1:N.RXNS)
               theta[,jj] <- rgamma(N,Suff[,2*jj-1],Suff[,2*jj])
       
       p.weights <- array(1, N)/N

       curt <- 0 
       i <- 1
       totalWeight <- 1
       saved.stats[loop,i,] <- saveStats(X, theta)

       while (curt < T-dt & trueX[i,2] >0) {
           # Sample from the posterior Gamma mixtures
           for (jj in 1:N.RXNS)
               theta[,jj] <- rgamma(N,Suff[,2*jj-1],Suff[,2*jj])
               
           # resample    
           p.weights <- predictiveLikelihood(X, Y[i+1,], theta, lambda, p.weights)
           newIndex <- resample.func(p.weights) 
         
           X <- X[newIndex,] 
           Suff <- Suff[newIndex,]
           lambda <- lambda[newIndex,]
           theta <- theta[newIndex,]
         
           out <- model.propagate.func(t(X), t(theta),lambda,curY=Y[i+1,],hyper=Suff)
           X <- t(out$X); dX <- t(out$dX); Suff <- out$hyper
           p.weights <- rep(1/N, N)

           curt <- (i-1)*dt

           if ( i %%15  == 6 & verbose =="HIST") # plot posterior of infectiousness parameter
           {
               # construct exact Gamma pdf on a grid
               gridx <- seq(0.45,1,by=0.005)
               gridy <- array(0, length(gridx))
               for (jj in 1:length(gridx))
                  gridy[jj] <- sum(dgamma(gridx[jj],Suff[,1],Suff[,2]))

               plot(gridx,gridy/N,type="l",col=col,main=sprintf("t=%d",curt),xlab='S->I Rate')
             
               abline(v=trueTheta[1], col="red")
           }
           totalWeight <- totalWeight*sum(p.weights)

        i <- i+1
        saved.stats[loop,i,] <- saveStats(X,theta)
        # end of main loop
      }
   }
   
   if (verbose=='CI')
      plot.ci(saved.stats[1,1:i,],trueX[1:i,],trueTheta,1,col)
 
   return( list(stat=saved.stats,trueX=trueX,Y=Y))
}

