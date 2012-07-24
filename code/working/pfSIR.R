

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

    if (is.null(Y)) # generate a scenario on the fly
    {
       scen <- generate.scenario(model.params,model.propagate.func,T)
       trueX <- scen$X
       Y <-scen$Y
    }
    
    for (loop in 1:LOOPN) {    # run the filter LOOPN times to understand MC variance if needed
    
       # Initialize particles
       X <- t(array(rep(initX,N),dim=c(N.RXNS,N)))
       #X[,2] <- rpois(N, initX[2])  # initial infecteds
       #X[,1] <- sum(initX) - X[,2]           # S_0 = N - I_0, R_0 =0, D_0 = 0
       pSamp <- t(array(rep(initP, N), dim=c(N.RXNS,N)))
       #pSamp[,1] <- rbeta(N, initP[1]*100/(1-initP[1]), 100)  # beta centered on initP[1]
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

       ######### MAIN LOOP OVER OBSERVATIONS ##############
       while (curt < T-dt) {
         if (aLW > 1 & is.null(Suff) == 0)  # Storvik filter: sample from the posterior Gamma mixture
           for (jj in 1:N.RXNS)
               theta[,jj] <- rgamma(N,Suff[,2*jj-1],Suff[,2*jj])
           
         else if (aLW < 1) # use the Liu-West thetas
            theta <- fixTheta
        
        # propagate particles 
        out <- model.propagate.func(t(X), t(theta),pSamp,curY=Y[i,],hyper=Suff)
        X <- t(out$X); dX <- t(out$dX); Suff <- out$hyper
        
        # update weights
        p.weights <- updateWeights(dX,Y[i,],pSamp,p.weights)
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
         pSamp <- pSamp[newIndex,]
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
   
   if (verbose != "NONE")
      kd <- build.density(X,theta,Suff)
    
   if (verbose=='CI')  # plot some CI over time 
      plot.ci(saved.stats,trueX[1:i,],trueTheta,1,col)
 
   return( list(stat=saved.stats,trueX=trueX,Y=Y,density=kd))
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
       
    if (is.null(curY) | is.nan(curY)) 
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
    hyper2 <- array(1, dim=c(N.RXNS,2,dim(X)[2])) # one.step.C takes in 4x2xN while hyper is Nx8
    
    out <- one.step.C(X, hyper2, theta, t(prop), sample=T)
    
    for (j in 1:dim(X)[2])
    {    
       #out <- simulate.one.step.C(X[,j],theta[,j],sum(X[,j]) )
       #dX[,j] <- out$dX
       #X[,j] <- out$X
       h[j,]  <- hazard.R(X[,j], sum(X[,j]))   
    }
    #browser()
    newX <- out$newX
    newX[is.nan(newX)] <- 0  # to take care of the case when some proportions are zero
       
    if (is.null(curY))
        Y <- rbinom(N.RXNS, newX[,1], prop)
    else
        Y <- curY
    #Y[is.nan(Y)] <- 0

    # update the hyperparameters
    if (!is.null(hyper) & !is.nan(Y[1])) 
      for (i in 1:N.RXNS) {
        hyper[,2*i-1] <- hyper[,2*i-1] + newX[i,] #Y[i]
        hyper[,2*i] <- hyper[,2*i] +  h[,i] # prop[,i]*
        
        hyper[,2*(i+N.RXNS)-1] <- hyper[,2*(i+N.RXNS)-1] + Y[i]
        hyper[,2*(i+N.RXNS)] <- hyper[,2*(i+N.RXNS)] + newX[i,] - Y[i]
      }
    
    return(list(X=out$X,dX=newX,hyper=hyper,Y=Y))
}

###############################################
# Generate a scenario of the SIRD model
# output is a path (X,Y)
generate.scenario <- function(model.params,model.propagate.func,T,seed=NULL)    
{

        X <- array(0, dim=c(T+1,N.RXNS))
        Y <- array(0,dim=c(T+1,N.RXNS))
       
         X[1,] <- model.params$initX
         theta <- model.params$trueTheta
         theta[1] <- theta[1]+0.1*rnorm(1);

       # Fix the true state trajectory for replicability
       if (!is.null(seed))
         set.seed(seed)
       # Construct the true process and the observations
       for (i in 1:T) 
       {
           out <- model.propagate.func(t(t(X[i,])),t(t(theta)),t(model.params$initP))
           
           Y[i+1,] <- out$Y
           X[i+1,1:dim(out$X)[1]] <- out$X
       }
       return(list(X=X,Y=Y,theta=theta))
}


#######################################################
# Update the weights of the particles using the binomial sampling
updateWeights <- function(dX, curY, prop, weights)
{
    for (jj in 1:N.RXNS)
      weights <- weights*dbinom(curY[jj], (dX[,jj]), prop[,jj])
      
    return(weights)

}

#######################################################
# Predictive likelihood of the next observation using Poisson approximations
predictiveLikelihood <- function(X, nextY, theta, prop, weights, Suff)
{
    h <- X  # just to set the size of h correctly
    
    for (j in 1:dim(X)[1])
       h[j,]  <- hazard.R(X[j,], sum(X[j,]))   
       
    for (jj in 1:N.RXNS) {
        weights <- weights*dpois(nextY[jj],prop[,jj]*theta[,jj]*h[,jj])
        #Analytic form for the predictive likelihood using Negative Binomial
        #pr <- Suff[,2*jj]/(prop[,jj]*h[,jj] + Suff[,2*jj])
        
        #ndx <- which (Suff[,2*jj-1] > 0 & pr > 0)
        #weights[ndx] <- weights[ndx]*dnbinom(nextY[jj],size=Suff[ndx,2*jj-1],prob=pr[ndx])
        #Use logs to make sure things do not blow up
        #GamTerm <- log(Suff[,2*jj])*Suff[,2*jj-1]- log(Suff[,2*jj]+prop[,jj]*h[,jj])*(nextY[jj]+Suff[,2*jj-1])
        #if (nextY[jj] > 0)
        #   GamTerm <- GamTerm + lgamma( nextY[jj] + Suff[,2*jj-1]) - lgamma(Suff[,2*jj-1]) - lgamma(nextY[jj]+1) + log(prop[,jj]*h[,jj])*nextY[jj]
        #weights <- weights*exp(GamTerm)  
    }
    #browser()
    return(weights)

}

#######################################################
# Summary statistics: for each X-coordinate and theta parameter
# save the 95% CI and the median generated by the particle cloud
saveStats <- function(X, theta)
{
    summ.stat <- array(0,dim=c(6*dim(X)[2],1))
    for (jj in 1:dim(X)[2])
    {
      #summ.stat[(jj-1)*6+1] <- mean(X[,jj])
      summ.stat[((jj-1)*6+1):((jj-1)*6+3)] <- quantile(X[,jj],c(0.5,0.025,0.975))

      #summ.stat[(jj-1)*6+4] <- mean(theta[,jj])
      summ.stat[((jj-1)*6+4):((jj-1)*6+6)] <- quantile(theta[,jj],c(0.5,0.025,0.975))
    }
    return(summ.stat)
}

#######################################################
# Construct the densities of the particles
build.density <- function(X, theta, Suff,sampP=NULL)
{
   N <- dim(X)[1]
   kd <- list()
   h1 <- hist(X[,1],breaks=seq(min(X[,1]),max(X[,1])+5, by=5),plot=F) # for S 
   kd[[1]] <- stepfun(h1$breaks, c(0,h1$counts/N/5,0)) # S
   
   h1 <- hist(X[,2],breaks=min(X[,2]):(max(X[,2])+1),plot=F) # for I 
   kd[[2]] <- stepfun(h1$breaks, c(0,h1$counts/N,0)) #(density(X[,2],bw=1) # I
   
   kd[[3]] = density(theta[,1]) # S->I
   kd[[4]] = density(theta[,2]) # I->R
   
   if (is.null(Suff) == F) {
     h1 <- array(0, length(kd[[3]]$x))
     for (jj in 1:length(kd[[3]]$x))
        h1[jj] <- sum(dgamma(kd[[3]]$x[jj],Suff[,1],Suff[,2]))
     kd[[3]]$y <- h1/N
   
     h1 <- array(0, length(kd[[4]]$x))
     for (jj in 1:length(kd[[4]]$x))
        h1[jj] <- sum(dgamma(kd[[4]]$x[jj],Suff[,3],Suff[,4]))
     kd[[4]]$y <- h1/N
   }
   
   if (is.null(sampP) == F) {
     kd[[5]] <- density(sampP[,1])
     kd[[6]] <- density(sampP[,2])
   }  
   
   return(kd)
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

    if (is.null(Y)) # generate a scenario
    {
        scen <- generate.scenario(model.params,model.propagate.func,T)
        Y <- scen$Y
        trueX <- scen$X 
        trueTheta <- scen$theta
    }     
    
    for (loop in 1:LOOPN) {
    
       # Initialize particles
       X <- t(array(rep(initX, N),dim=c(N.RXNS,N)))
       #X[,2] <- rpois(N, initX[2])  # initial infecteds
       #X[,1] <- sum(initX) - X[,2]           # S_0 = N - I_0, R_0 =0, D_0 = 0
       pSamp <- t(array(rep(initP, N), dim=c(N.RXNS,N)))
       #pSamp[,1] <- rbeta(N, initP[1]*100/(1-initP[1]), 100)  # beta centered on initP[1]
       Suff <- t(array(rep(hyperPrior,N), dim=c(4*N.RXNS,N))) # sufficient conjugate statistics
       theta <- array(0,dim=c(N,N.RXNS))
       for (jj in 1:N.RXNS)
               theta[,jj] <- rgamma(N,Suff[,2*jj-1],Suff[,2*jj])
       for (jj in 1:N.RXNS)
               pSamp[,jj] <- rbeta(N,Suff[,2*(jj+N.RXNS)-1],Suff[,2*(jj+N.RXNS)])
       
       p.weights <- array(1, N)/N

       curt <- 0 
       i <- 1
       totalWeight <- 1
       saved.stats[loop,i,] <- saveStats(X, theta)

       ####### MAIN LOOP OVER OBSERVATIONS
       while (curt < T-dt) {
           #browser()
           
           # Sample from the posterior mixtures
           for (jj in 1:N.RXNS)
               theta[,jj] <- rgamma(N,Suff[,2*jj-1],Suff[,2*jj])
           for (jj in 1:N.RXNS)
               pSamp[,jj] <- rbeta(N,Suff[,2*(jj+N.RXNS)-1],Suff[,2*(jj+N.RXNS)])
            
               
           # resample 
           if (!is.nan(Y[i,1])) # else missing observation for that date
           {  
             p.weights <- predictiveLikelihood(X, Y[i,], theta, pSamp, p.weights,Suff)
             newIndex <- resample.func(p.weights) 
         
             X <- X[newIndex,] 
             Suff <- Suff[newIndex,]
             pSamp <- pSamp[newIndex,]
             theta <- theta[newIndex,]
             p.weights <- rep(1/N, N)
           }
           
           out <- model.propagate.func(t(X), t(theta),pSamp,curY=Y[i,],hyper=Suff)
           X <- t(out$X); Suff <- out$hyper

           curt <- (i-1)*dt

           if ( i %%15  == 6 & verbose =="HIST") # plot posterior of infectiousness parameter
           {
               # construct exact Gamma pdf on a grid
               gridx <- seq(initP[1]-0.1,initP[1]+0.1,by=0.001)
               gridy <- array(0, length(gridx))
               for (jj in 1:length(gridx))
                  gridy[jj] <- sum(dbeta(gridx[jj],Suff[,9],Suff[,10]))
               
               #gridx <- seq(trueTheta[1]-0.25,trueTheta[1]+0.25,by=0.005)
               #gridy <- array(0, length(gridx))
               #for (jj in 1:length(gridx))
               #   gridy[jj] <- sum(dgamma(gridx[jj],Suff[,1],Suff[,2]))
               #browser()

               par(mfg=c(1,ceiling(i/15))) 
               lines(gridx,gridy/N,type="l",col=col,main=sprintf("t=%d",curt),xlab='S->I Rate',yaxt="n")
             
               abline(v=initP[1], col="red")
           }
           totalWeight <- totalWeight*sum(p.weights)

        i <- i+1
        saved.stats[loop,i,] <- saveStats(X,theta)
        # end of main loop
      }
   }
  
   if (verbose != "NONE")
     kd <- build.density(X,theta,Suff,pSamp)
   
   
   #browser()
   
   # high-quality quantile computation
   for (jj in 1:N.RXNS) {
      theta <- rgamma(max(N,25000),Suff[,2*jj-1],Suff[,2*jj])
      saved.stats[loop, i, ((jj-1)*6+4):((jj-1)*6+6)] <- quantile(theta,c(0.5,0.025,0.975))
   }      

   if (verbose=='CI') 
      plot.ci(saved.stats,trueX[1:i,],trueTheta,1,col)
      
 
   return( list(stat=saved.stats,trueX=trueX,Y=Y,theta=trueTheta,density=kd))
}

