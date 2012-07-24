####################################################
# Plotting Commands related to the particle learning algorithms
# Plot the 95% CI over time: default is to show beta/gamma only
# which.dim -- what to plot
#
plot.ci <- function(saved.stats,trueX,trueTheta,sir.plotCI=0,col="blue",plot.diff=0, which.dim=c(1,2),in.legend=NULL,ltype=1)
{
   par(mfcol=c(2,length(which.dim)),mar=c(4,3,2,1)+0.1, oma=c(0,0,0,1))
   XLab <- c("S", "I", "R", "D")
   TLab <- c("S->I", "I->R", "S->R","I->D")
   n.comps <- dim(saved.stats)[1]
   col.list <- c(col,2:n.comps)
   
   for (jj in which.dim)
   {
      len <- dim(trueX)[1]
      if (plot.diff==1) {  # absolute difference
         for (k in 1:n.comps) {
            if (k==1)  
               plot(1:len,(saved.stats[k,,(6*jj-5)]-trueX[,jj]),col=col.list[k],lwd=2,type="l",
                   ylab=XLab[jj], xlab='',lty=ltype[k])
            else 
               lines(1:len,(saved.stats[k,,(6*jj-5)]-trueX[,jj]),col=col.list[k],lwd=2,lty=ltype[k])
            
            lines(1:len,(saved.stats[k,,(6*jj-4)]-trueX[,jj]),col=col.list[k],lty=3)
            lines(1:len,(saved.stats[k,,(6*jj-3)]-trueX[,jj]),col=col.list[k],lty=3)
          
         }
      }
      else if (plot.diff==2) {  # percent-difference
         for (k in 1:n.comps) {
            if (k==1) 
               plot(1:len,log(saved.stats[k,,(6*jj-5)]/trueX[,jj]),col=col.list[k],lwd=2,type="l",
                   ylab=XLab[jj], xlab='',lty=ltype[k])
            else
               lines(1:len,log(saved.stats[k,,(6*jj-5)]/trueX[,jj]),col=col.list[k],lwd=2,lty=ltype[k])
      
            lines(1:len,log(saved.stats[k,,(6*jj-4)]/trueX[,jj]),col=col.list[k],lty=3)
            lines(1:len,log(saved.stats[k,,(6*jj-3)]/trueX[,jj]),col=col.list[k],lty=3)
          
         }
      }     
      else { # plot absolute values
         for (k in 1:n.comps) {
            if (k==1)
               plot(1:len,saved.stats[k,,(6*jj-5)],col=col.list[k],lwd=2,type="l",ylab='',main=XLab[jj],lty=ltype[k],
                  xlab='',ylim=c(min(saved.stats[k,,(6*jj-4)]),max(max(saved.stats[k,,(6*jj-5)]),0.7*max(saved.stats[k,,(6*jj-3)]))))
            else
               lines(1:len,saved.stats[k,,(6*jj-5)],col=col.list[k],lwd=2,lty=ltype[k])
            
            lines(1:len,saved.stats[k,,(6*jj-4)],col=col.list[k],lty=3)
            lines(1:len,saved.stats[k,,(6*jj-3)],col=col.list[k],lty=3)

            if (!is.null(trueX))
               lines(1:len,trueX[,jj],col ='red',xlim=c(0,len),lwd=1.5)
            if (n.comps == 1)
              legend("topright", c("Truth", "Median", "95% CI"), col=c("red",col,col), lty=c(ltype[1],ltype[1],3), lwd=c(1.5,2,1))
          }
      }
      if (n.comps > 1 & !is.null(in.legend))
         legend("topright", in.legend, lty=ltype, col=col.list, lwd=2)

      for (k in 1:n.comps) {
         if (k==1) {
            plot(1:len,saved.stats[k,,(6*jj-2)],col=col.list[k],xlab='Period',ylab='',main=TLab[jj],type="l",lty=ltype[k],
               ylim=c(min(saved.stats[k,,(6*jj-1)])*1.02, max(saved.stats[k,,6*jj]))*0.98 , lwd=2.5)
      #c(trueTheta[jj]*0.75,trueTheta[jj]*1.25)
            abline(h=trueTheta[jj],col='red',lwd=1.5)  
         }
         else 
            lines(1:len,saved.stats[k,,(6*jj-2)],col=col.list[k], lwd=2.5,lty=ltype[k])  
          
         if (sir.plotCI == 1) {
            lines(1:len,saved.stats[k,,(6*jj-1)],col=col.list[k],lty=3)
            lines(1:len,saved.stats[k,,6*jj],col=col.list[k],lty=3)
         }
      }   

   }
}

plot.ci.mcError <- function(in.stats, which.var,col1="#ddddddaa", col2="blue")
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

    plm2 <- rbind(apply(in.stats[,,ndx1],2,min),apply(in.stats[,,ndx2],2,max),apply(in.stats[,,ndx1],2,mean),apply(in.stats[,,ndx2],2,mean))

    #polygon(c(1:len,len:1),c(plm2[1,1:len],plm2[2,len:1]),col=col1,border=NA)
    lines(1:len,plm2[1,],col=col2,lwd=0.5,lty=2)
    lines(1:len,plm2[2,],col=col2,lwd=0.5,lty=2)

    lines(1:len,plm2[3,],col=col2,lwd=2,lty=2)
    lines(1:len,plm2[4,],col=col2,lwd=2,lty=2)
    polygon(c(1:len,len:1),c(plm[1,1:len],plm[2,len:1]),col=col1,border=NA)
    lines(1:len,plm[3,],col=col2,lwd=3)

}

