# - - - - - - - - - - - - - - - - - 
# Analysis of tracing
# - - - - - - - - - - - - - - - - - 

library(tidyverse)

# setwd("~/Documents/GitHub/2020-ncov-forward-back-tracing/")


# Define functions --------------------------------------------------------

tertiary_cases <- function(R=1.2, # reproduction number
                            k=0.3, # dispersion
                            q=1, # probability cluster identified
                            p=0.5, # probability contact traced
                            c1=0.5, # proportional reduction from tracing index case
                            c2=0.5, # proportional reduction from tracing contact of cluster 
                            d=0.5 # probability detection of cases
                           ){
  tertiary_avert <- function(R, k, q, p, c1, c2, d){
      R^2*p*c1 + (1-(1-q*p)*(1-d))*p*c2*R^3*(1+1/k)
  }
  
  # No tracing
  baseline_tertiary <- tertiary_avert(R,k,1,1,1,1,1)
  
  # Cases averted with forward tracing
  forward_tertiary_avert <- tertiary_avert(R,k,0,p,c1,c2,d)

  # Cases averted with forward + backward tracing
  forward_backward_tertiary_avert <- tertiary_avert(R,k,q,p,c1,c2,d)

  list(base=baseline_tertiary,
    forward_avert=forward_tertiary_avert,
    fwd_back_avert=forward_backward_tertiary_avert,
    backward_gain=forward_backward_tertiary_avert-forward_tertiary_avert)
}


#secondary_cases()

# Plot functions ----------------------------------------------------------

#kk <- seq(0.1,0.5,0.1)

plotbypc <- function(R=1.2, k=0.3, q=1, d = 0.5, ps=seq(0.0,1,0.1), cs=seq(0.2,1,0.2), noplot=F){
    psmat <- matrix(ps,length(ps),length(cs))
    csmat <- matrix(cs,length(ps),length(cs),byrow=T)
    tertiary <- tertiary_cases(R,k,q,psmat,csmat,csmat,d)
    if(noplot){return(tertiary)} # skip plotting
    
    cols1=terrain.colors(length(cs)+2)[1:length(cs)]
    cols2=terrain.colors(length(cs)+2)[1:length(cs)]
    cols3=terrain.colors(length(cs)+2)[1:length(cs)]
    
    matplot(tertiary$forward_avert,x=ps,type="o",pch=19,lty=1,col=cols1,lwd=2,ylab="cases averted (forward tracing)",xlab="success rate of contact tracing",ylim=c(0,tertiary$base*1.02))
    legend("topleft",legend=rev(paste("adherence = ",cs)),pch=19,col=rev(cols1),bty="n")
    abline(h=tertiary$base,lwd=2,lty=2)
    matplot(tertiary$fwd_back_avert,x=ps,type="o",pch=19,lty=1,col=cols2,lwd=2,ylab="cases averted (forward + backward tracing)",xlab="success rate of contact tracing",ylim=c(0,tertiary$base*1.02))
    legend("topleft",legend=rev(paste("adherence = ",cs)),pch=19,col=rev(cols1),bty="n")
    abline(h=tertiary$base,lwd=2,lty=2)
    matplot(tertiary$backward_gain,x=ps,type="o",pch=19,lty=1,col=cols3,lwd=2,ylab="difference",xlab="success rate of contact tracing",ylim=c(0,tertiary$base*1.02))
    legend("topleft",legend=rev(paste("adherence = ",cs)),pch=19,col=rev(cols1),bty="n")
    abline(h=tertiary$base,lwd=2,lty=2)
    return(tertiary)
}
