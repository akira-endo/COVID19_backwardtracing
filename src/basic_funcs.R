# - - - - - - - - - - - - - - - - - 
# Analysis of tracing
# - - - - - - - - - - - - - - - - - 

# Define functions --------------------------------------------------------

tertiary_cases <- function(R=1.2, # reproduction number
                            k=0.3, # dispersion
                            b=1, # probability cluster identified
                            q=0.5, # probability contact traced
                            c1=0.5, # proportional reduction from tracing index case
                            c2=0.5, # proportional reduction from tracing contact of cluster 
                            d=0.5 # probability detection of cases
                           ){
  tertiary_avert <- function(R, k, b, q, c1, c2, d){
      R^2*q*c1 + (1-(1-b*q)*(1-d))*q*c2*R^3*(1+1/k)
  }
  
  # No tracing
  baseline_tertiary <- tertiary_avert(R,k,1,1,1,1,1)
  
  # Cases averted with forward tracing
  forward_tertiary_avert <- tertiary_avert(R,k,0,q,c1,c2,d)

  # Cases averted with forward + backward tracing
  forward_backward_tertiary_avert <- tertiary_avert(R,k,b,q,c1,c2,d)

  list(base=baseline_tertiary,
    forward_avert=forward_tertiary_avert,
    fwd_back_avert=forward_backward_tertiary_avert,
    backward_gain=forward_backward_tertiary_avert-forward_tertiary_avert)
}


# Plot functions ----------------------------------------------------------


plotbypc <- function(R=1.2, k=0.3, b=1, d = 0.5, qs=seq(0.0,1,0.1), cs=seq(0.2,1,0.2), noplot=F, panel=rep("",3)){
    qsmat <- matrix(qs,length(qs),length(cs))
    csmat <- matrix(cs,length(qs),length(cs),byrow=T)
    tertiary <- tertiary_cases(R,k,b,qsmat,csmat,csmat,d)
    if(noplot){return(tertiary)} # skip plotting
    
    cols1=terrain.colors(length(cs)+2)[1:length(cs)]
    cols2=terrain.colors(length(cs)+2)[1:length(cs)]
    cols3=terrain.colors(length(cs)+2)[1:length(cs)]
    
    pars = paste0(" R = ",R,", k = ",k,", b = ",format(b,nsmall=1),", d = ",d,paste0(rep(" ",10),collapse=""))
     matplot(tertiary$forward_avert/tertiary$base,x=qs,type="o",pch=19,lty=1,col=cols1,lwd=2,ylab="effectiveness (forward)",xlab="proportion of contacts traced (q)",ylim=c(0,1.02),main=c(pars,paste(c(panel[1],rep(" ",45)),collapse="")))
    legend("topleft",legend=rev(paste(c(character(4),"reduction if traced (c) ="),cs)),pch=19,col=rev(cols1),bty="n")
    matplot(tertiary$fwd_back_avert/tertiary$base,x=qs,type="o",pch=19,lty=1,col=cols2,lwd=2,ylab="effectiveness (forward + backward)",xlab="proportion of contacts traced (q)",ylim=c(0,1.02),main=c("",paste(c(panel[2],rep(" ",45)),collapse="")))
    legend("topleft",legend=rev(paste(c(character(4),"reduction if traced (c) ="),cs)),pch=19,col=rev(cols1),bty="n")
    matplot(tertiary$backward_gain/tertiary$base,x=qs,type="o",pch=19,lty=1,col=cols3,lwd=2,ylab="increment with backward tracing",xlab="proportion of contacts traced (q)",ylim=c(0,1.02),main=c("", paste(c(panel[3],rep(" ",45)),collapse="")))
    legend("topleft",legend=rev(paste(c(character(4),"reduction if traced (c) ="),cs)),pch=19,col=rev(cols1),bty="n")
    return(tertiary)
}

plotbyRk <- function(Rs = c(1.2, 1.5, 2, 2.5), ks=seq(0.1,0.5,0.1), b = 1, d= 0.5, q = 0.8, c = 0.8, noplot = F, panel = rep("",3)){
    ksmat <- matrix(ks,length(ks),length(Rs))
    Rsmat <- matrix(Rs,length(ks),length(Rs),byrow=T)
    tertiary <- tertiary_cases(Rsmat,ksmat,b,q,c,c,d)
    if(noplot){return(tertiary)} # skip plotting
    ymax=max(tertiary$fwd_back_avert)#max(tertiary$base)
    
    cols1=topo.colors(length(Rs)+10)[c(1,3+1:(length(Rs)-1))]
    cols2=topo.colors(length(Rs)+10)[c(1,3+1:(length(Rs)-1))]
    cols3=topo.colors(length(Rs)+10)[c(1,3+1:(length(Rs)-1))]
    
    pars = paste0(" q = ",q,", c = ",c,", b = ",format(b,nsmall=1),", d = ",d,paste0(rep(" ",10),collapse=""))
     matplot(tertiary$forward_avert,x=ks,type="o",pch=19,lty=1,col=rev(cols1),lwd=2,ylab="cases averted (forward)",xlab="overdispersion parameter (k)",ylim=c(0,ymax*1.02),main=c(pars,paste(c(panel[1],rep(" ",45)),collapse="")))
    legend("topright",legend=rev(paste(c(character(length(Rs)-1),"R ="),format(Rs,nsmall=1))),pch=19,col=(cols1),bty="n")
    #abline(h=tertiary$base,lwd=2,lty=2)
    matplot(tertiary$fwd_back_avert,x=ks,type="o",pch=19,lty=1,col=rev(cols2),lwd=2,ylab="cases averted (forward + backward)",xlab="overdispersion parameter (k)",ylim=c(0,ymax*1.02),main=c("",paste(c(panel[2],rep(" ",45)),collapse="")))
    legend("topright",legend=rev(paste(c(character(length(Rs)-1),"R ="),format(Rs,nsmall=1))),pch=19,col=(cols2),bty="n")
    #abline(h=tertiary$base,lwd=2,lty=2)
    matplot(tertiary$backward_gain,x=ks,type="o",pch=19,lty=1,col=rev(cols3),lwd=2,ylab="increment with backward tracing",xlab="overdispersion parameter (k)",ylim=c(0,ymax*1.02),main=c("", paste(c(panel[3],rep(" ",45)),collapse="")))
    legend("topright",legend=rev(paste(c(character(length(Rs)-1),"R ="),format(Rs,nsmall=1))),pch=19,col=(cols2),bty="n")
    #abline(h=tertiary$base,lwd=2,lty=2)
    return(tertiary)
}