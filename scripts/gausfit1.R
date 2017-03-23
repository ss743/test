library("minpack.lm")
source("functions.R")

gausfit <- function(input,bereich,weighted=FALSE,sig0=0,N0=0){ #--- Fitten der Exponentialfunktion
  
  thegaussian <- y ~ C + N/(sqrt(2*pi)*sig)*exp(-(x-mu)^2/(2*sig^2))
  
  daten=subset(input,x>=bereich[1] & x<= bereich[2])
  ymin=min(daten$y)
  if(N0==0){
    ymax=max(daten$y)
  } else {
    ymax=N0
  }
  mu0 =daten$x[which.max(daten$y)]
  if(sig0==0)
  {
    #sig0=(daten$x[bereich[2]]-daten$x[bereich[1]])/3
    sig0=(bereich[2]-bereich[1])/6
  }
  N0=(sqrt(2*pi)*sig0)*(ymax-ymin)
  #cat(paste("\nC=",ymin,"\nN=",N0,"\nmu=",mu0,"\nsigma=",sig0,sep=""))
  #plot (function(x){ymin + N0/(sqrt(2*pi)*sig0)*exp(-(x-mu0)^2/(2*sig0^2))},bereich[1],bereich[2],add=TRUE,col="yellow")
  if(weighted)
    err=daten$sy
  else
    err=1
  
  fit = nlsLM(thegaussian,daten,weights=1/err^2,start=list(C=ymin,N=ymax,mu=mu0,sig=sig0))
  
  chiquadratndf=sum(summary(fit)[[2]]^2)/summary(fit)[[4]][[2]]
  fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
  
  return(fitdata)
  
}
sixgaus <- function(input,bereich,weighted=FALSE,mu01,mu02,mu03,mu04,mu05,mu06){
  thegaussian <- y ~ C + N1*exp(-(x-mu1)^2/(2*sig1^2)) + N2*exp(-(x-mu2)^2/(2*sig2^2)) + N3*exp(-(x-mu3)^2/(2*sig3^2)) + N4*exp(-(x-mu4)^2/(2*sig4^2)) + N5*exp(-(x-mu5)^2/(2*sig5^2)) + N6*exp(-(x-mu6)^2/(2*sig6^2))
  daten=subset(input,x>=bereich[1] & x<= bereich[2])
  ymin=min(daten$y)
  ymax=max(daten$y)
  sig0=(bereich[2]-bereich[1])/36
  N0=(sqrt(2*pi)*sig0)*(ymax-ymin)
  
  if(weighted)
    err=daten$sy
  else
    err=daten$y/daten$y
  
  fit = nlsLM(thegaussian,daten,weights=1/err^2,start=list(C=ymin,N1=ymax,N2=ymax,N3=ymax,N4=ymax,N5=ymax,N6=ymax,mu1=mu01,mu2=mu02,mu3=mu03,mu4=mu04,mu5=mu05,mu6=mu06,sig1=sig0,sig2=sig0,sig3=sig0,sig4=sig0,sig5=sig0,sig6=sig0))
  
  chiquadratndf=sum(summary(fit)[[2]]^2)/summary(fit)[[4]][[2]]
  fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
  
  return(fitdata)
  
}

plotgaus <- function(fitdata,bereich,col="red",lwd=1,lty=1){ #--- Plotten der gefitteten Gaußfunktion in vorhandenen Graph
  
  N<-fitdata["N","Estimate"]
  C<-fitdata["C","Estimate"]
  mu<-fitdata["mu","Estimate"]
  sig<-fitdata["sig","Estimate"]
  
  plot (function(x){C + N/(sqrt(2*pi)*sig)*exp(-(x-mu)^2/(2*sig^2))},bereich[1],bereich[2],add=TRUE,col=col,n=10000,lwd=lwd,lty=lty)
  
}

plotsix <- function(fitdata,bereich,col="red",lwd=1,lty=1){ #--- Plotten der gefitteten Gaußfunktion in vorhandenen Graph
  
  C<-fitdata["C","Estimate"]
  N1<-fitdata["N1","Estimate"]
  mu1<-fitdata["mu1","Estimate"]
  sig1<-fitdata["sig1","Estimate"]
  N2<-fitdata["N2","Estimate"]
  mu2<-fitdata["mu2","Estimate"]
  sig2<-fitdata["sig2","Estimate"]
  N3<-fitdata["N3","Estimate"]
  mu3<-fitdata["mu3","Estimate"]
  sig3<-fitdata["sig3","Estimate"]
  N4<-fitdata["N4","Estimate"]
  mu4<-fitdata["mu4","Estimate"]
  sig4<-fitdata["sig4","Estimate"]
  N5<-fitdata["N5","Estimate"]
  mu5<-fitdata["mu5","Estimate"]
  sig5<-fitdata["sig5","Estimate"]
  N6<-fitdata["N6","Estimate"]
  mu6<-fitdata["mu6","Estimate"]
  sig6<-fitdata["sig6","Estimate"]
  
  plot (function(x){C + N1*exp(-(x-mu1)^2/(2*sig1^2)) + N2*exp(-(x-mu2)^2/(2*sig2^2)) + N3*exp(-(x-mu3)^2/(2*sig3^2)) + N4*exp(-(x-mu4)^2/(2*sig4^2)) + N5*exp(-(x-mu5)^2/(2*sig5^2)) + N6*exp(-(x-mu6)^2/(2*sig6^2))},bereich[1],bereich[2],add=TRUE,col=col,n=10000,lwd=lwd,lty=lty)
  
}


printfitdata <- function(fitdata,title=""){ #--- Ausgabe der Gaußfit-Daten
  
  mu<-fitdata["mu","Estimate"]
  smu<-fitdata["mu","Std. Error"]
  A<-fitdata["N","Estimate"]
  sA<-fitdata["N","Std. Error"]
  sig<-fitdata["sig","Estimate"]
  ssig<-fitdata["sig","Std. Error"]
  
  cat(paste("\\text{",title,"}",sep=""))
  cat("\\\\\n")
  
  cat(" \\mu = ")
  cat(mu)
  cat("\\pm")
  cat(smu)
  cat("\\\\\n")
  
  cat(" \\sigma = ")
  cat(sig)
  cat("\\pm")
  cat(ssig)
  cat("\\\\\n")
  
  cat(" A = ")
  cat(A)
  cat("\\pm")
  cat(sA)
  cat("\\\\\n")
  

}
