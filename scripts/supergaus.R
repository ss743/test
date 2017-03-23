library("minpack.lm")
source("functions.R")

supergausfit <- function(input,bereich,weighted=FALSE,mu=0,sig0=0,N0=0,N01=2000000,sig01=80,mu01=750){ #--- Fitten der Exponentialfunktion
  
  thegaussian <- y ~ C + N/(sqrt(2*pi)*sig)*exp(-(x-mu)^2/(2*sig^2))+ N1/(sqrt(2*pi)*sig1)*exp(-(x-mu1)^2/(2*sig1^2))
  
  daten=subset(input,x>=bereich[1] & x<= bereich[2])
  ymin=min(daten$y)
  if(N0==0){
    ymax=max(daten$y)
  } else {
    ymax=N0
  }
  if(mu==0){
    mu0 =daten$x[which.max(daten$y)]
  } else {
    mu0=mu
  }
  if(sig0==0)
  {
    #sig0=(daten$x[bereich[2]]-daten$x[bereich[1]])/3
    sig0=(bereich[2]-bereich[1])/6
  }
  N0=(sqrt(2*pi)*sig0)*(ymax-ymin)
#  cat(paste("\nC=",ymin,"\nN=",N0,"\nmu=",mu0,"\nsigma=",sig0,sep=""))
#  cat(paste("\nN1=",N0,"\nmu=",mu0,"\nsigma=",sig0,sep=""))
  #plot (function(x){ymin + N0/(sqrt(2*pi)*sig0)*exp(-(x-mu0)^2/(2*sig0^2)) + N01/(sqrt(2*pi)*sig01)*exp(-(x-mu01)^2/(2*sig01^2))},bereich[1],bereich[2],add=TRUE,col="green")
  if(weighted)
  {
    err=daten$sy
    fit = nlsLM(thegaussian,daten,weights=1/err^2,start=list(C=ymin,N=N0,mu=mu0,sig=sig0, N1=N0, mu1=mu0, sig1=sig0))
  }
  else
    fit = nlsLM(thegaussian,daten,start=list(C=ymin,N=ymax,mu=mu0,sig=sig0, N1=N0, mu1=mu0, sig1=sig0))
  
  
  #chiquadratndf=sum(residuals(fit)^2/abs(fitted(fit)))/summary(fit)$df[2]
  chiquadratndf=sum(summary(fit)[[2]]^2)/summary(fit)[[4]][[2]]
  fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
  return(fitdata)  
}

plotsupergaus <- function(fitdata,bereich,col="red"){ #--- Plotten der gefitteten Gaußfunktion in vorhandenen Graph
  
  N<-fitdata["N","Estimate"]
  C<-fitdata["C","Estimate"]
  mu<-fitdata["mu","Estimate"]
  sig<-fitdata["sig","Estimate"]
  
  N1<-fitdata["N1","Estimate"]
  mu1<-fitdata["mu1","Estimate"]
  sig1<-fitdata["sig1","Estimate"]
  
  
  plot (function(x){C + N/(sqrt(2*pi)*sig)*exp(-(x-mu)^2/(2*sig^2)) + N1/(sqrt(2*pi)*sig1)*exp(-(x-mu1)^2/(2*sig1^2))},bereich[1],bereich[2],add=TRUE,col=col)
  
}

printsupergausdata <- function(fitdata,title=""){ #--- Ausgabe der Gaußfit-Daten
  
  mu<-fitdata["mu","Estimate"]
  smu<-fitdata["mu","Std. Error"]
  mu1<-fitdata["mu1","Estimate"]
  smu1<-fitdata["mu1","Std. Error"]
  A<-fitdata["N","Estimate"]
  sA<-fitdata["N","Std. Error"]
  C<-fitdata["C","Estimate"]
  sC<-fitdata["C","Std. Error"]
  sig1<-fitdata["sig1","Estimate"]
  ssig1<-fitdata["sig1","Std. Error"]
  sig<-fitdata["sig","Estimate"]
  ssig<-fitdata["sig","Std. Error"]
  A1<-fitdata["N1","Estimate"]
  sA1<-fitdata["N1","Std. Error"]
  
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

  cat(" \\mu1 = ")
  cat(mu1)
  cat("\\pm")
  cat(smu1)
  cat("\\\\\n")
  
  cat(" \\sigma1 = ")
  cat(sig1)
  cat("\\pm")
  cat(ssig1)
  cat("\\\\\n")
  
  cat(" A1 = ")
  cat(A1)
  cat("\\pm")
  cat(sA1)
  cat("\\\\\n")
  
  cat(" C = ")
  cat(C)
  cat("\\pm")
  cat(sC)
  cat("\\\\\n")
  
}
