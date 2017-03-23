source("functions.R")

gausfit <- function(input,bereich=c(1,length(input)),weighted=FALSE,sig0=0,N0=0){ #--- Fitten der Gaußfunktion
  
  thegaussian <- y ~ C + N*exp(-(x-mu)^2/(2*sig^2))
  
  daten=input[bereich[1]:bereich[2],]
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
    sig0=(bereich[2]-bereich[1])/3
  }
  err=daten$sy

  #plot(function(x){ymin + ymax*exp(-(x-mu0)^2/(2*sig0^2))},bereich[1],bereich[2],add=TRUE,col="green")
  
    
  if(weighted)
    fit = nls(thegaussian,daten,weights=1/err^2,start=list(C=ymin,N=ymax,mu=mu0,sig=sig0))
  else
    fit = nls(thegaussian,daten,start=list(C=ymin,N=ymax,mu=mu0,sig=sig0))
  
  return(fit)
  
}




  




plotgaus <- function(fit,bereich,log="",col="red"){ #--- Plotten der gefitteten Gaußfunktion in vorhandenen Graph
  fitdata <- summary(fit)$parameters
  
  N<-fitdata["N","Estimate"]
  C<-fitdata["C","Estimate"]
  mu<-fitdata["mu","Estimate"]
  sig<-fitdata["sig","Estimate"]
  
  plot (function(x){C + N*exp(-(x-mu)^2/(2*sig^2))},bereich[1],bereich[2],add=TRUE,col=col,log=log,lwd=2)
  
}




printfitdata <- function(fit,title=""){ #--- Ausgabe der Gaußfit-Daten
  fitdata <- summary(fit)$parameters

  N<-fitdata["N","Estimate"]
  sN<-fitdata["N","Std. Error"]
  C<-fitdata["C","Estimate"]
  sC<-fitdata["C","Std. Error"]
  
  mu<-fitdata["mu","Estimate"]
  smu<-fitdata["mu","Std. Error"]
  sigma<-fitdata["sig","Estimate"]
  ssigma<-fitdata["sig","Std. Error"]
  chisquare<-sum(((summary(fit))[[2]])^2)/(summary(fit)[[4]][[2]])
  
  cat(title)
  cat("&$")

  rN=roundfunc(c(N,sN))
  N=rN[1]
  sN=rN[2]
  
  rC=roundfunc(c(C,sC))
  C=rC[1]
  sC=rC[2]
  
  rmu=roundfunc(c(mu,smu))
  mu=rmu[1]
  smu=rmu[2]
  
  rsigma=roundfunc(c(sigma,ssigma))
  sig=rsigma[1]
  ssig=rsigma[2]
  
  #cat(" N    = ")
  cat(N)
  cat("\\pm")
  cat(sN)
  cat("$&$")
  
  #cat(" C    = ")
  cat(C)
  cat("\\pm")
  cat(sC)
  cat("$&$")
  
    
  #cat(" mu    = ")
  cat(mu)
  cat("\\pm")
  cat(smu)
  cat("$&$")
  
  #cat(" sigma = ")
  cat(sig)
  cat("\\pm")
  cat(ssig)
  cat("$&$")
  
  cat(round(chisquare,2))
  cat("$\\\\\n")

}

starttable <- function(){
  cat("\\begin{table}[h!]\n\\footnotesize\\centering\n\\begin{tabular}{|c||c|c|c|c||c|}\n\\hline\nEnergie / keV&$N/\\mathrm{s^{-1}}$&$C/\\mathrm{s^{-1}}$&$\\mu/\\mathrm{keV}$&$\\sigma/\\mathrm{keV}$&$\\chi^2$ / ndf\\\\\\hline\\hline")
}
starttableKanal <- function(){
  cat("\\begin{table}[h!]\n\\footnotesize\\centering\n\\begin{tabular}{|c||c|c|c|c||c|}\n\\hline\nEnergie / Kanal&$N/\\mathrm{s^{-1}}$&$C/\\mathrm{s^{-1}}$&$\\mu/\\mathrm{Kanal}$&$\\sigma/\\mathrm{Kanal}$&$\\chi^2$ / ndf\\\\\\hline\\hline")
}


endtable <- function(caption="",label=""){
  txtlabel=paste("\\label{",label,"}",sep="")
  if(label=="")
    txtlabel=""
  
  cat(paste("\\hline\\end{tabular}\n\\caption{",caption,txtlabel,"}\n\\end{table}\n",sep=""))
  
}



getresult <- function(fit){
  fitdata <- summary(fit)$parameters
  
  mu<-fitdata["mu","Estimate"]
  smu<-fitdata["mu","Std. Error"]
  sig<-fitdata["sig","Estimate"]
  ssig<-fitdata["sig","Std. Error"]
  chisquare<-sum(((summary(fit))[[2]])^2)/(summary(fit)[[4]][[2]])
  
  return(c(roundfunc(c(mu,abs(sig/2))),chisquare))
  
}
