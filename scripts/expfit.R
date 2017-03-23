expfit <- function(input, bereich, weighted=FALSE){ #--- Fitten der Exponentialfunktion
  
  theexponential <- y ~ A + C *exp(lambda*x)
  
  #daten=input[bereich[1]:bereich[2],]
  daten=subset(input,x>=bereich[1] & x <= bereich[2])
  
  #---Startwerte
  ymin=min(daten$y)
  xmin=daten$x[which.min(daten$y)]
  ymax=max(daten$y)
  xmax=daten$x[which.max(daten$y)]
  lambda_est=(log(ymin)-log(ymax))/(xmin-xmax)
  C_est=ymin*exp(-lambda_est*xmin)
  A_est=0
  err=daten$sy
  startvalues=list(C=C_est,A=A_est,lambda=lambda_est)
  #---Startwerte
  
  #---Durchführen des Fits
  nlc<-nls.control(maxiter=5000)
  if(weighted)
    fit = nls(theexponential,daten,weights=1/err^2,start=startvalues,control=nlc)
  else
    fit = nls(theexponential,daten,start=startvalues,control=nlc)
  
  chiquadratndf=sum(summary(fit)[[2]]^2)/summary(fit)[[4]][[2]]
  fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
  
  return(fitdata)
  
}

dexpfit <- function(input, bereich, weighted=FALSE){ #--- Fitten der Exponentialfunktion
  
  theexponential <- y ~ C *exp(lambda*x) + D*exp(mu*x)
  
  #daten=input[bereich[1]:bereich[2],]
  daten=subset(input,x>=bereich[1] & x <= bereich[2])
  
  #---Startwerte
  ymin=min(daten$y)
  xmin=daten$x[which.min(daten$y)]
  ymax=max(daten$y)
  xmax=daten$x[which.max(daten$y)]
  lambda_est=(log(ymin)-log(ymax))/(xmin-xmax)
  C_est=ymin*exp(-lambda_est*xmin)
  D_est=C_est
  mu_est=-1
  A_est=0
  err=daten$sy
  startvalues=list(C=C_est,lambda=lambda_est,D=D_est,mu=mu_est)
  #---Startwerte
  
  #plot (function(x){A_est + C_est *exp(lambda_est*x) + D_est*exp(mu_est*x)},bereich[1],bereich[2],add=TRUE,col="green")
  
  
  #---Durchführen des Fits
  nlc<-nls.control(maxiter=5000)
  if(weighted)
    fit = nls(theexponential,daten,weights=1/err^2,start=startvalues,control=nlc)
  else
    fit = nls(theexponential,daten,start=startvalues,control=nlc)
  
  
  chiquadratndf=sum(summary(fit)[[2]]^2)/summary(fit)[[4]][[2]]
  fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
  
  return(fitdata)
  
}

plotexp <- function(fitdata,bereich){ #--- Plotten der gefitteten Exponentialfunktion in vorhandenen Graph
  
  lambda<-fitdata["lambda","Estimate"]
  C<-fitdata["C","Estimate"]
  A<-fitdata["A","Estimate"]
  
  plot (function(x){A + C *exp(lambda*x)},bereich[1],bereich[2],add=TRUE,col="red")
  
}

plotdexp <- function(fitdata,bereich){ #--- Plotten der gefitteten Exponentialfunktion in vorhandenen Graph
  
  lambda<-fitdata["lambda","Estimate"]
  mu<-fitdata["mu","Estimate"]
  C<-fitdata["C","Estimate"]
  D<-fitdata["D","Estimate"]
  #A<-fitdata["A","Estimate"]
  
  plot (function(x){C *exp(lambda*x) + D*exp(mu*x)},bereich[1],bereich[2],add=TRUE,col="red")
  
}


printexpdata <- function(fitdata,title="",factor=1,error=0){ #--- Ausgabe der Exponentialfit-Daten
  
  lambda<-fitdata["lambda","Estimate"]/factor
  slambda<-lambda*sqrt((fitdata["lambda","Std. Error"]/fitdata["lambda","Estimate"])^2+(error/factor)^2)

  cat(title)
  cat("\n")
  
  cat(" Zerfallskonstante    = ")
  cat(lambda)
  cat("+-")
  cat(slambda)
  cat("\n")
  
  tau=1/lambda
  stau=1/lambda*slambda/lambda
  
  cat(" Mittlere Lebensdauer = ")
  cat(tau)
  cat("+-")
  cat(stau)
  cat("\n")
  
  T=log(2)*tau
  sT=log(2)*stau

  cat(" Halbwertszeit        = ")
  cat(T)
  cat("+-")
  cat(sT)
  cat("\n")
    
}

printdexpdata <- function(fitdata,title="",factor=1,error=0){ #--- Ausgabe der Exponentialfit-Daten
  
  lambda<-fitdata["lambda","Estimate"]/factor
  slambda<-lambda*sqrt((fitdata["lambda","Std. Error"]/fitdata["lambda","Estimate"])^2+(error/factor)^2)
  C<-fitdata["C","Estimate"]/factor
  sC<-C*sqrt((fitdata["C","Std. Error"]/fitdata["C","Estimate"])^2+(error/factor)^2)
  D<-fitdata["D","Estimate"]/factor
  sD<-D*sqrt((fitdata["D","Std. Error"]/fitdata["D","Estimate"])^2+(error/factor)^2)
  
  cat(title)
  cat("\n")
  cat(" Erste Funktion:")
  cat("\n")
  
  cat(" Faktor               = ")
  cat(C)
  cat("+-")
  cat(sC)
  cat("\n")
  
  cat(" Zerfallskonstante    = ")
  cat(lambda)
  cat("+-")
  cat(slambda)
  cat("\n")
  
  tau=1/lambda
  stau=1/lambda*slambda/lambda
  
  cat(" Mittlere Lebensdauer = ")
  cat(tau)
  cat("+-")
  cat(stau)
  cat("\n")
  
  T=log(2)*tau
  sT=log(2)*stau
  
  cat(" Halbwertszeit        = ")
  cat(T)
  cat("+-")
  cat(sT)
  cat("\n")
  
  cat("\n")
  cat(" Zweite Funktion:")
  cat("\n")
  
  cat(" Faktor               = ")
  cat(D)
  cat("+-")
  cat(sD)
  cat("\n")
  
  
  lambda<-fitdata["mu","Estimate"]/factor
  slambda<-lambda*sqrt((fitdata["mu","Std. Error"]/fitdata["mu","Estimate"])^2+(error/factor)^2)
  
  
  cat(" Zerfallskonstante    = ")
  cat(lambda)
  cat("+-")
  cat(slambda)
  cat("\n")
  
  tau=1/lambda
  stau=1/lambda*slambda/lambda
  
  cat(" Mittlere Lebensdauer = ")
  cat(tau)
  cat("+-")
  cat(stau)
  cat("\n")
  
  T=log(2)*tau
  sT=log(2)*stau
  
  cat(" Halbwertszeit        = ")
  cat(T)
  cat("+-")
  cat(sT)
  cat("\n")
  
  
}
