library("minpack.lm")
library(RcppFaddeeva)#install.packages('RcppFaddeeva')

voigtfit <- function(input,A0,C0,mu0,sigma0,gamma0,weighted=FALSE) {
  
  y <- input$y
  x <- input$x
  
  thevoigt <- y ~ A * Voigt(x,mu,sigma,gamma) + C
  
  if(weighted)
    err=input$sy
  else
    err=1*input$y/input$y
  
  try({
    fit=nlsLM(thevoigt,input,weights=1/err^2,start=list(A=A0,C=C0,mu=mu0,sigma=sigma0,gamma=gamma0),control=nls.control(maxiter=100))
    chiquadratndf=sum(summary(fit)[[2]]^2)/summary(fit)[[4]][[2]]
    fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
    return(fitdata)
  })
}

sixvoigt <- function(input,A01,A02,A03,A04,A05,A06,C0,mu01,mu02,mu03,mu04,mu05,mu06,sigma01,sigma02,sigma03,sigma04,sigma05,sigma06,gamma01,gamma02,gamma03,gamma04,gamma05,gamma06,weighted=FALSE) {
  
  y <- input$y
  x <- input$x
  
  thevoigt <- y ~ A1 * Voigt(x,mu1,sigma1,gamma1) + A2 * Voigt(x,mu2,sigma2,gamma2) + A3 * Voigt(x,mu3,sigma3,gamma3) + A4 * Voigt(x,mu4,sigma4,gamma4) + A5 * Voigt(x,mu5,sigma5,gamma5) + A6 * Voigt(x,mu6,sigma6,gamma6) + C
  
  if(weighted)
    err=input$sy
  else
    err=1*input$y/input$y
  print(C0)
  plot(function(x){A01 * Voigt(x,mu01,sigma01,gamma01) + A02 * Voigt(x,mu02,sigma02,gamma02) + A03 * Voigt(x,mu03,sigma03,gamma03) + A04 * Voigt(x,mu04,sigma04,gamma04) + A05 * Voigt(x,mu05,sigma05,gamma05) + A06 * Voigt(x,mu06,sigma06,gamma06) + C0},-8.05,8.05,add=TRUE,col="yellow",n=10000)
  
  try({
    fit=nlsLM(thevoigt,input,weights=1/err^2,start=list(A1=A01,A2=A02,A3=A03,A4=A04,A5=A05,A6=A06,C=C0,mu1=mu01,mu2=mu02,mu3=mu03,mu4=mu04,mu5=mu05,mu6=mu06,sigma1=sigma01,sigma2=sigma02,sigma3=sigma03,sigma4=sigma04,sigma5=sigma05,sigma6=sigma06,gamma1=gamma01,gamma2=gamma02,gamma3=gamma03,gamma4=gamma04,gamma5=gamma05,gamma6=gamma06), control = list(maxiter = 100))
    chiquadratndf=sum(summary(fit)[[2]]^2)/summary(fit)[[4]][[2]]
    fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
    return(fitdata)
  })
}


plotvoigt<-function(fitdata,bereich,col="red",lwd=1,lty=1){
  A<-fitdata["A","Estimate"]
  C<-fitdata["C","Estimate"]
  mu<-fitdata["mu","Estimate"]
  sigma<-fitdata["sigma","Estimate"]
  gamma<-fitdata["gamma","Estimate"]
  
  try({plot(function(x){A * Voigt(x,mu,sigma,gamma) + C},bereich[1],bereich[2],add=TRUE,col=col,n=10000,lwd=lwd,lty=lty)})
  #print(optimize(function(x){A * Voigt(x,mu,sigma,gamma) + C},interval=c(-6,6)))
  #print(optimize(function(x){A * Voigt(x,mu,sigma,gamma) + C},interval=c(-6,6),maximum=TRUE))
  print(A * Voigt(mu,mu,sigma,gamma) + C)
}

plotsixvoigt<-function(fitdata,bereich,col="red",lwd=1,lty=1){
  C<-fitdata["C","Estimate"]
  A1<-fitdata["A1","Estimate"]
  A2<-fitdata["A2","Estimate"]
  A3<-fitdata["A3","Estimate"]
  A4<-fitdata["A4","Estimate"]
  A5<-fitdata["A5","Estimate"]
  A6<-fitdata["A6","Estimate"]
  mu1<-fitdata["mu1","Estimate"]
  mu2<-fitdata["mu2","Estimate"]
  mu3<-fitdata["mu3","Estimate"]
  mu4<-fitdata["mu4","Estimate"]
  mu5<-fitdata["mu5","Estimate"]
  mu6<-fitdata["mu6","Estimate"]
  sigma1<-fitdata["sigma1","Estimate"]
  sigma2<-fitdata["sigma2","Estimate"]
  sigma3<-fitdata["sigma3","Estimate"]
  sigma4<-fitdata["sigma4","Estimate"]
  sigma5<-fitdata["sigma5","Estimate"]
  sigma6<-fitdata["sigma6","Estimate"]
  gamma1<-fitdata["gamma1","Estimate"]
  gamma2<-fitdata["gamma2","Estimate"]
  gamma3<-fitdata["gamma3","Estimate"]
  gamma4<-fitdata["gamma4","Estimate"]
  gamma5<-fitdata["gamma5","Estimate"]
  gamma6<-fitdata["gamma6","Estimate"]
  
  try({plot(function(x){A1 * Voigt(x,mu1,sigma1,gamma1) + A2 * Voigt(x,mu2,sigma2,gamma2) + A3 * Voigt(x,mu3,sigma3,gamma3) + A4 * Voigt(x,mu4,sigma4,gamma4) + A5 * Voigt(x,mu5,sigma5,gamma5) + A6 * Voigt(x,mu6,sigma6,gamma6) + C},bereich[1],bereich[2],add=TRUE,col=col,n=10000,lwd=lwd,lty=lty)})
}
