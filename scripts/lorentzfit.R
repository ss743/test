library("minpack.lm")

lorentzfit<-function(input,neg=FALSE,weighted=FALSE){
  
  #lorentz <- y~D+C/2*(2-(2*(x-omega)*tau)^2/(1+(2*(x-omega)*tau)^2))
  lorentz <- y~D  + C*1/(2*pi)*tau/((x-omega)^2+(1/2*tau)^2)
  
  D0=max(input$y)
  C0=min(input$y)-max(input$y)
  omega0=input$x[which.max(input$y)]
  tau0=10^(0.18)
  
  if(neg){
    C0=max(input$y)-min(input$y)
    D0=max(input$y)
    omega0=input$x[which.min(input$y)]
    
  }
  
  
  #print(input)
  #print(c(D0,C0,omega0,tau0))
  
  #plot(function(x){D0+C0*((2*(x-omega0)*tau0)^2/(1+(2*(x-omega0)*tau0)^2))},-6,6,add=TRUE,col="green")
  if(weighted)
    err=input$sy
  else
    err=1*input$y/input$y
  
  try({
    fit=nlsLM(lorentz,input,weights=1/err^2,start=list(D=D0,C=C0,omega=omega0,tau=tau0))
    #chiquadratndf=sum(residuals(fit)^2/abs(fitted(fit)))/summary(fit)$df[2]
    chiquadratndf=sum(summary(fit)[[2]]^2)/summary(fit)[[4]][[2]]
    fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
    return(fitdata)
  })
  
  return(NULL)
}

sixlorentz<-function(input,neg=FALSE,weighted=FALSE,omega01,omega02,omega03,omega04,omega05,omega06){
  
  lorentz <- y~D  + C1*1/(2*pi)*tau1/((x-omega1)^2+(1/2*tau1)^2)  + C2*1/(2*pi)*tau2/((x-omega2)^2+(1/2*tau2)^2)  + C3*1/(2*pi)*tau3/((x-omega3)^2+(1/2*tau3)^2)  + C4*1/(2*pi)*tau4/((x-omega4)^2+(1/2*tau4)^2)  + C5*1/(2*pi)*tau5/((x-omega5)^2+(1/2*tau5)^2)  + C6*1/(2*pi)*tau6/((x-omega6)^2+(1/2*tau6)^2)
  #+C1/2*(2-(2*(x-omega1)*tau1)^2/(1+(2*(x-omega1)*tau1)^2))+C2/2*(2-(2*(x-omega2)*tau2)^2/(1+(2*(x-omega2)*tau2)^2))+C3/2*(2-(2*(x-omega3)*tau3)^2/(1+(2*(x-omega3)*tau3)^2))+C4/2*(2-(2*(x-omega4)*tau4)^2/(1+(2*(x-omega4)*tau4)^2))+C5/2*(2-(2*(x-omega5)*tau5)^2/(1+(2*(x-omega5)*tau5)^2))+C6/2*(2-(2*(x-omega6)*tau6)^2/(1+(2*(x-omega6)*tau6)^2))
  
  
  D0=max(input$y)
  C0=min(input$y)-max(input$y)
  omega0=input$x[which.max(input$y)]
  tau0=10^(0.18)
  
  if(neg){
    C0=max(input$y)-min(input$y)
    D0=max(input$y)-0.0045
    omega0=input$x[which.min(input$y)]
    
  }
  
  #print(input)
  #print(c(D0,C0,omega0,tau0))
  
  #plot(function(x){D0+C0*((2*(x-omega0)*tau0)^2/(1+(2*(x-omega0)*tau0)^2))},-6,6,add=TRUE,col="green")
  if(weighted)
    err=input$sy
  else
    err=1*input$y/input$y
  
  try({
    fit=nlsLM(lorentz,input,weights=1/err^2,start=list(D=D0,C1=C0,C2=C0,C3=C0,C4=C0,C5=C0,C6=C0,omega1=omega01,omega2=omega02,omega3=omega03,omega4=omega04,omega5=omega05,omega6=omega06,tau1=tau0,tau2=tau0,tau3=tau0,tau4=tau0,tau5=tau0,tau6=tau0))
    #chiquadratndf=sum(residuals(fit)^2/abs(fitted(fit)))/summary(fit)$df[2]
    chiquadratndf=sum(summary(fit)[[2]]^2)/summary(fit)[[4]][[2]]
    fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
    return(fitdata)
  })
  
  return(NULL)
}
dispersionsfit<-function(input){
  
  disp <- y~D+C*((2*(x-omega)*tau)/(1+(2*(x-omega)*tau)^2)^2)
  
  D0=(min(input$y)+max(input$y))/2
  C0=max(input$y)-min(input$y)
  omega0=(input$x[which.max(input$y)]+input$x[which.min(input$y)])/2
  tau0=10^(-7)
  

  #print(c(D0,C0,omega0,tau0))
  
  #plot(function(x){D0+C0*((2*(x-omega0)*tau0)/(1+(2*(x-omega0)*tau0)^2)^2)},-6*10^(7),6*10^(7),add=TRUE,col="green")
  
  try({
    fit=nls(disp,input,start=list(D=D0,C=C0,omega=omega0,tau=tau0))
    chiquadratndf=sum(residuals(fit)^2/abs(fitted(fit)))/summary(fit)$df[2]
    fitdata=rbind(summary(fit)$parameters,c(chiquadratndf,0,0,0))
    return(fitdata)
  })
  
  return(NULL)
}


plotlorentz<-function(fitdata,bereich,lwd=1,lty=1){
  D<-fitdata["D","Estimate"]
  C<-fitdata["C","Estimate"]
  omega<-fitdata["omega","Estimate"]
  tau<-fitdata["tau","Estimate"]
  
  try({plot(function(x){D + C*1/(2*pi)*tau/((x-omega)^2+(1/2*tau)^2)},-6,6,add=TRUE,col="red",n=10000,lwd=lwd,lty=lty)})
}

plotsixlorentz<-function(fitdata,bereich,lwd=1,lty=1){
  D<-fitdata["D","Estimate"]
  C1<-fitdata["C1","Estimate"]
  omega1<-fitdata["omega1","Estimate"]
  tau1<-fitdata["tau1","Estimate"]
  C2<-fitdata["C2","Estimate"]
  omega2<-fitdata["omega2","Estimate"]
  tau2<-fitdata["tau2","Estimate"]
  C3<-fitdata["C3","Estimate"]
  omega3<-fitdata["omega3","Estimate"]
  tau3<-fitdata["tau3","Estimate"]
  C4<-fitdata["C4","Estimate"]
  omega4<-fitdata["omega4","Estimate"]
  tau4<-fitdata["tau4","Estimate"]
  C5<-fitdata["C5","Estimate"]
  omega5<-fitdata["omega5","Estimate"]
  tau5<-fitdata["tau5","Estimate"]
  C6<-fitdata["C6","Estimate"]
  omega6<-fitdata["omega6","Estimate"]
  tau6<-fitdata["tau6","Estimate"]
  
  try({plot(function(x){D  + C1*1/(2*pi)*tau1/((x-omega1)^2+(1/2*tau1)^2)  + C2*1/(2*pi)*tau2/((x-omega2)^2+(1/2*tau2)^2)  + C3*1/(2*pi)*tau3/((x-omega3)^2+(1/2*tau3)^2)  + C4*1/(2*pi)*tau4/((x-omega4)^2+(1/2*tau4)^2)  + C5*1/(2*pi)*tau5/((x-omega5)^2+(1/2*tau5)^2)  + C6*1/(2*pi)*tau6/((x-omega6)^2+(1/2*tau6)^2)},bereich[1],bereich[2],add=TRUE,col="red",n=10000,lwd=lwd,lty=lty)})
}


plotdisp<-function(fitdata,bereich){
  D<-fitdata["D","Estimate"]
  C<-fitdata["C","Estimate"]
  omega<-fitdata["omega","Estimate"]
  tau<-fitdata["tau","Estimate"]
  
  try({plot(function(x){D+C*((2*(x-omega)*tau)/(1+(2*(x-omega)*tau)^2)^2)},-6*10^(7),6*10^(7),add=TRUE,col="red")})
}

getlorentzvalue<-function(fitdata,x){
  D<-fitdata["D","Estimate"]
  C<-fitdata["C","Estimate"]
  omega<-fitdata["omega","Estimate"]
  tau<-fitdata["tau","Estimate"]

  return(D+C/2*(2-(2*(x-omega)*tau)^2/(1+(2*(x-omega)*tau)^2)))    
}

getdispvalue<-function(fitdata,x){
  D<-fitdata["D","Estimate"]
  C<-fitdata["C","Estimate"]
  omega<-fitdata["omega","Estimate"]
  tau<-fitdata["tau","Estimate"]
  
  return(D+C*((2*(x-omega)*tau)/(1+(2*(x-omega)*tau)^2)^2))    
}

getlorentzformula<-function(fitdata){
  D0<-fitdata["D","Estimate"]
  C0<-fitdata["C","Estimate"]
  omega0<-fitdata["omega","Estimate"]
  tau0<-fitdata["tau","Estimate"]
  sD<-fitdata["D","Std. Error"]
  sC<-fitdata["C","Std. Error"]
  somega<-fitdata["omega","Std. Error"]
  stau<-fitdata["tau","Std. Error"]
  
  D=roundfunc(c(D0,sD))
  C=roundfunc(c(C0,sC))
  omega=roundfunc(c(omega0,somega))
  tau=roundfunc(c(tau0,stau))
  
  txt=paste("(",D[1],"+-",D[2],") V + (",C[1]/2,"+-",C[2]/2,") V * (2-(2(omega_L-(",omega[1],"+-",omega[2],")Hz)*(",tau[1],"+-",tau[2],")s)^2/(1+(2(omega_L-(",omega[1],"+-",omega[2],")Hz)*(",tau[1],"+-",tau[2],")s)^2))",sep="")
  
  return(txt)
}

getdispformula<-function(fitdata){
  D0<-fitdata["D","Estimate"]
  C0<-fitdata["C","Estimate"]
  omega0<-fitdata["omega","Estimate"]
  tau0<-fitdata["tau","Estimate"]
  sD<-fitdata["D","Std. Error"]
  sC<-fitdata["C","Std. Error"]
  somega<-fitdata["omega","Std. Error"]
  stau<-fitdata["tau","Std. Error"]
  
  D=roundfunc(c(D0,sD))
  C=roundfunc(c(C0,sC))
  omega=roundfunc(c(omega0,somega))
  tau=roundfunc(c(tau0,stau))
  
  txt=paste("(",D[1],"+-",D[2],") V + (",C[1],"+-",C[2],") V * (2-(2(omega_L-(",omega[1],"+-",omega[2],")Hz)*(",tau[1],"+-",tau[2],")s)/(1+(2(omega_L-(",omega[1],"+-",omega[2],")Hz)*(",tau[1],"+-",tau[2],")s)^2)^2)",sep="")
  
  return(txt)
}