library(RcppFaddeeva)#install.packages('RcppFaddeeva')

source("readFiles.R")
source("functions.R")
source("expfit.R")
source("lorentzfit.R")
source("gausfit1.R")
source("voigtfit.R")

par(mar=c(5,5,1,1))

data = readTXT(paste("sechslinien/sechslinien0",sep=""))

data[[2]]=data[[2]]*10^-3
rate=data[[3]]/data[[2]]
sdata=sqrt(data[[3]])
srate=rate/sdata

untergrund=18.5
suntergrund=0.3
rate=rate-untergrund
srate=sqrt(srate^2+suntergrund^2)

k=1.243
sk=0.010
srate=rate*k*sqrt((srate/rate)^2+(sk/k)^2)
rate=rate*k

x=data[[1]]

drawCI(x,rate,srate,xlab=expression(v / mms^-1),ylab=expression(r / s^-1))
 fit=sixlorentz(data.frame(x=x,y=rate,sy=srate),neg=TRUE,weighted=TRUE,-5.23,-2.97,-0.51,0.92,3.24,5.43)
 plotsixlorentz(fit,c(-8.05,8.05),lwd=1.5,lty=3)
 fit2=sixgaus(data.frame(x=x,y=rate,sy=srate),c(-8.05,8.05),weighted=TRUE,-5,-3,-0.5,1,3,5)
 plotsix(fit2,c(-8.05,8.05),col="blue",lwd=1.5,lty=2)
 A01=fit2['N1','Estimate']
 A02=fit2['N2','Estimate']
 A03=fit2['N3','Estimate']
 A04=fit2['N4','Estimate']
 A05=fit2['N5','Estimate']
 A06=fit2['N6','Estimate']
 C0=28.6#fit2['C','Estimate']
 mu01=fit2['mu1','Estimate']
 mu02=fit2['mu2','Estimate']
 mu03=fit2['mu3','Estimate']
 mu04=fit2['mu4','Estimate']
 mu05=fit2['mu5','Estimate']
 mu06=fit2['mu6','Estimate']
 sigma01=fit2['sig1','Estimate']
 sigma02=fit2['sig2','Estimate']
 sigma03=fit2['sig3','Estimate']
 sigma04=fit2['sig4','Estimate']
 sigma05=fit2['sig5','Estimate']
 sigma06=fit2['sig6','Estimate']
 gamma=0.1
 gamma01=fit['tau1','Estimate']#gamma#fit['omega1','Estimate']
 gamma02=fit['tau2','Estimate']#gamma#fit['omega2','Estimate']
 gamma03=fit['tau3','Estimate']#gamma#fit['omega3','Estimate']
 gamma04=fit['tau4','Estimate']#gamma#fit['omega4','Estimate']
 gamma05=fit['tau5','Estimate']#gamma#fit['omega5','Estimate']
 gamma06=fit['tau6','Estimate']#gamma#fit['omega6','Estimate']
 #gamma01=gamma#fit['omega1','Estimate']
 #gamma02=gamma#fit['omega2','Estimate']
 #gamma03=gamma#fit['omega3','Estimate']
 #gamma04=gamma#fit['omega4','Estimate']
 #gamma05=gamma#fit['omega5','Estimate']
 #gamma06=gamma#fit['omega6','Estimate']
 fit3=sixvoigt(data.frame(x=x,y=rate,sy=srate),A01,A02,A03,A04,A05,A06,C0,mu01,mu02,mu03,mu04,mu05,mu06,sigma01,sigma02,sigma03,sigma04,sigma05,sigma06,gamma01,gamma02,gamma03,gamma04,gamma05,gamma06,weighted=TRUE)
 plotsixvoigt(fit3,c(-8.05,8.05),col="green",lwd=2.5,lty=1)

 legend(0,10.8,c("Lorentzfit","Gaussfit","Voigtfit"),col=c("red","blue","green"),lty=c(3,2,1),lwd=c(1.5,1.5,2.5))
 
 mu1=c(fit3['mu1','Estimate'],fit3['mu1','Std. Error'])
 mu2=c(fit3['mu2','Estimate'],fit3['mu2','Std. Error'])
 mu3=c(fit3['mu3','Estimate'],fit3['mu3','Std. Error'])
 mu4=c(fit3['mu4','Estimate'],fit3['mu4','Std. Error'])
 mu5=c(fit3['mu5','Estimate'],fit3['mu5','Std. Error'])
 mu6=c(fit3['mu6','Estimate'],fit3['mu6','Std. Error'])
 isomerie1=c()
 isomerie2=c()
 isomerie3=c()
 isomerie1[1]=toenergy((mu4[1]+mu3[1])/2)
 isomerie2[1]=toenergy((mu5[1]+mu2[1])/2)
 isomerie3[1]=toenergy((mu6[1]+mu1[1])/2)
 isomerie1[2]=toenergy(sqrt(mu4[2]^2+mu3[2]^2)/2)
 isomerie2[2]=toenergy(sqrt(mu5[2]^2+mu2[2]^2)/2)
 isomerie3[2]=toenergy(sqrt(mu6[2]^2+mu1[2]^2)/2)
 
 isomerie=c()
 isomerie[1]=(isomerie1[1]/isomerie1[2]^2+isomerie2[1]/isomerie2[2]^2+isomerie3[1]/isomerie3[2]^2)/(1/isomerie1[2]^2+1/isomerie2[2]^2+1/isomerie3[2]^2)
 isomerie[2]=sqrt(1/(1/isomerie1[2]^2+1/isomerie2[2]^2+1/isomerie3[2]^2))
 

 E1=roundfunc(substract(toenergy(mu1),isomerie))
 E2=roundfunc(substract(toenergy(mu2),isomerie))
 E3=roundfunc(substract(toenergy(mu3),isomerie))
 E4=roundfunc(substract(toenergy(mu4),isomerie))
 E5=roundfunc(substract(toenergy(mu5),isomerie))
 E6=roundfunc(substract(toenergy(mu6),isomerie))
 
 E23=divide(E2,E3)
 E13=divide(E1,E3)
 E12=divide(E1,E2)
 E54=divide(E5,E4)
 E64=divide(E6,E4)
 E65=divide(E6,E5)
 
 mug=c(0.090604,0.000009)#muN
 one=c(1,0)
 mua23=multiplicate(3*one,multiplicate(mug,divide(substract(one,E23),add(one,E23))))
 mua54=multiplicate(3*one,multiplicate(mug,divide(substract(one,E54),add(one,E54))))
 
 mua12=multiplicate(one,multiplicate(mug,divide(substract(E12,one),substract(E12/3,one))))
 mua65=multiplicate(one,multiplicate(mug,divide(substract(E65,one),substract(E65/3,one))))

 mua13=multiplicate(one,multiplicate(mug,divide(substract(one,E13),add(one,1/3*E13))))
 mua64=multiplicate(one,multiplicate(mug,divide(substract(one,E64),add(one,1/3*E64))))
 
 mua1=roundfunc(mean(mua23,mua54))
 mua2=roundfunc(mean(mua13,mua64))
 mua3=roundfunc(mean(mua12,mua65))
 
 mua=roundfunc(mean3(mua1,mua2,mua3))
 
 hquer=6.582*10^(-16)
 e=1.602*10^(-19)
 mproton=1.672*10^(-27)
 mun=e*hquer/(2*mproton)
 
 B1=divide(E1,multiplicate(substract(mua,mug),mun*one))
 B6=divide(E6,multiplicate(substract(mug,mua),mun*one))
 B2=divide(E2,multiplicate(substract(1/3*mua,mug),mun*one))
 B5=divide(E5,multiplicate(substract(mug,1/3*mua),mun*one))
 B3=substract(c(0,0),divide(E3,multiplicate(add(1/3*mua,mug),mun*one)))
 B4=divide(E4,multiplicate(add(mug,1/3*mua),mun*one))
 
 B16=roundfunc(mean(B1,B6))
 B25=roundfunc(mean(B2,B5))
 B34=roundfunc(mean(B3,B4))
 
 B=roundfunc(mean3(B16,B25,B34))