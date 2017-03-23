library(RcppFaddeeva)#install.packages('RcppFaddeeva')

source("readFiles.R")
source("functions.R")
source("expfit.R")
source("lorentzfit.R")
source("gausfit1.R")
source("voigtfit.R")

par(mar=c(5,5,1,1))

data = readTXT(paste("einlinien/einlinien2",sep=""))

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
fit=lorentzfit(data.frame(x=x,y=rate,sy=srate),neg=TRUE,weighted=TRUE)
plotlorentz(fit,c(-6,6),lwd=1.5,lty=3)
fit2=gausfit(data.frame(x=x,y=rate,sy=srate),c(-6,6),weighted=TRUE,N0=-5)
plotgaus(fit2,c(-6,6),col="blue",lwd=1.5,lty=2)
A0=fit2['N','Estimate']
C0=fit2['C','Estimate']
mu0=fit2['mu','Estimate']
sigma0=fit2['sig','Estimate']
gamma0=fit['omega','Estimate']
fit3=voigtfit(data.frame(x=x,y=rate,sy=srate),A0,C0,mu0,sigma0,gamma0,weighted=TRUE)
plotvoigt(fit3,c(-6,6),col="green",lwd=2.5,lty=1)

legend(2,10,c("Lorentzfit","Gaussfit","Voigtfit"),col=c("red","blue","green"),lty=c(3,2,1),lwd=c(1.5,1.5,2.5))

sigma=fit3['sigma','Estimate']
ssigma=fit3['sigma','Std. Error']
gamma=fit3['gamma','Estimate']
sgamma=fit3['gamma','Std. Error']

fg=2*sigma*sqrt(2*log(2,exp(1)))
fl=gamma

fv=0.5346*fl+sqrt(0.2166*fl^2+fg^2)

sfg=2*ssigma*sqrt(2*log(2,exp(1)))
sfl=sgamma
sfv=fv*sqrt(0.02^2+(sfg/fg)^2+(sfl/fl)^2)

sigmaG=fit2['sig','Estimate']
ssigmaG=fit2['sig','Std. Error']
gammaL=fit['tau','Estimate']
sgammaL=fit['tau','Std. Error']

fG=2*sigmaG*sqrt(2*log(2,exp(1)))
fL=gammaL

sfG=2*ssigmaG*sqrt(2*log(2,exp(1)))
sfL=sgammaL
fv2=0.5346*fL+sqrt(0.2166*fL^2+fG^2)
sfv2=fv2*sqrt(0.02^2+(sfG/fG)^2+(sfL/fL)^2)

GAMMA0=roundfunc(toenergy(c(fL,sfL)))

W=c(1.78,0.06)
GAMMA=c()
GAMMA[1]=GAMMA0[1]/4/W[1]
GAMMA[2]=GAMMA[1]*sqrt((GAMMA0[2]/GAMMA0[1])^2+(W[2]/W[1])^2)

#GAMMA=roundfunc(toenergy(c(fl,sfl)))

hquer=6.582*10^(-16)#eV
tau=c()
tau[1]=hquer/GAMMA[1]
tau[2]=tau[1]*GAMMA[2]/GAMMA[1]

Thalb=tau*log(2,exp(1))

#drawCIlim(x,rate,srate,xlab=expression(v / mms^-1),ylab=expression(r / s^-1),xlim=c(-1,1),ylim=c(8.5,13))
#plotlorentz(fit,c(-6,6),lwd=1.5,lty=3)
#plotgaus(fit2,c(-6,6),col="blue",lwd=1.5,lty=2)
#plotvoigt(fit3,c(-6,6),col="green",lwd=2.5,lty=1)

mu=c(fit3['mu','Estimate'],fit3['mu','Std. Error'])
isomerie=toenergy(mu)