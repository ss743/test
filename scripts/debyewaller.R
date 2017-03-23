library(Bessel)
source('functions.R')

Zinf=c(12.84,0.08)
Zmu=c(9.06,0.10)
Z=c()
Z[1]=(Zinf[1]-Zmu[1])/Zinf[1]
Z[2]=Zmu[1]/Zinf[1]*sqrt((Zinf[2]/Zinf[1])^2+(Zmu[2]/Zmu[1])^2)
i=complex(imaginary=1)
y=c()
y[1]=1/Re(1-exp(-TA[1]/2)*BesselJ(i*TA[1]/2,0))
y[2]=Re(TA[2]/2*exp(-TA[1]/2)*(-BesselJ(i*TA[1]/2,0)+i*BesselJ(i*TA[1]/2,-1))/(1-exp(-TA[1]/2)*BesselJ(i*TA[1]/2,0))^2)
fQ=multiplicate(Z,y)
#fQ[2]=fQ[1]*sqrt((Z[2]/Z[1])^2)