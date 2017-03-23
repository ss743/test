setwd("C:\\Users\\Saskia\\Documents\\Physik\\FP\\fp2-moessbauer\\scripts")
source("readFiles.R")
source("functions.R")
source("gausfit1.R")
source("supergaus.R")

par(mar=c(5,5,1,1))
Tb = readTKA("eichung/TbEichung")[[1]]
Rb = readTKA("eichung/RbEichung")[[1]]
Mo = readTKA("eichung/MoEichung")[[1]]
Cu = readTKA("eichung/CuEichung")[[1]]
Ag = readTKA("eichung/AgEichung")[[1]]
Ba = readTKA("eichung/BaEichung")[[1]]

Tb_times=Tb[1:2]
Rb_times=Rb[1:2]
Mo_times=Mo[1:2]
Cu_times=Cu[1:2]
Ag_times=Ag[1:2]
Ba_times=Ba[1:2]


Tb[1:2]=0
Rb[1:2]=0
Mo[1:2]=0
Cu[1:2]=0
Ag[1:2]=0
Ba[1:2]=0

x=c(1:2048)

drawCI(x, Tb, sqrt(Tb), scol="black", xlab = "channel", ylab = "Counts")
bereich=c(450, 1000)
fit5=supergausfit(data.frame(x=x, y=Tb, sy=sqrt(Tb)), bereich,mu=593,sig0=55,N0=22500,N01=2500000,sig01=80,mu01=740, weighted=TRUE) 
plotsupergaus(fit5, bereich, col="red")
printsupergausdata(fit5)
drawCI(x, Rb, sqrt(Rb), scol="black", xlab = "channel", ylab = "Counts")
fit4=gausfit(data.frame(x=x,y=Rb, sy=sqrt(Rb)),c(100,218), weighted=TRUE)
plotgaus(fit4,c(100,218),col="red")
drawCI(x, Mo, sqrt(Mo), scol="black", xlab = "channel", ylab = "Counts")
fit3=gausfit(data.frame(x=x,y=Mo, sy=sqrt(Mo)),c(135,290), weighted=TRUE)
plotgaus(fit3,c(135,290),col="red")

#drawCI(x, Cu, sqrt(Cu), scol="black")

drawCI(x, Ag, sqrt(Ag), scol="black", xlab = "channel", ylab = "Counts")
fit2=gausfit(data.frame(x=x,y=Ag, sy=sqrt(Ag)),c(200,475), weighted=TRUE)
plotgaus(fit2,c(200,475),col="red")

drawCI(x, Ba, sqrt(Ba), scol="black", xlab = "channel", ylab = "Counts")
fit1=gausfit(data.frame(x=x,y=Ba, sy=sqrt(Ba)),c(300,550), weighted=TRUE)
plotgaus(fit1,c(300,600),col="red")
