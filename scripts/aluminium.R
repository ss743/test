source("readFiles.R")
source("functions.R")
source("expfit.R")

par(mar=c(5,5,1,1))

alu_names = c("0-0","0-5","1-0","1-5","2-0","2-5","3-0","3-5","4-0","4-5","5-0","5-5","6-0","6-5","7-0","7-5","8-0","8-5","9-0","9-5","10-0","10-5","11-0","11-5","12-0")
n = length(alu_names)
alu = array(dim=c(2048,n))
#alu[,1] = readTKA("aluminium/Alu_0-0mm")[[1]]

m=c(0,1,1,1,1,1,1,2,1,2,2,2,2,2,2,3,2,3,3,3,3,3,3,4,4)
s_di=0.1

s_d=sqrt(m)*s_di

for(i in c(1:n)){
  alu[,i] = readTKA(paste("aluminium/Alu_",alu_names[i],"mm",sep=""))[[1]]
}

alu_times=alu[1:2,]
alu[1:2,]=0

alu_rates=alu/alu_times[1,]

alu_sum2=c()
alu_sum=c()

for(i in c(1:n)){
  alu_sum[i]=sum(alu_rates[,i])
  alu_sum2[i]=sum(alu[,i])
}

alu_err=alu_sum/sqrt(alu_sum2)

x=c(1:n)*0.5-0.5
bereich=c(0,12)
drawCIx(x,y=alu_sum,sy=alu_err,sx=s_d,scol="black",xlab=expression(d/mm),ylab=expression(r/s^-1))
fit=dexpfit(data.frame(x=x,y=alu_sum,sy=alu_err),bereich,weighted=TRUE)
plotdexp(fit,c(0,12))
printdexpdata(fit)

C<-fit["C","Estimate"]
sC<-C*sqrt((fit["C","Std. Error"]/fit["C","Estimate"])^2)

untergrund=roundfunc(c(C,sC))

cat("\nUntergrund:\nr_U = (")
cat(paste(untergrund[1]," +- ",untergrund[2],") 1/s",sep=""))
