source("readFiles.R")
source("functions.R")
source("expfit.R")

par(mar=c(5,5,1,1))
empty = readTKA("plexiglas/leer")[[1]]
plexi = readTKA("plexiglas/Plexi")[[1]]

empty_times=empty[1:2]
plexi_times=plexi[1:2]

empty[1:2]=0
plexi[1:2]=0

rate_e=sum(empty)/empty_times[1]
rate_p=sum(plexi)/plexi_times[1]
s_e=rate_e/sqrt(sum(empty))
s_p=rate_p/sqrt(sum(plexi))

ratio=roundfunc(c(rate_p/rate_e,rate_p/rate_e*sqrt((s_e/rate_e)^2+(s_p/rate_p)^2)))
ratio2=roundfunc(c(rate_e/rate_p,rate_e/rate_p*sqrt((s_e/rate_e)^2+(s_p/rate_p)^2)))

cat(paste("Verhaeltnis r_p/r_e = ",ratio[1]," +- ",ratio[2],sep=""))
cat(paste("\n      k_p = r_e/r_p = ",ratio2[1]," +- ",ratio2[2],sep=""))