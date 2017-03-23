drawCI <- function(x,y,sy,xlab="x",ylab="y",col="black",scol="darkgrey"){
  plot(x,y,pch=4,cex=0.6,bty="l",col=col,xlab=xlab,ylab=ylab)
  arrows(x,y,x,y-sy,cex=0.6,pch=4,bty="l",col=scol,length=0.05,angle=90)
  arrows(x,y,x,y+sy,cex=0.6,pch=4,bty="l",col=scol,length=0.05,angle=90)
  points(x,y,cex=0.6,pch=4,col=col)
  grid()
}
drawCIlim <- function(x,y,sy,xlab="x",ylab="y",col="black",scol="darkgrey",xlim=c(0,60),ylim=c(0,800)){
  plot(x,y,pch=4,cex=0.6,bty="l",col=col,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
  arrows(x,y,x,y-sy,cex=0.6,pch=4,bty="l",col=scol,length=0.05,angle=90)
  arrows(x,y,x,y+sy,cex=0.6,pch=4,bty="l",col=scol,length=0.05,angle=90)
  points(x,y,cex=0.6,pch=4,col=col)
  grid()
}
drawCIx <- function(x,y,sy,sx,xlab="x",ylab="y",col="black",scol="darkgrey"){
  plot(x,y,pch=4,cex=0.6,bty="l",col=col,xlab=xlab,ylab=ylab)
  arrows(x,y,x,y-sy,cex=0.6,pch=4,bty="l",col=scol,length=0.05,angle=90)
  arrows(x,y,x,y+sy,cex=0.6,pch=4,bty="l",col=scol,length=0.05,angle=90)
  arrows(x,y,x-sx,y,cex=0.6,pch=4,bty="l",col=scol,length=0.05,angle=90)
  arrows(x,y,x+sx,y,cex=0.6,pch=4,bty="l",col=scol,length=0.05,angle=90)
  points(x,y,cex=0.6,pch=4,col=col)
  grid()
}


draw <- function(x,y,xlab="x",ylab="y",col="black",scol="darkgrey"){
  plot(x,y,pch=4,cex=0.6,bty="l",col=col,xlab=xlab,ylab=ylab)
  grid()
}

roundfunc <- function(vals){
  x=vals[1]
  xerr=vals[2]
  #print(vals)
  n=0
  for(i in -20:20){
    a=round(xerr,i)*10^i
    #print(a)
    if(a==1){
      n=i+1
      return(c(round(x,n),round(xerr,n)))
    }
    if(a==2){
      if(xerr*10^i<1.95){
        n=i+1
      } else {
        n=i
      }
      return(c(round(x,n),round(xerr,n)))
    }
    if(a>2){
      n=i
      return(c(round(x,n),round(xerr,n)))
    }
  }
  return(vals)
  
}

toenergy <- function(val){
  c=2.998*10^8*10^3#mm
  E0=14.412*10^3
  
  return(E0*val/c)
}

add <- function(vec1,vec2){
  res=c()
  res[1]=vec1[1]+vec2[1]
  res[2]=sqrt(vec1[2]^2+vec2[2]^2)
  return(res)
}
substract <- function(vec1,vec2){
  res=c()
  res[1]=vec1[1]-vec2[1]
  res[2]=sqrt(vec1[2]^2+vec2[2]^2)
  return(res)
}
divide <- function(vec1,vec2){
  res=c()
  res[1]=vec1[1]/vec2[1]
  res[2]=abs(res[1])*sqrt((vec1[2]/vec1[1])^2+(vec2[2]/vec2[1])^2)
  return(res)
}
multiplicate <- function(vec1,vec2){
  res=c()
  res[1]=vec1[1]*vec2[1]
  res[2]=abs(res[1])*sqrt((vec1[2]/vec1[1])^2+(vec2[2]/vec2[1])^2)
  return(res)
}

mean <- function(vec1,vec2){
  res=c()
  res[1]=(vec1[1]/vec1[2]^2+vec2[1]/vec2[2]^2)/(1/vec1[2]^2+1/vec2[2]^2)
  res[2]=sqrt(1/(1/vec1[2]^2+1/vec2[2]^2))
  
  return(res)
}

mean3 <- function(vec1,vec2,vec3){
  res=c()
  res[1]=(vec1[1]/vec1[2]^2+vec2[1]/vec2[2]^2+vec3[1]/vec3[2]^2)/(1/vec1[2]^2+1/vec2[2]^2+1/vec3[2]^2)
  res[2]=sqrt(1/(1/vec1[2]^2+1/vec2[2]^2+1/vec3[2]^2))
  
  return(res)
}