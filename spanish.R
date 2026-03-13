library(fda.usc)
data(aemet)

temp = aemet$temp$data
wind = aemet$wind.speed$data
prec = aemet$logprec$data

tempdata=list()
winddata=list()
precdata=list()
for(i in 1:73){
  temp[i,] = smooth.spline(temp[i,],spar=0.5)$y
  wind[i,] = smooth.spline(wind[i,],spar=0.5)$y
  prec[i,] = smooth.spline(prec[i,],spar=0.5)$y
  tempdata[[i]] = temp[i,]
  winddata[[i]] = wind[i,]
  precdata[[i]] = prec[i,]
  
}
location = cbind(aemet$df$longitude,aemet$df$latitude)
location = scale_location(location)


restemp = fit_sdf(tempdata, location,r = 5, K = 10, degree = 3, diff_order = 2, dia = 0.3, gamma = 0.2,
                n_iter = 30, step_A = 1, step_Theta = 1e-3, Ucap = 30,
                cv = F, lambda1 = 2, lambda2 = 0.2 ,
                n_folds = 5,plotA = F,interp_lambda = NULL,verbose = TRUE)

restemp2 = fit_sdf(tempdata, location,r = 5, K = 10, degree = 3, diff_order = 2, dia = 0.3, gamma = 0.2,
                  n_iter = 30, step_A = 1, step_Theta = 1e-3, Ucap = 30,
                  cv = F, lambda1 = 2, lambda2 = 0 ,
                  n_folds = 5,plotA = F,interp_lambda = NULL,verbose = TRUE)

set.panel(1,4)
lon=aemet$df$longitude
lat=aemet$df$latitude
plotclu = function(x, y=0, z=0){
  lowx=min(x)
  hix=max(x)
  if(length(dim(y))!=0){
    lowy=min(y)
    hiy=max(y)
  }
  else{
    lowy=lowx+1
    hiy=hix-1
  }
  
  if(length(dim(z))!=0){
    lowz = min(z)
    hiz = max(z)
  }
  else{
    lowz = lowx+1
    hiz = hix-1
  }
  low = min(lowx, lowy, lowz)
  hi = max(hix, hiy, hiz)
  
  
  n = dim(x)[1]
  plot(x[1,],type='l', col=1, ylim=c(low, hi),xlab='',ylab='',xaxt='n',yaxt='n',main='SRCSE')
  title(xlab='Time: day', ylab='Temprature', 
        line=2, cex.lab=1.2)
  axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
  grid()
  for(i in 2:n){
    lines(x[i,],col=1)
  }
  if(length(dim(y))!=0){
    for(i in 1:dim(y)[1]){
      lines(y[i,],col=2)
    }
  }
  if(length(dim(z))!=0){
    for(i in 1:dim(z)[1]){
      lines(z[i,],col=3)
    }
  }
}
d = dist(t(restemp$sdfhat))
h = hclust(d, method = 'ward.D')
res = cutree(h,2)
plot(lon,lat,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n',main= 'SRCSE')
title(xlab='Lon', ylab='Lat', 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
plotclu(temp[which(res==1),],temp[which(res==2),])


plotclu = function(x, y=0, z=0){
  lowx=min(x)
  hix=max(x)
  if(length(dim(y))!=0){
    lowy=min(y)
    hiy=max(y)
  }
  else{
    lowy=lowx+1
    hiy=hix-1
  }
  
  if(length(dim(z))!=0){
    lowz = min(z)
    hiz = max(z)
  }
  else{
    lowz = lowx+1
    hiz = hix-1
  }
  low = min(lowx, lowy, lowz)
  hi = max(hix, hiy, hiz)
  
  
  n = dim(x)[1]
  plot(x[1,],type='l', col=1, ylim=c(low, hi),xlab='',ylab='',xaxt='n',yaxt='n',main='SE')
  title(xlab='Time: day', ylab='Temprature', 
        line=2, cex.lab=1.2)
  axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
  grid()
  for(i in 2:n){
    lines(x[i,],col=1)
  }
  if(length(dim(y))!=0){
    for(i in 1:dim(y)[1]){
      lines(y[i,],col=2)
    }
  }
  if(length(dim(z))!=0){
    for(i in 1:dim(z)[1]){
      lines(z[i,],col=3)
    }
  }
}
d = dist(t(restemp2$sdfhat))
h = hclust(d, method = 'ward.D')
res = cutree(h,2)
res = 3-res
plot(lon,lat,col=res,pch=20,cex=1.5
     ,xlab='',ylab='',xaxt='n',yaxt='n',main = 'SE')
title(xlab='Lon', ylab='Lat', 
      line=2, cex.lab=1.2)
axis(side = 1, tck = -0.02) ;axis(side = 2, tck = -0.02)
grid()
plotclu(temp[which(res==1),],temp[which(res==2),])


