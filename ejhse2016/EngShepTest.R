setwd("F:\\Dropbox\\Научка\\hedge ratio\\r_hedge") # Megatron

save.image("F:\\Dropbox\\Научка\\hedge ratio\\r_hedge\\mforecast.RData")

install.packages("rmgarch")
require(rmgarch)
install.packages("R.matlab")
require(R.matlab)

# read data from *.mat
data <- readMat("F:\\Dropbox\\Научка\\hedge ratio\\hedge_matlab\\start.mat")
y<-data$y
test<-list();
res<-matrix(0,length(y),4)
for (i in 1:length(y)) {
  test[[i]]<-DCCtest(y[[i]][,2:3], garchOrder = c(1,1), n.lags = c(1,2,5,10))
  res[i,]<-test[[i]]$p.value
}
rm(i)

writeMat("F:\\Dropbox\\Научка\\hedge ratio\\hedge_matlab\\tmp.mat",
         dm=dm)

# DCC
uniGarchSpec <- ugarchspec(mean.model=list(armaOrder=c(0,0)),
                                             variance.model=list(model="fGARCH",submodel="TGARCH",
                                                                 garchOrder=c(1,1)))
multiGarchSpec <- multispec(replicate(2, uniGarchSpec))
dccSpec <- dccspec(multiGarchSpec,dccOrder=c(1,1), model="DCC",distribution="mvt")
forc <- list()
for (i in 1:10) {
dccFit <- dccfit(dccSpec, y[[i]][,2:3], out.sample = 49)
dccForecast <- dccforecast(dccFit, n.ahead=1,n.roll=49)
s<-rcov(dccForecast)
s<-simplify2array(s)
forc[[i]] <- s
rm(s)
}

rm(i,uniGarchSpec, multiGarchSpec, dccSpec, dccFit, dccForecast)

# GO-GARCH
forc <- list()
ggSpec <-gogarchspec(mean.model = list(model = "AR"),
                     variance.model = list(model = "fGARCH", garchOrder = c(1,1), submodel = "TGARCH"),
                     distribution.model = "mvnorm")
for (i in 1:10) {
ggFit <- gogarchfit(ggSpec,  y[[i]][,2:3], out.sample = 49)
ggForc <- gogarchforecast(ggFit, n.ahead = 1, n.roll = 49)
s<-rcov(ggForc)
s<-simplify2array(s)
forc[[i]] <- s
rm(s)
}
rm(i,ggSpec,ggFit,ggForc)
forcGG <- forc
rm(forc)

# create an array
# frcDcc <- vector("list", 10)
frcDcc <- array(0,dim=c(2,2,50,10))
for (j in 1:10) {
  for (i in 1:50) {
    frcDcc[,,i,j] <- as.matrix(forcDcc[[j]][,,,i])  
  }
}
rm(i,j)


loss2 <- function(s,sf) {
  loss2 <- rep(0,dim(sf)[3])
  for (i in 1:dim(sf)[3]) {
    s. <- sf[,,i]^(-1)%*%s[,,i]
    loss2[i] <- sum(diag(s.))-log(norm(s.))-ncol(s)  
  }
  return(loss2)
}

sTrue <- vector("list",10)
for (i in 1:10) {
  sTrue[[i]] <- array(data=NA,dim=c(2,2,50))
  ly <- nrow(y[[i]])-50
  for (j in 1:50) {    
#     browser()
    sTrue[[i]][,,j] <- t(y[[i]][(ly+j-10):(ly+j),2:3])%*%y[[i]][(ly+j-10):(ly+j),2:3]
  }
}
rm(i,j,ly)

lGg <- vector("list",10)
for (i in 1:10) {
  lGg[[i]]<-loss2(sTrue[[i]],frcGg[,,,i])
}
rm(i)

# Diebold-Mariano
d<-vector("list",10)
for (i in 1:10) {
d[[i]]<-lDcc[[i]]-lGg[[i]]
}
rm(i)
dm <-sapply(d,mean)

require(sandwich)
v <- sapply(d,lrvar, type = "Andrews")

dm<-dm/sqrt(v)

