# hook up packages
require(xts) # time series objects
require(rugarch) # univariate garch
require(rmgarch) # multivariate garch
# require(parmaMy) # portfolio optimization with 3-moment utility
require(parma) # portfolio optimization with 3-moment utility
require(abind) # abind()
require(PerformanceAnalytics)
# require(moments) # skewness()
options(repos = c(CRAN = "http://cran.rstudio.com")) # when Rstudio not want to install packages

install.packages("F:/r_packages_sources/parma_1.5-3.tar.gz", repos = NULL, type="source")

setwd('d:\\Google Диск\\Nauka\\mvs_portfolio\\mvs_portfolio')
setwd('f:/1work/mvs_portfolio_local/mvs_portfolio/mvs_portfolio/')
load('mvs12.RData') 
# load("mvs_start.RData")
source('mvs_functions.r')


e1 <- new.env()
load(file = 'mvs8.Rdata', envir = e1)

ls(envir = e1)
da_yah <- get('da_yah', envir = e1)

# read returns of assets ####
names <- c('aflt', 'avaz', 'gazpn', 'gmkn', 'lkoh', 'rtkm', 'sber', 'sngs', 'tatn')
ret <- readAll(names) # from Finam, Yahoo only from 2010
ret <- lapply(ret, function(x) {data.frame(date=x[-1,1], close=diff(log(x)))})
ret <- lapply(ret, function(x) {xts(x = x[,'close'], order.by = x[,'date'])})
ret1 <- xts(NULL)
for (i1 in 1:length(ret)) {
  ret1 <- merge.xts(ret1, ret[[i1]])
}
rm(i1)
ret <- ret1[complete.cases(ret1),] # remove rows with NA
rm(ret1)
ret <- diff(log(ret))
ret <- ret[complete.cases(ret),]
colnames(ret) <- names  
  
n <- ncol(ret) # # of assetes 
ts <- nrow(ret) # length of series (monthly data)
oos <- 6 # out of sample, half a year
isa <- ts-oos # length of in-sample

# descriptive statistics ####
ds <- NULL
# require(e1071)
for (i1 in 1:n) {
  ds <- cbind(ds, summStat(ret6[,i1])) # function from hedge_ratio work (aprilConf.R)
}
names <- colnames(ret6)
names <- gsub('\\.ME', '', names)
colnames(ds) <- names
ds <- ds[-1,]
rm(i1)

chart.Boxplot(ret_old)
# ret <- ret_old[ret_old$sber > -2 & ret_old$avaz > -2,]
ret <- ret_old # adjust for stocks split
ret$sber['2007-08-01'] <- ret$sber['2007-08-01']/1000
ret$avaz['2007-11-01'] <- ret$avaz['2007-11-01']/100
chart.Boxplot(ret)

# mgarch ####
oos <- 50 # out of sample
mgspec <- gogarchspec(mean.model = list(model='constant', lag.max=3, lag.criterion='HQ'),
                      variance.model = list(model='sGARCH', garchOrder=c(2,2)),
                      distribution.model='mvnorm')
gfit <- gogarchfit(mgspec, ret5, out.sample = oos)

mgspec.sk <- gogarchspec(mean.model = list(model='AR', lag.max=3, lag.criterion='AIC'),
                         variance.model = list(model='sGARCH', garchOrder=c(1,1)),
                         distribution.model='manig')

gfit.sk <- gogarchfit(mgspec.sk, ret5, out.sample = oos)#,solver.control = list(rho=1, delta=1e-7, tol=1e-6, trace=0))
t1 <- coef(gfit.sk)
sum(abs(t1['skew',]) > t1['shape',]) # wrong estimates

# utility maximization ####

# portfolio controls
n <- ncol(ret5)
ts <- nrow(ret5) # length of series (monthly data)
isa <- ts-oos # length of in-sample
LB <- rep(-.9, n) # -.9 for short
UB <- rep(.9, n)
bu <- 1 # leverage constraint (bu=1 is full investment constraint)
ra <- c(3, 5, 10, 20, 40, 80) # absolute risk aversion # c(5, 10) #
ra <- c(0.1, 1.5, 5, 10, 20, 40)
ra <- seq(1e-6, 2.5, length.out = 4)
ra <- seq(1e-6, 5, length.out = 20)
ra <- c(0.1, 0.5, 1, 2, 5, 10)
# 2 mom: good: ra <- 80; 3 mom:  ra= from 3 to 5 is ok; 
# it seems that short is worse for 3-mom

ggplot(data=data.frame(x=1:length(ra), ra=ra), aes(x=x, y=ra)) + 
  geom_bar(stat="identity") + 
  theme_bw()


# CARA, 2 Moment Approximation, in-sample, long/short
# fitted and forecasted moments
gmean <- fitted(gfit)
gcov <- rcov(gfit)
ret_mom = fmoments(Data = ret5, n.ahead = 1, roll = oos, spec = mgspec, 
                   solver = "hybrid", cluster = NULL, gfun = "tanh", maxiter1 = 100000, 
                   rseed = 500)
# output variables
w.go <- array(NA, dim = c(ts, n, length(ra))) # weights, for both in- and out-of-sample
mv.go <- array(NA, dim = c(ts, 4, length(ra))) # mean, sd, utility of portfolio; parmastatus
dimnames(mv.go) <- list(c(rep('isa', isa), rep('oos', oos)), 
                        c('portf.mean', 'portf.sd', 'util', 'status'), ra)
# portfolio estimation
for (i1 in 1:length(ra)) {
  putil.go <- vector(mode='list', length = ts)
  for(i in 1:isa){
    putil.go[[i]] = parmautility(U = "CARA", method = "moment",
                                 M1 = gmean[i,], M2 =  gcov[,,i],
                                 RA = ra[i1], budget = bu, 
                                 LB = LB, UB = UB)
    w.go[i,,i1] = weights(putil.go[[i]])
    mv.go[i,'portf.mean',i1] <- parmareward(putil.go[[i]])
    mv.go[i,'portf.sd',i1] <- sqrt(w.go[i,,i1] %*% gcov[,,i] %*% w.go[i,,i1])
  }
  mv.go[1:isa,'util',i1] <- cara2(pm1 = mv.go[1:isa,1,i1], pm2 = mv.go[1:isa,2,i1]^2, l=ra)
  rm(i)
  
# CARA, 2 Moment Approximation, out-of-sample
  for(i in (isa+1):ts){
    putil.go[[i]] = parmautility(U = "CARA", method = "moment", 
                                 M1 = ret_mom@moments$forecastMu[i-isa,], 
                                 M2 =  ret_mom@moments$forecastCov[,,i-isa],
                                 RA = ra[i1], budget = bu, 
                                 LB = LB, UB = UB)
    w.go[i,,i1] = weights(putil.go[[i]])
    mv.go[i,'portf.mean',i1] <- parmareward(putil.go[[i]])
    mv.go[i,'portf.sd',i1] <- sqrt(w.go[i,,i1] %*% ret_mom@moments$forecastCov[,,i-isa] %*% w.go[i,,i1])
  }
  mv.go[(isa+1):ts,'util',i1] <- cara2(pm1 = mv.go[(isa+1):ts,1,i1], pm2 = mv.go[(isa+1):ts,2,i1]^2, l=ra[i1])
  mv.go[,'status',i1] <- unlist(lapply(putil.go, parmastatus))
}
rm(i, i1, putil.go)

plot(as.numeric(mv.go[,'status',]), type='l')

# CARA, 3 Moment Approximation, in-sample, long/short

# fitted and forecasted moments
gmean.sk <- fitted(gfit.sk)#[,1:n]
gcov.sk <- rcov(gfit.sk)#[1:(n),1:(n),]
# 3rd and 4th momrnts are standardized by default
# gskew <- abind(rcoskew(gfit.sk, from=1, to=100, standardize=T)[1:n,1:(n^2),], 
               # rcoskew(gfit.sk, from=101, to=isa, standardize=T)[1:n,1:(n^2),])
# gkurt <- abind(rcokurt(gfit.sk, from=1, to=100)[1:n,1:(n^3),], 
               # rcokurt(gfit.sk, from=101, to=isa)[1:n,1:(n^3),]) #array(0, dim = c(n, n^3, ts))
gskew <- NULL
for (i in 1:floor(isa/100)) { 
  gskew <- abind(gskew, rcoskew(gfit.sk, from=1+(i-1)*100, to=i*100, standardize=T)[1:n,1:(n^2),])
}
gskew <- abind(gskew, rcoskew(gfit.sk, from=1+i*100, to=isa, standardize=T)[1:n,1:(n^2),])
rm(i)

gkurt <- array(dim = c(n,n**3,isa))
for (i in 1:floor(isa/100)) { 
  gkurt[,,(1+(i-1)*100):(i*100)] <- rcokurt(gfit.sk, from=1+(i-1)*100, to=i*100)[1:n,1:(n^2),]
}
# gkurt <- abind(gkurt, rcokurt(gfit.sk, from=1+i*100, to=isa, standardize=T)[1:n,1:(n^3),])
gkurt[,,(1+i*100):isa] <- rcokurt(gfit.sk, from=1+i*100, to=isa, standardize=T)
rm(i)
# gkurt <- abind(rcokurt(gfit.sk, from=1, to=100, standardize=T)[1:n,1:(n^3),], gkurt)

ret_mom.sk = fmoments(Data = ret5, n.ahead = 1, roll = oos, spec = mgspec.sk)#, 
                      # solver = "hybrid", cluster = NULL, gfun = "tanh", maxiter1 = 100000, 
                      # rseed = 500) 
# output variables
w.go.sk <- array(NA, dim = c(ts, n, length(ra))) # weights, for both in- and out-of-sample
mv.go.sk <- array(NA, dim = c(ts, 6, length(ra))) # mean, sd, utility of portfolio; parmastatus
dimnames(mv.go.sk) <- list(c(rep('isa', isa), rep('oos', oos)), 
                        c('portf.mean', 'portf.sd', 'portf.sk', 'portf.kurt', 'util', 'status'),
                        ra)

for (i1 in 1:length(ra)) {
  putil.go.sk <- vector(mode='list', length = ts)
  for(i in 1:isa){
    putil.go.sk[[i]] = parmautility(U = "CARA", method = "moment", 
                                    M1 = gmean.sk[i,], M2 =  gcov.sk[,,i], 
                                    M3 = gskew[,,i], M4 = gkurt[,,i],
                                    RA = ra[i1], budget = bu, 
                                    LB = LB, UB = UB)
    w.go.sk[i,,i1] = weights(putil.go.sk[[i]])
    mv.go.sk[i,'portf.mean',i1] <- parmareward(putil.go.sk[[i]])
    mv.go.sk[i,'portf.sd',i1] <- sqrt(w.go.sk[i,,i1] %*% gcov.sk[,,i] %*% w.go.sk[i,,i1])
    mv.go.sk[i,'portf.sk',i1] <- w.go.sk[i,,i1] %*% gskew[,,i] %*% (w.go.sk[i,,i1] %x% w.go.sk[i,,i1])
    mv.go.sk[i,'portf.kurt',i1] <- w.go.sk[i,,i1] %*% gkurt[,,i] %*% (w.go.sk[i,,i1] %x% (w.go.sk[i,,i1] %x% w.go.sk[i,,i1]))
  }
  rm(i)
  # utility with unnormalized moments
  mv.go.sk[1:isa,'util',i1] <- cara4(pm1 = mv.go.sk[1:isa,'portf.mean',i1], pm2 = mv.go.sk[1:isa,'portf.sd',i1]^2,
                                          pm3 = mv.go.sk[1:isa,'portf.sk',i1], l=ra[i1]) # pm4 = mv.go.sk[1:isa,'portf.kurt',i1],
    
  # normalized moments
#   mv.go.sk[1:isa,'portf.sk',i1] <- mv.go.sk[1:isa,'portf.sk',i1]/(mv.go.sk[1:isa,'portf.sd',i1]^3)
#   mv.go.sk[1:isa,'portf.kurt',i1] <- mv.go.sk[1:isa,'portf.kurt',i1]/(mv.go.sk[1:isa,'portf.sd',i1]^4)
  
  
# CARA, 3 Moment Approximation, out-of-sample
  for(i in (isa+1):ts){
    putil.go.sk[[i]] = parmautility(U = "CARA", method = "moment", 
                                    M1 = ret_mom.sk@moments$forecastMu[i-isa,], 
                                    M2 =  ret_mom.sk@moments$forecastCov[,,i-isa], 
                                    M3 = ret_mom.sk@moments$forecastM3[,,i-isa], 
                                    M4 = ret_mom.sk@moments$forecastM4[,,i-isa], 
                                    RA = ra[i1], budget = bu, 
                                    LB = LB, UB = UB)
    w.go.sk[i,,i1] = weights(putil.go.sk[[i]])
    mv.go.sk[i,'portf.mean',i1] <- parmareward(putil.go.sk[[i]])
    mv.go.sk[i,'portf.sd',i1] <- sqrt(w.go.sk[i,,i1] %*% ret_mom.sk@moments$forecastCov[,,i-isa] %*% w.go.sk[i,,i1])
    mv.go.sk[i,'portf.sk',i1] <- w.go.sk[i,,i1] %*% ret_mom.sk@moments$forecastM3[,,i-isa] %*% (w.go.sk[i,,i1] %x% w.go.sk[i,,i1])
    mv.go.sk[i,'portf.kurt',i1] <- w.go.sk[i,,i1] %*% ret_mom.sk@moments$forecastM4[,,i-isa] %*% (w.go.sk[i,,i1] %x% (w.go.sk[i,,i1] %x% w.go.sk[i,,i1]))
  }
  rm(i)
  # utility with unnormalized moments
  mv.go.sk[(isa+1):ts,'util',i1] <- cara4(pm1 = mv.go.sk[(isa+1):ts,'portf.mean',i1], pm2 = mv.go.sk[(isa+1):ts,'portf.sd',i1]^2,
                                          pm3 = mv.go.sk[(isa+1):ts,'portf.sk',i1], l=ra[i1]) # pm4 = mv.go.sk[(isa+1):ts,'portf.kurt',i1],
  
  # normalized moments
  # mv.go.sk[(isa+1):ts,'portf.sk',i1] <- mv.go.sk[(isa+1):ts,'portf.sk',i1]/(mv.go.sk[(isa+1):ts,'portf.sd',i1]^3)
  # mv.go.sk[(isa+1):ts,'portf.kurt',i1] <- mv.go.sk[(isa+1):ts,'portf.kurt',i1]/(mv.go.sk[(isa+1):ts,'portf.sd',i1]^4)
  mv.go.sk[,'status',i1] <- unlist(lapply(putil.go.sk, parmastatus))
}
rm(i1, LB, UB, putil.go.sk)

plot(mv.go.sk[,'portf.mean',5], type='l')
plot(as.numeric(mv.go.sk[,'status',5]), type='l')
table(as.numeric(mv.go.sk[,'status',]))/sum(table(as.numeric(mv.go.sk[,'status',])))

# forecast measures ####
# forecast accuracy measures
# here for univariate
# what about multivariate?
sqrt(mean((ret_mom@moments$forecastMu[1:oos,1] - ret6[(isa+1):ts,1])^2))
sqrt(mean((ret_mom.sk@moments$forecastMu[1:oos,1] - ret6[(isa+1):ts,1])^2))

mean(abs(ret_mom@moments$forecastMu[1:oos,1] - ret6[(isa+1):ts,1]))
mean(abs(ret_mom.sk@moments$forecastMu[1:oos,1] - ret6[(isa+1):ts,1]))

# mean(abs(100*(ret_mom@moments$forecastMu[1:oos,1] - ret6[(isa+1):ts,1])/ret6[(isa+1):ts,1]))
# mean(abs(100*(ret_mom.sk@moments$forecastMu[1:oos,] - ret6[(isa+1):ts,])/ret6[(isa+1):ts,]))

# risk measures ####
tech <- 0 # obs for risk measures calculations
riskcomp <- riskcomp.sk <- array(NA, c(ts-tech, 5, length(ra)))
dimnames(riskcomp) <- dimnames(riskcomp.sk) <- list(c(rep('isa', ts-tech-oos), rep('oos', oos)),
                           c('mad', 'sd', "cvar", "cdar", 'minimax'), ra)
for (i1 in 1:length(ra)) {
  riskcomp[,,i1] <- riskm(r = ret5[(tech+1):ts,], w = w.go[(tech+1):ts,,which(ra[i1]==ra)], oos = oos)
  riskcomp.sk[,,i1] <- riskm(r = ret5[(tech+1):ts,], w = w.go.sk[(tech+1):ts,,which(ra[i1]==ra)], oos = oos)
}
rm(i1)

sk2mom <- matrix(NA, ts, length(ra))
for (i1 in 1:length(ra)) {
    sk2mom[,i1] <- skewpt(ret, w.go[,,which(ra[i1]==ra)]) # portfolio skewness in 2 mom optimization
}
rm(i1)

# for review ##############
e1 <- new.env()
load(file = 'mvs3.RData', envir = e1)
ls(e1)
riskcomp <- get('riskcomp', envir = e1)

dimnames(riskcomp)[[2]]


# in-sample, no good picture
plot(riskcomp[,1,1], type='l', col='green',  ylab='', ylim=c(0,2)) # ,
par(new=T)
plot(riskcomp.sk[,1,1], type='l',  col='blue',  ylab='', ylim=c(0,2)) # ylim=c(0,0.2),
par(new=F)

hist(riskcomp[61:160,3,1],  main='',xlim=c(0.2,0.6),ylim=c(0,50), xlab='', col=rgb(red = 0, green = 1, blue = 0, alpha = 0.5)) #  
abline(v = median(riskcomp[61:160,3,1]), lty=2, col='green')
par(new=T)
hist(riskcomp.sk[61:160,3,1], main='',xlim=c(0.2,0.6),ylim=c(0,50), xlab='', col=rgb(red = 0, green = 0, blue = 1, alpha = 0.5)) # 
abline(v = median(riskcomp.sk[61:160,3,1]), lty=2, col='blue')
par(new=F)

# mad, good picture for lambda=c(1,2)
plot(riskcomp.sk[51:56,1,2], type='l', ylim=c(0,0.2), ylab='') 
par(new=T)
plot(riskcomp[51:56,1,2], type='b', ylim=c(0,0.2), ylab='')
par(new=F)

# sd, good picture for lambda=c(1,2)
plot(riskcomp.sk[51:56,2,4], type='l', ylim=c(0,0.1), ylab='')
par(new=T)
plot(riskcomp[51:56,2,4], type='b', ylim=c(0,0.1), ylab='')
par(new=F)

# cvar, good picture for lambda=c(1,2)
plot(riskcomp.sk[51:56,3,2], type='l', ylim=c(0,0.28), ylab='')
par(new=T)
plot(riskcomp[51:56,3,2], type='b', ylim=c(0,0.28), ylab='')
par(new=F)

# cdar, good picture for lambda=c(1,2)
plot(riskcomp.sk[51:56,4,1], type='l', ylim=c(0,3), ylab='')
par(new=T)
plot(riskcomp[51:56,4,1], type='b', ylim=c(0,3), ylab='')
par(new=F)

# average weights
tab <- matrix(NA, n*2, length(ra))
for (i1 in 1:length(ra)) {
  # tab[1:n,i1] <- sprintf('%0.3f', colMeans(w.go[,,which(ra[i1]==ra)]))
  tab[(n+1):(2*n),i1] <- sprintf('%0.3f', colMeans(w.go.sk[,,which(ra[i1]==ra)]))
}
write.table(tab, 'clipboard', quote = F)
###################

save.image(file = 'mvs13.RData')
save(list='da_yah', file='mvs8.RData')

# indifference map
# View(colMeans(w.go))
# write.table(colMeans(w.go), file='clipboard', sep='\t', dec=',')
