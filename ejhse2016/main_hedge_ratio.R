setwd("F:\\Dropbox\\Научка\\hedge ratio\\r_hedge") # Megatron
# setwd("C:/Users/Одмин/Dropbox/Научка/hedge ratio/r_hedge/") # home station
getwd()

require(rmgarch)

names <- c("gazp", "gmkn", "lkoh", "nvtk", "rosn", "sber", "sngs", "trnf", "vtbr")

source("aprilConf.R")

dataFromTxt <- readAll(names)
mergedData <- mergeData(dataFromTxt)
prices <- accountLotSize(mergedData)
prices <- accountSberSplit(prices)
rm("dataFromTxt", "mergedData")

logret <- calcLogRet(prices)

ds.prices <- summStatAll(prices) 
ds.logret <- summStatAll(logret)

# to copy ds to Excel
write.table(ds.prices[[1]], "clipboard", sep ="\t", dec = ",") # stocks
write.table(ds.prices[[2]], "clipboard", sep ="\t", dec = ",") # futures
write.table(ds.logret[[1]], "clipboard", sep ="\t", dec = ",") # stocks
write.table(ds.logret[[2]], "clipboard", sep ="\t", dec = ",") # futures

save.image("aprilConf_start.RData")

# createBoxplots(prices, F) # not informative due to the means differ very much
createBoxplots(logret, T) # plots are to be saved manually

# require(rmgarch)
# require(sandwich)

# make constant correlation test for 1st and 2nd lags
constCorr <- testConstCorr(logret)
write.table(constCorr, "clipboard", sep ="\t", dec = ",")

ggResult <- estimGoGarch(logret)
adccResult <- estimAdcc(logret)
copResult <- estimCopGarch(logret) 
save.image("aprilConf_estimRes2.RData")

ll <- extractLL(dcc = adccResult, cop = copResult, gg = ggResult)
write.table(ll, "clipboard", sep ="\t", dec = ",")
apply(ll, 1, function(x) {x==max(x)}) # to see where LL is max

bic <- extractBic(dcc = adccResult, cop = copResult, gg = ggResult)
write.table(bic, "clipboard", sep ="\t", dec = ",")
apply(bic, 1, function(x) {x==min(x)})

vc <- ExtractVarCov(dcc = adccResult, cop = copResult, gg = ggResult, forecast = F)
hr <- calcHedgeRat(vc)
hret <- CalcHedgedRet(logret, hr, CalcConstHedgeRat(logret))
eff <- calcEff(logret, hret)

save.image('aprilConf_eff.RData')

# plot effectiveness in EXCEL!!!! SATANA!!!
write.table(eff, "clipboard", sep ="\t", dec = ",")
apply(eff, 2, function(x) {x==max(x)}) # to see which model is effective for each portfolio

profit <- apply(hret[[6]],2,sum) # profit as sum of returns
write.table(profit, "clipboard", sep ="\t", dec = ",")

PlotKernelDens(logret, hret, 5, 'kernel_rosn') # kernel estimates of returns

# Diebold-Mariano test
losses2 <- CalcLossFun(dccF = adccResult[[2]], copF = copResult[[2]], 
                       ggF = ggResult[[2]], lr = logret) # calculate loss function, type 2
diffs2 <- CalcLossDiffs(losses2)
require(sandwich)
dm2Mean <- CalcDiebMarTest(diffs2, F)
dm2Med <- CalcDiebMarTest(diffs2, T)
dm2MeanProb <- CalcProb(dm2Mean)
dm2MedProb <- CalcProb(dm2Med)
dm2MeanProb<0.05 # show significant results
dm2MedProb<0.05 # show significant results
# results are almost equal
write.table(dm2Mean, "clipboard", sep ="\t", dec = ",")
write.table(dm2MeanProb[,3], "clipboard", sep ="\t", dec = ",")


vcForc <- ExtractVarCov(dcc = adccResult, cop = copResult, gg = ggResult, forecast = TRUE)
hrForc <- calcHedgeRat(vcForc)
hretForc <- CalcHedgedRet(logret, hrForc, CalcConstHedgeRat(logret))
effForc <- calcEff(logret, hretForc)

### further there is SRACH!!! >:|
### SRACH is everywhere!!!!

# try to estimate relative risk aversion
# see Cotter, Hanly, A utility based approach to energy hedging, 2012
micex <- readTxt('MICEXINDEXCF')
require(xts)
micret <- xts(x = diff(log(micex[,2])), order.by = micex[-1,1])
mic.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 0),
                                             submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
                       mean.model = list(armaOrder = c(1, 0), include.mean = TRUE, archm = TRUE,
                                         archpow = 2, arfima = FALSE, external.regressors = NULL, archex = FALSE),
                       distribution.model = "sstd")
mic.ufit <- ugarchfit(mic.spec, micret, out.sample = 0, solver = "nloptr")
mic.ufit@fit$robust.matcoef
infocriteria(mic.ufit)
write.table(mic.ufit@fit$robust.matcoef, "clipboard", sep ="\t", dec = ",")
risk.av <- mic.ufit@fit$robust.matcoef['archm',1]

futRetFitted <- ExtractFittedFutRet(dcc = adccResult, cop = copResult, gg = ggResult)
hrMeanVarForc <- calcHedgeRatMeanVar(vcForc, futRetFitted, risk.av) # 1e20
hrMeanVarForcAver <- lapply(hrMeanVarForc, colMeans)
hrMeanVarForcAver <- matrix(unlist(hrMeanVarForcAver),3,9)
dimnames(hrMeanVarForcAver) <- list(dimnames(hrMeanVarForc[[1]])[[2]], names)
write.table(hrMeanVarForcAver, "clipboard", sep ="\t", dec = ",")

hretMeanVarForc <- CalcHedgedRet(logret, hrMeanVarForc, CalcConstHedgeRat(logret))

hist(hretMeanVarForc[['sngs']][,3])
write.table(hretMeanVarForc[[2]][,2:3],"clipboard", sep ="\t", dec = ",")
# calculation of skewness
require(moments)
hret.skew <- skewness(hretMeanVarForc[[2]][,2])
n.hret <- length(hretMeanVarForc[[2]][,2])
hret.skew.signif <- hret.skew*sqrt((6*n.hret)*(n.hret-1)/(n.hret-2)/(n.hret+1)/(n.hret+3))
hret.skew.signif # critical vale is |2|

hret.ds <- lapply(hretMeanVarForc,var)
hret.ds <- matrix(unlist(hret.ds),4,9)
write.table(hret.ds,"clipboard", sep ="\t", dec = ",")

write.table(colMeans(hretMeanVarForc[[2]]),"clipboard", sep ="\t", dec = ",")
write.table(apply(hretMeanVarForc[[2]],2,var),"clipboard", sep ="\t", dec = ",")
write.table(effMeanVarForc,"clipboard", sep ="\t", dec = ",")

# calculation of VaR
# later; needed GARCH for hret

effMeanVarForc <- calcEff(logret, hretMeanVarForc)
profitMeanVarForc <- CalcProfit(logret, hretForc) # profit as sum of returns
write.table(profitMeanVarForc, "clipboard", sep ="\t", dec = ",")

write.table(effForc, "clipboard", sep ="\t", dec = ",")
apply(effForc, 2, function(x) {x==max(x)}) # to see which model is effective for each portfolio

profitForc <- apply(hretForc[[8]],2,sum) # profit as sum of returns
write.table(profitForc, "clipboard", sep ="\t", dec = ",")

PlotKernelDens(logret, hretForc, 5, 'kernel_rosn_forc') # kernel estimates of returns

# extract coefs with se
ggMatcoefGmkn <- rbind(cbind(ggResult[[1]][[2]]@mfit$ufit@fit[[1]]@fit$robust.matcoef[,1],
                             ggResult[[1]][[2]]@mfit$ufit@fit[[1]]@fit$robust.matcoef[,4]),
                       cbind(ggResult[[1]][[2]]@mfit$ufit@fit[[2]]@fit$robust.matcoef[,1],
                             ggResult[[1]][[2]]@mfit$ufit@fit[[2]]@fit$robust.matcoef[,4]))
dimnames(ggMatcoefGmkn)[[2]] <- c('gg.coef', 'gg.prob')
copMatcoefGmkn <- cbind(adccResult[[1]][[2]]@mfit$matcoef[,1],
                              adccResult[[1]][[2]]@mfit$matcoef[,4])                        
dimnames(copMatcoefGmkn)[[2]] <- c('cop.coef', 'cop.prob')
write.table(ggMatcoefGmkn, "clipboard", sep ="\t", dec = ",")
write.table(copMatcoefGmkn, "clipboard", sep ="\t", dec = ",")
likelihood(ggResult[[1]][[2]]); likelihood(copResult[[1]][[2]])

View(printGgCoefsSe())

n.as <- 4
mse.fut(ggResult[[1]][[n.as]],logret[[n.as]])
mse.fut(copResult[[1]][[n.as]],logret[[n.as]])
rm(n.as) # gg is worse :(((

source("aprilConf.R")

save.image('aprilConf_estimRes4.RData')

# todo: 
# - plot const and non-const hedge ratios
# - why ols and mgarch give the same results?
# - eliminate SRACH!!!!!
 
# setwd("F:\\Dropbox\\Научка\\hedge ratio\\r_hedge")
load('aprilConf_estimRes4.RData')

