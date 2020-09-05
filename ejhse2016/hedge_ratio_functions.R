readTxt <- function(ticker) {
  # read data from *.txt got from Finam.ru (FEFO)
  fileName <- paste0("F:/Dropbox/Научка/hedge ratio/data/для Апрельской/", ticker, "_070101_141001.txt")
  readedData <- read.table(fileName, check.names = FALSE, header=F,
                           skip=1, sep=",", dec=".", 
                           colClasses = c("character", rep("NULL", 4), "double", "NULL"))
  names(readedData) <- c("date", "close")
  readedData$date <- as.Date(readedData$date, format="%d%m%y")
  return(readedData)
}

readAll <- function(names) {
  # read data from 9 files in cycle
  l <- length(names)
  readyData <- vector(mode = "list", length = 2*l)  # vector for the result
  tickerAll <- NULL # vector of ticker names
  for (i1 in 1:(l*2)) {
    if (i1 < l+1) {
      ticker <- toupper(names[i1]) # ticker for stock
    }
    else {
      ticker <- paste0("SPFB.", toupper(names[i1-9])) # ticker for futures
    }
    tickerAll <- c(tickerAll, ticker)
    readedData <- readTxt(ticker)

    readyData[[i1]] <- readedData
  }
  names(readyData) <- tickerAll
  readyData # explicit return() is a little bit slower
}

# merge stock and futeres prices
mergeData <- function(data) {
  l <- length(data)/2
  readyData <- vector(mode = "list", length = l)  # vector for the result
  names(readyData) <- tolower(head(names(data),l))
  for (i1 in 1:l) {
#     browser()
    mergedData <- merge(data[[i1]], data[[i1+9]], by=1, 
                        suffixes = c(".stock",".fut"))
    readyData[[i1]] <- mergedData
  }
  readyData
}

accountSberSplit <- function(prices) {
  splitDate <- as.Date(x = "18/07/07", format="%d/%m/%y")
  ind <- prices$sber$date<splitDate
  prices$sber[ind, "close.stock"] <- prices$sber[ind, "close.stock"]/1e3
  prices
}

accountLotSize <- function(data) {
  lotSize <- c(100, 10, 10, 100, 100, 1, 1000, 1, 1e5) # lot size for futures
  for (i1 in 1:length(data)) {
    if(names(data)[i1]=="sber") {
      changeLotSizeDate <- as.Date("10/12/07", format="%d/%m/%y")
      ind <- prices$sber$date<=changeLotSizeDate
      data[[i1]][ind, "close.fut"] <- data[[i1]][ind, "close.fut"]/1e3
      data[[i1]][!ind, "close.fut"] <- data[[i1]][!ind, "close.fut"]/1e2
    }
    else {
      data[[i1]]$close.fut <- data[[i1]]$close.fut/lotSize[i1] 
    }
  }
  data
}

summStat <- function(data.vec) {
  ds.names <- c("N", "Min.", "1stQ.", "Mean", "Median", "3rdQ.", 
                "Max.", "St.dev.", "Skewn.", "Kurt.")
  ds <- matrix(nrow = length(ds.names), ncol = 1,
               dimnames=list(ds.names, "stat"))
  ds["N",] <- length(data.vec)
  ds["Min.",] <- min(data.vec)
  ds["1stQ.",] <- quantile(x = data.vec, probs=0.25)
  ds["Mean",] <- mean(data.vec)
  ds["Median",] <- median(data.vec)
  ds["3rdQ.",] <- quantile(x = data.vec, probs=0.75)
  ds["Max.",] <- max(data.vec)
  ds["St.dev.",] <- sd(data.vec)
  ds["Skewn.",] <- skewness(data.vec) # Gauss has skewness=0
  ds["Kurt.",] <- kurtosis(data.vec) # Gauss has kurtosis=3
  ds
}

summStatAll <- function(data) {
  l <- length(data)
  ds.names <- c("N", "Min.", "1stQ.", "Mean", "Median", "3rdQ.", 
                "Max.", "St.dev.", "Skewn.", "Kurt.") # range?
  ds.stock <- matrix(nrow = length(ds.names), ncol = l, 
                dimnames = list(ds.names, names(data)))
  ds.fut <- ds.stock
  library(moments)
  for (i1 in 1:l) {
    ds.stock[,i1] <- summStat(data[[i1]][,"close.stock"])
    ds.fut[,i1] <- summStat(data[[i1]][,"close.fut"])
  }
  detach("package:moments")
  ds <- list(ds.stock = ds.stock, ds.fut = ds.fut)
}

calcLogRet <- function(prices) {
  l <- length(prices)
  lr <- vector("list", l)
  names(lr) <- names(prices)
  for (i1 in 1:l) {
    lr1 <- prices[[i1]][-1,]
    names(lr1) <- names(prices[[i1]])
    lr1[,1] <- prices[[i1]][-1,"date"]
    lr1[,2] <- diff(log(prices[[i1]][,"close.stock"]))
    lr1[,3] <- diff(log(prices[[i1]][,"close.fut"]))
    lr[[i1]] <- lr1
  }
  lr
}

createBoxplots <- function(data, logret) {
  l <- length(data)
  lengths <- sapply(data,function(x) {nrow(x)})
  dataBoxplot.stock <- matrix(nrow=max(lengths), ncol=l)
  dataBoxplot.fut <- dataBoxplot.stock
  for (i1 in 1:l) {
    dataBoxplot.stock[1:lengths[i1],i1] <- data[[i1]][,"close.stock"]
    dataBoxplot.fut[1:lengths[i1],i1] <- data[[i1]][,"close.fut"]                                                
  }
  boxplot.matrix(dataBoxplot.stock, axes=FALSE)
  axis(1, at=1:l, labels=names(data))
  axis(2)
  title(main="Акции", xlab="Тикер", ylab=ifelse(logret,"Лог. доходность","Цена"))
  boxplot.matrix(dataBoxplot.fut, axes=FALSE)
  axis(1, at=1:l, labels=names(data))
  axis(2)
  title(main="Фьючерсы", xlab="Тикер", ylab=ifelse(logret,"Лог. доходность","Цена"))
}

testConstCorr <- function(ret) {
  # Engle, Sheppard, 2001. 
  # ret - log returns
  # function does 1 and 2 lag testing
  l <- length(ret)
  res<-matrix(0,l,2)
  for (i1 in 1:l) {
    test<-DCCtest(ret[[i1]][,2:3], garchOrder = c(1,1), n.lags = c(1,2))
    res[i1,]<-test$p.value # CHECK THIS IN SOURCES!!!!
  }
  res
}

estimGoGarch <- function(logret) {
  fit <- list()
  forc <- list()
  l <- length(logret)
  ggSpec <- gogarchspec(mean.model = list(model = "AR", lag.max=3, lag.criterion="HQ"),
                       variance.model = list(model = "fGARCH", 
                                             garchOrder = c(1,1), 
                                             submodel = "GARCH"),
                       distribution.model = "magh")
  for (i1 in 1:l) {
    os <- floor(nrow(logret[[i1]])/3)
    ggFit <- gogarchfit(ggSpec,  logret[[i1]][,2:3], out.sample = os)
    ggForc <- gogarchforecast(ggFit, n.ahead = 1, n.roll = os)
    fit <- c(fit, ggFit)
    forc <- c(forc, ggForc)
  }
  names(fit) <- names(logret)[1:length(fit)]
  names(forc) <- names(logret)[1:length(forc)]
  res <- list(fit=fit, forc=forc)
}

estimAdcc <- function(logret) {
  fit <- list()
  forc <- list()
  l <- length(logret)
  #   groups <- c(1, 3, 1, 1, 1, 2, 1, 1, 2) # 1 - oil&gas, 2 - finance, 3 - metallurgy
#   groupsPerm <- c(1, 3, 4, 5, 7, 8, 6, 9, 2)
#   logret <- logret[groupsPerm]
  
  uniGarchSpec <- ugarchspec(mean.model=list(armaOrder=c(1,0)),
                             variance.model=list(model="fGARCH",submodel="GARCH",
                                                 garchOrder=c(1,1)))
  multiGarchSpec <- multispec(replicate(2, uniGarchSpec))
  dccSpec <- dccspec(multiGarchSpec,dccOrder=c(1,1), model="aDCC", 
                     distribution="mvt")
  for (i1 in 1:l) {
    os <- floor(nrow(logret[[i1]])/3)
    adccFit <- dccfit(dccSpec, logret[[i1]][,2:3], out.sample = os)
    adccForc <- dccforecast(adccFit, n.ahead=1,n.roll=os)
    fit <- c(fit, adccFit)
    forc <- c(forc, adccForc)
  }
  names(fit) <- names(logret)[1:length(fit)]
  names(forc) <- names(logret)[1:length(forc)]
  res <- list(fit=fit, forc=forc)
}

estimCopGarch <- function(logret) {
  fit <- list()
  forc <- list()
  l <- length(logret)
  uniGarchSpec <- ugarchspec(mean.model=list(armaOrder=c(1,0)),
                             variance.model=list(model="fGARCH",submodel="GARCH",
                                                 garchOrder=c(1,1)))
  multiGarchSpec <- multispec(replicate(2, uniGarchSpec))
  copSpec <- cgarchspec(multiGarchSpec,
                        dccOrder = c(1, 1), asymmetric = T,
                        distribution.model = list(copula = "mvt",
                                                  method = "ML", 
                                                  time.varying = T,
                                                  transformation = "empirical")) # parametric
  for (i1 in 1:l) {
    pair <- logret[[i1]][,2:3]
    os <- floor(nrow(pair)/3)
    copFit <- cgarchfit(copSpec, pair, out.sample = os, fit.control = list(eval.se=FALSE))
    copForc <- copGarch1AheadForecast(copFit, os, copSpec)

    fit <- c(fit, copFit)
    forc <- c(forc, copForc)
    print(i1)
  }
  names(fit) <- names(logret)[1:length(fit)]
  names(forc) <- names(logret)[1:length(forc)]
  res <- list(fit=fit, forc=forc)
}

copGarch1AheadForecast <- function(fit3, os, copSpec) {
  # fit3 - copFit
  logret1 <- fit3@model$modeldata$data
  T = dim(logret1)[1]-os
  nc <- dim(logret1)[2]
  simMu = simS = filtMu = filtS = matrix(NA, ncol = nc, nrow = os)
  simCor = simC = filtC = filtCor = array(NA, dim = c(nc,nc,os))
  colSd = function(x) apply(x, 2, "sd")
  specx = copSpec
  for(i in 1:nc) specx@umodel$fixed.pars[[i]] = as.list(fit3@model$mpars[fit3@model$midx[,i]==1,i])
  setfixed(specx)<-as.list(fit3@model$mpars[fit3@model$midx[,nc+1]==1,nc+1])
  
  for(i in 1:os){
    if(i==1){
      # arma = c(1,0) therefore need 1 lag
      lag <- 1
      presigma = matrix(tail(sigma(fit3), lag), ncol = nc)
      prereturns = matrix(unlist(logret1[(T-1):T, ]), ncol = nc, nrow = lag)
      preresiduals = matrix(tail(residuals(fit3),lag), ncol = nc, nrow = lag)
      preR = last(rcor(fit3))[,,1]
      diag(preR) = 1
      preQ = fit3@mfit$Qt[[length(fit3@mfit$Qt)]]
      preZ = tail(fit3@mfit$Z, 1)
      
      tmp = cgarchfilter(specx, logret1[1:(T+1), ], filter.control = list(n.old = T))
      filtMu[i,] = tail(fitted(tmp), 1)
      filtS[i,] = tail(sigma(tmp), 1)
      filtC[,,i] = last(rcov(tmp))[,,1]
      filtCor[,,i] = last(rcor(tmp))[,,1]
    } 
    else {
      presigma = matrix(tail(sigma(tmp), lag), ncol = nc)
      prereturns = matrix(unlist(logret1[(T+i-2):(T+i-1), ]), ncol = nc, nrow = lag)
      preresiduals = matrix(tail(residuals(tmp),lag), ncol = nc, nrow = lag)
      preR = last(rcor(tmp))[,,1]
      diag(preR) = 1
      preQ = tmp@mfilter$Qt[[length(tmp@mfilter$Qt)]]
      preZ = tail(tmp@mfilter$Z, 1)
      
      tmp = cgarchfilter(specx, logret1[1:(T+i), ], filter.control = list(n.old = T))			
      filtMu[i,] = tail(fitted(tmp), 1)
      filtS[i,] = tail(sigma(tmp), 1)
      filtC[,,i] = last(rcov(tmp))[,,1]
      filtCor[,,i] = last(rcor(tmp))[,,1]
    }
    sim3 = cgarchsim(fit3, n.sim = 1, m.sim = 2000, startMethod = "sample", preR = preR, preQ = preQ, preZ = preZ,
                     prereturns = prereturns, presigma = presigma, preresiduals = preresiduals)
    simx = t(sapply(sim3@msim$simX, FUN = function(x) x[1,]))
    simMu[i,] = colMeans(simx) # conditional mean
    # Note: There is no uncertainty for the 1-ahead simulation of cov
    simC[,,i] = sim3@msim$simH[[1]][,,1] # conditional covariance
    simS[i,] = sqrt(diag(simC[,,i])) # conditional sigma
    simCor[,,i] = sim3@msim$simR[[1]][,,1] # conditional correlation
#     print(i) 
    simRes <- list(simMu = simMu, simS = simS, simC = simC, simCor = simCor)
  }
  simRes
}

extractLL <- function(dcc, cop, gg) {
  # dcc, cop, gg  - lists with estimation results (both fit and forecast)
  l <- length(dcc[[1]])
  ll.dcc = ll.cop = ll.gg <- NULL
  for (i1 in 1:l) {
    ll.dcc <- rbind(ll.dcc, likelihood(dcc[[1]][[i1]]))
    ll.cop <- rbind(ll.cop, likelihood(cop[[1]][[i1]]))
    ll.gg <- rbind(ll.gg,likelihood(gg[[1]][[i1]]))
  }
  ll <- cbind(ll.dcc, ll.cop, ll.gg)
  dimnames(ll)[[1]] <- names(dcc[[1]])
  dimnames(ll)[[2]] <- c("dcc", "cop", "gg")
  ll
}

# here BIC is nObs-standardized, reference?? find in rugarch manuals
extractBic <- function(dcc, cop, gg) {
  # dcc, cop, gg  - lists with estimation results (both fit and forecast)
  l <- length(dcc[[1]])
  bic.dcc =  bic.cop = bic.gg <- NULL
  for (i1 in 1:l) {
#     browser()
    bic.dcc <- rbind(bic.dcc, infocriteria(dcc[[1]][[i1]])[2])
    bic.cop <- rbind(bic.cop, calcBic(cop[[1]][[i1]]))
    bic.gg <- rbind(bic.gg, calcBic(gg[[1]][[i1]]))
  }
  bic <- cbind(bic.dcc, bic.cop, bic.gg)
  dimnames(bic)[[1]] <- names(dcc[[1]])
  dimnames(bic)[[2]] <- c("dcc", "cop", "gg")
  bic
}

calcBic <- function(obj) {
  # obj - rmgarch object
  nobs <- nrow(fitted(obj))
  npars <- length(coef(obj))
  ll <- likelihood(obj)
  bic1 <- (-2*ll)/nobs + npars * log(nobs)/nobs
  bic1
}

ExtractVarCov <- function(dcc, cop, gg, forecast) { #todo: here it should be "..."
  # forecast - logical; if TRUE extract forecasted volatility
  l <- length(dcc[[1]])
  vc <- vector(mode="list", length = l)
  names(vc) <- names(dcc[[1]])
  if(!forecast) {
  models <- list(dcc[[1]], cop[[1]], gg[[1]]) # fit objects
  }
  else {
    models <- list(dcc[[2]], cop[[2]], gg[[2]]) # forecast objects
  }
  for (i1 in 1:l) {
    vc1p <- vector("list", 3)
    for (i2 in 1:3) {
      # 3 - number of models
      vc1p[[i2]] <- VarCov1Portf1Model(models[[i2]], i1, forecast)
      names(vc1p) <- c("dcc", "cop", "gg")
    }
    vc[[i1]] <- vc1p
  }
  vc
}

VarCov1Portf1Model <- function(model, stockInd, forecast) {
  # variance-covariance for one portfolio and one model
  # model - one of the dcc, cop or gg
  # stockInd - index of stock
  # forecast - logical; if TRUE extract forecasted volatility
  if(forecast && typeof(model[[stockInd]])=='double') {
    # logical expression means that model is copula garch, which has no rcov() method
    var.stock <- model[[stockInd-1+3*stockInd]][1,1,] # stock variance
    var.fut <- model[[stockInd-1+3*stockInd]][2,2,] # futures variance
    cov.stock.fut <- model[[stockInd-1+3*stockInd]][1,2,] # stock-futures covariance
  }
  else {
    trim <- length(rcov(model[[stockInd]])) # trim the last element, because copula forecast has no T+1 forecast
    var.stock <- Convert2Array(rcov(model[[stockInd]]))[1,1,-trim] # stock variance
    var.fut <- Convert2Array(rcov(model[[stockInd]]))[2,2,-trim] # futures variance
    cov.stock.fut <- Convert2Array(rcov(model[[stockInd]]))[1,2,-trim] # stock-futures covariance
  }
  vc1p1m <- cbind(var.stock, var.fut, cov.stock.fut)
  names(vc1p1m) <- c("var.stock", "var.fut", "cov.stock.fut")
  vc1p1m
}

ExtractFittedFutRet <- function(dcc, cop, gg) {
  # extract conditional returns of hedging asset (futures)
  # dcc, cop, gg - estimated models (list of 2 each)
  
  l <- length(dcc[[1]])  
  futRetFit <- vector(mode='list', length = l)
  names(futRetFit) <- names(dcc[[1]])
  models <- list(dcc[[1]], cop[[1]], gg[[1]]) # fit objects
  for (i1 in 1:l) {
    nr <- nrow(fitted(models[[1]][[i1]]))
    frf1p <- matrix(0,nr,length(models))
    for (i2 in 1:length(models)) {
      frf1p[,i2] <- fitted(models[[i2]][[i1]])[,'close.fut'] # futures' returns only
      dimnames(frf1p)[[2]] <- c("dcc", "cop", "gg")
    }
    futRetFit[[i1]] <- frf1p
  }
  futRetFit
}

calcHedgeRat <- function(vc) {
  # calculate hedge ratio
  l <- length(vc)
  nmod <- length(vc[[1]]) # number of models
  hr <- vector("list", l)
  for (i1 in 1:l) {
    hr1p <- matrix(NA, nrow=dim(vc[[i1]][[1]])[1], ncol=nmod) # 1p - one portfolio
    for (i2 in 1:nmod) {
#       browser()
      hr1p[,i2] <- vc[[i1]][[i2]][,"cov.stock.fut"]/vc[[i1]][[i2]][,"var.fut"]
    }
    dimnames(hr1p)[[2]] <- names(vc[[1]])
    hr[[i1]] <- hr1p
  }
  names(hr) <- names(vc)
  hr
}

calcHedgeRatMeanVar <- function(vc, fr, rra) {
  # calculate hedge ratio taking into account both mean and variance
  # through 2-moment utility function
  # see Wahab, 1995
  
  # vc - covariance matrix of portfolio
  # fr - conditional returns of hedging asset (futures)
  # rra - relative risk aversion
  
  l <- length(vc)
  nmod <- length(vc[[1]]) # number of models
  hr <- vector("list", l)
  for (i1 in 1:l) {
    tos <- dim(vc[[i1]][[1]])[1] # length of out-of-sample
    hr1p <- matrix(NA, nrow=tos, ncol=nmod) # 1p - one portfolio
    for (i2 in 1:nmod) {
      #       browser()
      hr1p[,i2] <- (-tail(fr[[i1]][,i2],tos)/(2*rra)+vc[[i1]][[i2]][,"cov.stock.fut"])/vc[[i1]][[i2]][,"var.fut"]
    }
    dimnames(hr1p)[[2]] <- names(vc[[1]])
    hr[[i1]] <- hr1p
  }
  names(hr) <- names(vc)
  hr
}

CalcConstHedgeRat <- function(logret) {
  # calculate constant hedge ratios by OLS
  constHr <- sapply(logret,FUN = function(x) {cc <- cov(x[,2:3]); 
                                              return(cc[1,2]/cc[2,2])})
}

CalcHedgedRet <- function(lr, hr, cHr) {
  # calculate returns of hedged position
  # lr - log returns
  # hr - hedge ratios
  # cHr - constant hedge ratios
  
  l <- length(lr) 
  nmod <- dim(hr[[1]])[[2]]+1 # number of multivariate models and OLS
  hret <- vector("list", l)
  for (i1 in 1:l) {
    tis <- dim(hr[[i1]])[1] # T in-sample, which is less that length(lr)
    hr[[i1]] <- cbind(hr[[i1]], rep(cHr[i1], length.out = tis))
    dimnames(hr[[i1]])[[2]][nmod] <- "ols"
    hret1p <- matrix(NA, nrow=tis, ncol=nmod) # 1p - one portfolio
    for (i2 in 1:nmod) {
      hret1p[,i2] <- lr[[i1]][1:tis,2]-hr[[i1]][,i2]*lr[[i1]][1:tis,3]
    }
    dimnames(hret1p)[[2]] <- dimnames(hr[[1]])[[2]]
    hret[[i1]] <- hret1p
  }
  names(hret) <- names(hr)
  hret    
}

calcEff <- function(lr, hret) { 
  # calculate effectiveness of hedging (hedged vs. unhedged)
  # lr - log returns
  # hret - hedged returns
  
  l <- length(lr)
  nmod <- dim(hret[[1]])[2] # number of models
#   browser()
  eff <- matrix(NA, nrow=nmod, ncol=l)
  dimnames(eff)[[1]] <- dimnames(hret[[1]])[[2]]
  dimnames(eff)[[2]] <- names(lr)
  for (i1 in 1:l) {
    tis <- dim(hret[[i1]])[1] # T in-sample, which is less that length(lr)
    for (i2 in 1:nmod) {
      eff[i2, i1] <- 1-var(hret[[i1]][,i2])/var(lr[[i1]][1:tis,2])
    }
  }
  eff
}

CalcProfit <- function(lr, hret) {
  # calculate profit of hedged and unhedged positions
  # see Penikas, 2011
  l <- length(lr)
  prf <- lapply(hret, function(x) {apply(x,2,sum)})
  prf <- matrix(unlist(prf), ncol(hret[[1]]), l)
  prf.unhed <- lapply(lr, function(x) {apply(tail(x[,2:3], floor(nrow(x)/3)),2,sum)})
  prf.unhed <- matrix(unlist(prf.unhed), ncol(lr[[1]]), l)
  prf <- rbind(prf, prf.unhed[1,]) # bind only stocks' sum return
  dimnames(prf) <- list(c(dimnames(hret[[1]])[[2]],'unhedged'),names(hret))
  prf
}

PlotKernelDens <- function(logret, hret, l, filename) {
  # l - number or ticker of a stock for example
  # filename - name of file with plot
#   name <- names(logret)[[l]]
  png(filename = paste0(filename, '.png'), width = 600, height = 350)
  par(lwd = 1.5)
  plot(density(hret[[l]][,'gg']), lty=2, xlab="", ylab="", main="",
       xlim = c(-0.12, 0.12))
  lines(density(hret[[l]][,'dcc']), lty=3)
  lines(density(hret[[l]][,'cop']), lty=4)
  lines(density(hret[[l]][,'ols']), lty=5)
  lines(density(logret[[l]][,2]), lty=1, lwd=2)
  title(main = '', xlab = "Логарифмическая доходность", 
        ylab="Плотность распределения")
  legend(x='topleft', 
         legend = c('GO-GARCH', 'ADCC', 'copula GARCH', 'OLS', 'Без хеджирования'),
         lty=c(2,3,4,5,1), lwd=c(1.5,1.5,1.5,1.5,2), bty='n', cex=0.8, adj=0.007)
  dev.off()
}

CalcLossFun <- function(dccF, copF, ggF, lr) {
  l <- length(lr)
  os <- rep(NA, times = l)
  trueVol <- list()
  losses <- vector("list", l)
  names(losses) <- names(lr)
  volF <- ExtractVolForecast(dccF, copF, ggF)
  dccVolF <- volF[[1]]
  copVolF <- volF[[2]]
  ggVolF <- volF[[3]]
  for (i2 in 1:l) {
    os[i2] <- dim(copF[[i2+3*(i2-1)]])[1]
    trueVol[[i2]] <- CalcTrueVolProxy(lr[[i2]], os[i2])
#         browser()
    losses[[i2]] <- rbind(CalcLoss2(trueVol[[i2]], dccVolF[[i2]]),
                          CalcLoss2(trueVol[[i2]], copVolF[[i2]]),
                          CalcLoss2(trueVol[[i2]], ggVolF[[i2]]))
    losses[[i2]] <- t(losses[[i2]])
    dimnames(losses[[i2]])[[1]] <- names(dccF)
    dimnames(losses[[i2]])[[2]] <- c("adcc", "cop", "gg")    
  }
  losses
}

ExtractVolForecast <- function(dccF, copF, ggF) {
  l <- length(dccF)
  os <- rep(NA, times = l)
  dccVolF = copVolF = ggVolF = trueVol <- list()
  for (i2 in 1:l) {
    os[i2] <- dim(copF[[i2+3*(i2-1)]])[1]
    dccVolF. <- rcov(dccF[[i2]])
    dccVolF[[i2]] <- Convert2Array(dccVolF.)
    copVolF[[i2]]  <- copF[[i2-1+3*(i2)]]
    ggVolF. <- rcov(ggF[[i2]])
    ggVolF[[i2]] <- Convert2Array(ggVolF.)
  }
  volF <- list(adccVolF = dccVolF, copVolF = copVolF, ggVolF = ggVolF)
#   names(volF) <- c('', '', '')
  volF
}

Convert2Array <- function(lst) {
  # lst - list, returned by rcov() function
  arr <- array(unlist(lst), dim = c(nrow(lst[[1]]), ncol(lst[[1]]), length(lst)))
  arr
}

CalcTrueVolProxy <- function(lr1, os1) {
  # calculate true volatility proxy
  
  # lr1 - logret for 1 portfolio (stock + futures)
  # os1 - out-of-sample length for 1 portfolio (stock + futures)
  
  trueVol <- array(data=NA, dim = c(2,2,os1))
  for (j in 1:os1) {
    ly <- nrow(lr1)-os1
    lr2 <- as.matrix(lr1[,2:3])
    trueVol[,,j] <- t(lr2[(ly+j-1):(ly+j),])%*%lr2[(ly+j-1):(ly+j),]
  }
  trueVol
}

CalcLoss2 <- function(trueV1, forcV1) {
  # loss2 - penalty for undervalued volatility
  # trueV1 - true volatility for 1 portfolio
  # forcV1 - forecasted volatility for 1 portfolio
  
  k <- dim(trueV1)[3]
  loss2 <- rep(NA, k)
  for (i1 in 1:k) {
    s. <- (forcV1[,,i1])^(-1) %*% trueV1[,,i1]
    loss2[i1] <- sum(diag(s.))-log(norm(s.))-ncol(forcV1)  
  }
  loss2
}

CalcLossDiffs <- function(losses) {
  # calculate loss function differences
  l <- length(losses)
  d <- vector('list', l) # list of matrices; each matrix contains loss function differences
  for (i1 in 1:l) {
    d1 <- matrix(NA, dim(losses[[i1]])[1], 3)
    dimnames(d1)[[2]] <- c('adcc-cop', 'adcc-gg', 'cop-gg')
    d1[,1] <- losses[[i1]][,1]-losses[[i1]][,2]
    d1[,2] <- losses[[i1]][,1]-losses[[i1]][,3]
    d1[,3] <- losses[[i1]][,2]-losses[[i1]][,3]
    d[[i1]] <- d1
    d1 <- NULL
  }
  d
}

CalcDiebMarTest <- function(d, median) {
  # compare mean values of loss function differences
  # d - list of loss function differences
  # median - logical; if TRUE Diebold-Mariano sign test is calculated
  l <- length(d)
  dm <- matrix(NA, l, 3)
  for (i1 in 1:l) {
    if(!median) {
      dMean <- apply(d[[i1]], 2, mean)
      dVar <- apply(d[[i1]], 2, lrvar, type = "Andrews")
      dm[i1,] <- dMean/sqrt(dVar)
    }
    else {
      dPlusSum <- apply(d[[i1]], 2, function(x) {sum(x>0)})
      os <- nrow(d[[i1]]) # number of observations in the out-of-sample subsample
      dm[i1,] <- (dPlusSum - 0.5*os)/sqrt(0.25*os)
    }
  }
  dimnames(dm) <- dimnames(d[[1]])
  dm
}

CalcProb <- function(dm) {
  dmProb <- pnorm(abs(c(dm)), lower.tail = F)
  dmProb <- matrix(dmProb,dim(dm)[1], dim(dm)[2])
  dimnames(dmProb) <- dimnames(dm)
  dmProb
}

# SRACH!!
printGgCoefsSe <- function() {
  ns <- matrix(NA, 6,9)
  for (i1 in 1:9) {
    if (ggResult[[1]][[i1]]@mfit$ufit@fit[[1]]@fit$robust.matcoef['skew',4]<0.05) {
      ns[1:2,i1] <- c(ggResult[[1]][[i1]]@mfit$ufit@fit[[1]]@fit$robust.matcoef['skew',1],
                      ggResult[[1]][[i1]]@mfit$ufit@fit[[1]]@fit$robust.matcoef['skew',4])
    }
    if (ggResult[[1]][[i1]]@mfit$ufit@fit[[1]]@fit$robust.matcoef['shape',4]<0.05) {
      ns[3:4,i1] <- c(ggResult[[1]][[i1]]@mfit$ufit@fit[[1]]@fit$robust.matcoef['shape',1],
                      ggResult[[1]][[i1]]@mfit$ufit@fit[[1]]@fit$robust.matcoef['shape',4])
    }
    ns[5:6,i1] <- c(ggResult[[1]][[i1]]@mfit$ufit@fit[[1]]@fit$robust.matcoef['ghlambda',1],
                    ggResult[[1]][[i1]]@mfit$ufit@fit[[1]]@fit$robust.matcoef['ghlambda',4])
  }
  rm(i1)
  dimnames(ns) <- list(c('skew', 'skew-prob', 'shape', 'shape-prob', 'ghlambda', 'ghlambda-prob'), names)
  ns
}

mse.fut <- function(fit.obj1, lr1) {
  # fit.obj1 - fit object for 1 asset (matrix)
  # lr1 - logret of 1 asset
  
  f <- fitted(fit.obj1)[,2] # futures' returns
  l <- nrow(f)
  lr1 <- tail(as.matrix(lr1[,3]),l)
  mse <- 1/l*sum((f-lr1)^2)
  mse
}

# coefs of unigarch in gg with se
# ggResult[[1]][[3]]@mfit$ufit@fit[[2]]@fit$robust.matcoef
# mean equation coefs
# ggResult[[1]][[3]]@mfit$arcoef
# variance equation coefs
# ggResult[[1]][[3]]@mfit$garchcoef
# garchcoef are the same as robust.matcoef but without se

# coefs of copgarch (se not estimated because of transformation='empirical'?)
# copResult[[1]][[3]]@mfit$matcoef
# all the params without se
# copResult[[1]][[3]]@model$mpars 

# but adcc is almost the same model!
# adccResult[[1]][[3]]@mfit$matcoef

#### t-distribution with degrees of freedom vector #####
t2.ll <- function(par, r) {
  # log likelihood for t-distribution with degrees of freedom vector 
  # for 2 assets!
  
  # par - vector of parameters
  # r - innovations, Txn, n=2
  
  v <- par[1:2] # vector of degrees of freedom, nx1, n=2
  alp <- par[3:5] # lower triangle of constant volatility matrix
  bet <- par[6:9] # arch effect
  gam <- par[10:13] # garch effect 
  
  # to matrices
  tmp <- matrix(0,2,2)
  tmp[lower.tri(tmp, T)] <- alp
  alp <- tmp %*% t(tmp)
  rm(tmp)
  bet <- matrix(bet, 2,2)
  gam <- matrix(gam, 2,2)
  
  hi <- array(NA, dim=c(2,2,nrow(r)))
  hi[,,1] <- var(r)
  lli <- rep(NA, nrow(r)) # addons of ll function
  for (i1 in 1:nrow(r)) {
    # bekk matrix
    # browser()
    if(i1>1) {
    hi[,,i1] <- alp + bet %*% t(r[i1,]) %*% r[i1,] %*% t(bet) + gam %*% hi[,,i1-1] %*% t(gam) }
    if(sum(eigen(hi[,,i1], only.values = T)$values<0)) {
      print(paste('volatility bad cond', i1))
    }
    
    pt <- sqrt((2*v[2]-3)/2)*
      matrix(c(sqrt(hi[1,1,i1]), hi[1,2,i1]*1/sqrt(hi[1,1,i1]),
               0, sqrt((v[1]-1)/(v[2]-1)*(hi[2,2,i1]-hi[1,2,i1]^2/hi[1,1,i1]))),2,2)
    if(sum(eigen(pt, only.values = T)$values<0)) {
      print(paste('pt bad cond', i1))
    }
    
    ai <- pt %*% t(pt)
    if(sum(eigen(ai, only.values = T)$values<0)) {
      print(paste('ai bad cond', i1))
    }
    
    lli[i1] <- lgamma(v[1]+1/2)+lgamma(v[2])-lgamma(v[1])-lgamma(v[2]-1/2)-
      1/2*log(det(ai))-(v[1]+1/2)*log(1+1/2* r[i1,] %*% solve(ai) %*% t(r[i1,]))+
      (v[1]-v[2])*log(1+1/2*r[i1,1]/ai[1,1])
  }
  -sum(lli)
}





