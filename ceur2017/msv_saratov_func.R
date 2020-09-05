readTxt <- function(ticker) {
  # read data from *.txt got from Finam.ru
  fileName <- paste0(ticker, "_000101_161231.txt")
  readedData <- read.table(fileName, check.names = FALSE, header=F,
                           skip=1, sep=",", dec=".", 
                           colClasses = c("character", rep("NULL", 4), "double", "NULL"))
  names(readedData) <- c("date", "close")
  # readedData$date <- )
  res <- xts(x = readedData$close, order.by = as.Date(readedData$date, format="%Y%m%d"))
  # res <- xts(x = readedData$close, order.by = as.POSIXct(readedData$date, format="%Y%m%d"))
  
  return(res)
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
      ticker <- paste0("SPFB.", toupper(names[i1-l])) # ticker for futures
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
    # mergedData <- merge(data[[i1]], data[[i1+l]], by=1, 
    #                     suffixes = c(".stock",".fut"))
    mergedData <- merge(data[[i1]], data[[i1+l]], join='inner',
                        suffixes = c('stock', 'fut'))
    readyData[[i1]] <- mergedData
  }
  readyData
}

accountLotSize <- function(data) {
  lotSize <- c(100, 100, 100000,  100, 10, 10000, 10, 1, 100, 100, 100, 100, 100, 1000, 100, 1, 100, 100000, 100) # lot size for futures
  for (i1 in 1:length(data)) {
    if(names(data)[i1]=="sber") {
      data[[i1]]['/2007-07-17',] <- data[[i1]]['/2007-07-17',]/1e3 # account for split
      # acoount for lot size
      data[[i1]]['2007-07-18/2007-12-10',2] <- data[[i1]]['2007-07-18/2007-12-10',2]/1e3
      data[[i1]]['2007-12-11/',2] <- data[[i1]]['2007-12-11/',2]/1e2
    } else {
      data[[i1]]$fut <- data[[i1]]$fut/lotSize[i1] 
    }
  }
  data
}

calcLogRet <- function(prices) {
  l <- length(prices)
  lr <- vector("list", l)
  names(lr) <- names(prices)
  for (i1 in 1:l) {
    lr[[i1]] <- diff(log(prices[[i1]]))[-1,]*1e2
  }
  lr
}

estimGoGarch <- function(logret) {
  fit <- list()
  forc <- list()
  l <- length(logret)
  ggSpec <- gogarchspec(mean.model = list(model = "AR", lag.max=3, lag.criterion="SC"),
                        variance.model = list(model = "sGARCH", 
                                              garchOrder = c(1,1)),
                        distribution.model = "mvnorm")
  for (i1 in 1:l) {
    os <- floor(nrow(logret[[i1]])/5)
    ggFit <- gogarchfit(ggSpec,  logret[[i1]], out.sample = os)
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
                             variance.model=list(model="sGARCH",
                                                 garchOrder=c(1,1)))
  multiGarchSpec <- multispec(replicate(2, uniGarchSpec))
  dccSpec <- dccspec(multiGarchSpec,dccOrder=c(1,1), model="aDCC", 
                     distribution="mvnorm")
  for (i1 in 1:l) {
    os <- floor(nrow(logret[[i1]])/5)
    adccFit <- dccfit(dccSpec, logret[[i1]], out.sample = os)
    adccForc <- dccforecast(adccFit, n.ahead=1,n.roll=os)
    fit <- c(fit, adccFit)
    forc <- c(forc, adccForc)
  }
  names(fit) <- names(logret)[1:length(fit)]
  names(forc) <- names(logret)[1:length(forc)]
  res <- list(fit=fit, forc=forc)
}

estimCopGarch <- function(logret) {
  l <- length(logret)
  fit <- list()
  forc <- vector('list', l)

  uniGarchSpec <- ugarchspec(mean.model=list(armaOrder=c(1,0)),
                             variance.model=list(model="fGARCH",submodel="GARCH",
                                                 garchOrder=c(1,1)))
  multiGarchSpec <- multispec(replicate(2, uniGarchSpec))
  copSpec <- cgarchspec(multiGarchSpec,
                        dccOrder = c(1, 1), asymmetric = T,
                        distribution.model = list(copula = "mvnorm",
                                                  method = "Kendall", 
                                                  time.varying = T,
                                                  transformation = "parametric")) # parametric
  for (i1 in 1:l) {
    pair <- logret[[i1]]
    os <- floor(nrow(pair)/5)
    copFit <- cgarchfit(copSpec, pair, out.sample = os, fit.control = list(eval.se=T))
    copForc <- copGarch1AheadForecast(copFit, os, copSpec)
    
    fit <- c(fit, copFit)
    forc[[i1]] <- copForc
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
    # simS[i,] = sqrt(diag(simC[,,i])) # conditional sigma
    # simCor[,,i] = sim3@msim$simR[[1]][,,1] # conditional correlation
    #     print(i) 
    simRes <- list(simMu = simMu, simC = simC) # simS = simS, , simCor = simCor
  }
  simRes
}

calcHedgeRat <- function(est, forecast, rra) {
  # calculate hedge ratio
  
  # est - list with estimation results
  # forecast - logical; if TRUE, calculate with forecasted volatility
  # rra - risk aversion coefficient
  
  m <- length(est)
  hr <- vector('list', m)
  names(hr) <- names(est)
  
  if(forecast) {
    for (i1 in 1:m) {
      tmp_vcov <- rcov(est[[i1]])
      tmp_fitt <- drop(fitted(est[[i1]]))['fut',]
      # index(tmp_fitt) <- as.POSIXct(format(index(tmp_fitt)))
      hr[[i1]] <- (sapply(tmp_vcov, 
                          function(x) {drop(x[1,2,1])}) - tmp_fitt/2/rra)/sapply(tmp_vcov, function(x) {drop(x[2,2,1])})
    }
  } else {
    for (i1 in 1:m) {
      tmp_vcov <- rcov(est[[i1]])
      tmp_fitt <- fitted(est[[i1]])[,'fut']
      index(tmp_fitt) <- as.POSIXct(format(index(tmp_fitt)))
      hr[[i1]] <- (as.xts(tmp_vcov[1,2,]) - tmp_fitt/2/rra)/as.xts(tmp_vcov[2,2,])
    }
  }
  hr
}

calcHedgeRatCgarch <- function(est, rra) {
  # calculate hedge ratio for copula-GARCH forecast
  
  # est - list with estimation results (here forecast only)
  # rra - risk aversion coefficient
  
  m <- length(est)
  hr <- vector('list', m)
  names(hr) <- names(est)
  
  for (i1 in 1:m) {
    tmp_vcov <- est[[i1]]$simC
    tmp_fitt <- est[[i1]]$simMu[,2]
    # index(tmp_fitt) <- as.POSIXct(format(index(tmp_fitt)))
    hr[[i1]] <- (tmp_vcov[1,2,] - tmp_fitt/2/rra)/tmp_vcov[2,2,]
  }
  hr
}

calcHedgeRatConst <- function(lr, est, forecast, rra) {
  # calculate hedge ratio for constant HR

  # lr - log returns
  # est - list with estimation results (here GO-GARCH)
  # forecast - logical; if TRUE, calculate with forecasted volatility
  # rra - risk aversion coefficient
  
  m <- length(est)
  hr <- vector('list', m)
  names(hr) <- names(est)
  
  if(forecast) {
    for (i1 in 1:m) {
      t_os <- length(drop(fitted(est[[i1]]))['fut',])
      t_all <- length(lr[[i1]][,1])
      # index(tmp_fitt) <- as.POSIXct(format(index(tmp_fitt)))
      hr[[i1]] <- (cov(lr[[i1]][(t_all-t_os+1):t_all,])[1,2] - mean(lr[[i1]][(t_all-t_os+1):t_all,2])/2/rra)/var(lr[[i1]][(t_all-t_os+1):t_all,2])
      }
  } else {
    for (i1 in 1:m) {
      t_is <- length(drop(fitted(est[[i1]]))[,'fut'])
      # index(tmp_fitt) <- as.POSIXct(format(index(tmp_fitt)))
      hr[[i1]] <- (cov(lr[[i1]][1:t_is,])[1,2] - mean(lr[[i1]][1:t_is,2])/2/rra)/var(lr[[i1]][1:t_is,2])
     }
  }
  hr
}

.calcHedgeRatMsv1 <- function(msv1, rra) {
  # calculate hedge ratio for MSV-t for 1 asset
  
  # msv - MSV estimation results for 1 asset (stanfit object)
  # rra - risk aversion coefficient
  
  tmp <- get_posterior_mean(msv1)[,3]
  Sigma21 <- tmp[grep('Sigma\\[\\d{1,3},2,1\\]', names(tmp))] # covariance
  Sigma22 <- tmp[grep('Sigma\\[\\d{1,3},2,2\\]', names(tmp))] # futures variance
  futMean <- tmp[grep('y_gen\\[\\d{1,4},2\\]', names(tmp))] # futures mean
  futMean <- tail(futMean, length(Sigma21))
  
  hr1 <- (Sigma21 - futMean/2/rra)/Sigma22
  hr1
}

calcHedgeRatMsv <- function(na, rra) {
  # calculate hedge ratio for MSV-t 
  
  # est - MSV estimation results
  # na - names of stocks
  # rra - risk aversion coefficient
  
  m <- length(na)
  hr <- vector('list', m)
  names(hr) <- na
  
  for (i1 in 1:m) {
    # load(paste0('f:/1work/stoch_vol_local/', 'res_separ1_', na[i1], '.RData'), 
    #      envir = .GlobalEnv)
    tmp <- get(x = paste0('est_msvt_', na[i1]))
    hr[[i1]] <- .calcHedgeRatMsv1(tmp, rra = rra)
    rm(tmp)
  }
  hr
}

CalcHedgedRet <- function(lr, hr) {
  # calculate returns of hedged position
  
  # lr - log returns
  # hr - hedge ratios
  
  l <- length(hr)
  hret <- vector('list', l)
  for (i1 in 1:l) {
    t_hr <- length(hr[[i1]])
    t_lr <- nrow(lr[[i1]])
    if(t_hr > (1/2*t_lr)) { # in-sample
      hret[[i1]] <- lr[[i1]][1:t_hr,1] - hr[[i1]]*lr[[i1]][1:t_hr,2]
    } else { # out-of-sample
      hret[[i1]] <- lr[[i1]][(t_lr-t_hr+1):t_lr,1] - hr[[i1]]*lr[[i1]][(t_lr-t_hr+1):t_lr,2]
    }
  }
  names(hret) <- names(hr)
  hret
}

CalcHedgedRetConst <- function(lr, est, forecast, hr) {
  # calculate returns of hedged position for constant HR
  
  # lr - log returns
  # est - list with estimation results (here GO-GARCH)
  # forecast - logical; if TRUE, calculate with forecasted volatility
  # hr - constant hedge ratios

  l <- length(lr)
  hret <- vector('list', l)
  
  if(forecast) {
    for (i1 in 1:l) {
      t_os <- length(drop(fitted(est[[i1]]))['fut',])
      t_all <- length(lr[[i1]][,1])
      # index(tmp_fitt) <- as.POSIXct(format(index(tmp_fitt)))
      hret[[i1]] <- lr[[i1]][(t_all-t_os+1):t_all,1] - hr[[i1]][1,1]*lr[[i1]][(t_all-t_os+1):t_all,2]
    }
  } else {
    for (i1 in 1:l) {
      t_is <- length(drop(fitted(est[[i1]]))[,'fut'])
      # index(tmp_fitt) <- as.POSIXct(format(index(tmp_fitt)))
      hret[[i1]] <- lr[[i1]][1:t_is,1] - hr[[i1]][1,1]*lr[[i1]][1:t_is,2]
      }
  }
  names(hret) <- names(hr)
  hret
}

# CalcHedgedRetMsv <- function(lr, hr) {
#   # calculate returns of hedged position for MSV model
#   
#   # lr - log returns
#   # hr - hedge ratios
# }


calcEff <- function(lr, hret) { 
  # calculate effectiveness of hedging (hedged vs. unhedged)
  # lr - log returns
  # hret - hedged returns
  
  l <- length(hret)
  eff <- vector('list', l)
  for (i1 in 1:l) {
    t_hr <- nrow(hret[[i1]])
    t_lr <- nrow(lr[[i1]])
    if(t_hr > 1/2*t_lr) { # in-sample
      eff[[i1]] <- 1-var(hret[[i1]])/var(lr[[i1]][1:t_hr,1])
    } else { # out-of-sample
      eff[[i1]] <- 1-var(hret[[i1]])/var(lr[[i1]][(t_lr-t_hr+1):t_lr,1])
    }
  }
  names(eff) <- names(hret)
  eff
}

calcMse <- function(forecast, ret, mod_flag) {
  # calculate MSE
    
  # forecast - list with forecast estimation results
  # ret - observed returns
  # mod_flag - string, name of model; values - c('adcc', 'cgarch', 'gg', 'msv')

  m <- length(forecast)
  mse <- rep(NA, m)
  names(mse) <- names(forecast)
  
  for (i1 in 1:m) {
    switch(mod_flag,
           adcc = {
             # tmp_mse <- rep(NA,nrow(forecast[[i1]]@model$residuals))
             # for (i2 in 1:length(tmp_mse)) {
             #   tmp_mse[i2] <- sum(forecast[[i1]]@model$residuals[i2,]^2) # adcc
             # }
             tmp_mse <- rep(NA,dim(forecast[[i1]]@mforecast$mu)[3])
             tmp_ret <- tail(ret[[i1]], length(tmp_mse))
             for (i2 in 1:length(tmp_mse)) {
               tmp_mse[i2] <- sum((forecast[[i1]]@mforecast$mu[,,i2] - tmp_ret[i2,])^2) # adcc
             }
           },
           cgarch = {
             tmp_mse <- rep(NA,nrow(forecast[[i1]]$simMu))
             tmp_ret <- tail(ret[[i1]], length(tmp_mse))
             for (i2 in 1:length(tmp_mse)) {
               tmp_mse[i2] <- sum((forecast[[i1]]$simMu[i2,] - tmp_ret[i2,])^2) # cgarch
             }
           },
           gg = {
               tmp_mse <- rep(NA,dim(forecast[[i1]]@mforecast$mu)[3])
               tmp_ret <- tail(ret[[i1]], length(tmp_mse))
               for (i2 in 1:length(tmp_mse)) {
                 tmp_mse[i2] <- sum((forecast[[i1]]@mforecast$mu[,,i2] - tmp_ret[i2,])^2) # gg
               }
             }
           )
    
    mse[i1] <- mean(tmp_mse) 
    rm(tmp_mse)
  }
  mse
}

calcMad <- function(forecast, ret, mod_flag) {
  # calculate mean absolute deviation
  
  # forecast - list with forecast estimation results
  # ret - observed returns
  # mod_flag - string, name of model; values - c('adcc', 'cgarch', 'gg', 'msv')
  
  m <- length(forecast)
  mse <- rep(NA, m)
  names(mse) <- names(forecast)
  
  for (i1 in 1:m) {
    switch(mod_flag,
           adcc = {
             # tmp_mse <- rep(NA,nrow(forecast[[i1]]@model$residuals))
             # for (i2 in 1:length(tmp_mse)) {
             #   tmp_mse[i2] <- sum(forecast[[i1]]@model$residuals[i2,]^2) # adcc
             # }
             tmp_mse <- rep(NA,dim(forecast[[i1]]@mforecast$mu)[3])
             tmp_ret <- tail(ret[[i1]], length(tmp_mse))
             for (i2 in 1:length(tmp_mse)) {
               tmp_mse[i2] <- sum(abs(forecast[[i1]]@mforecast$mu[,,i2] - tmp_ret[i2,])) # adcc
             }
           },
           cgarch = {
             tmp_mse <- rep(NA,nrow(forecast[[i1]]$simMu))
             tmp_ret <- tail(ret[[i1]], length(tmp_mse))
             for (i2 in 1:length(tmp_mse)) {
               tmp_mse[i2] <- sum(abs(forecast[[i1]]$simMu[i2,] - tmp_ret[i2,])) # cgarch
             }
           },
           gg = {
             tmp_mse <- rep(NA,dim(forecast[[i1]]@mforecast$mu)[3])
             tmp_ret <- tail(ret[[i1]], length(tmp_mse))
             for (i2 in 1:length(tmp_mse)) {
               tmp_mse[i2] <- sum(abs(forecast[[i1]]@mforecast$mu[,,i2] - tmp_ret[i2,])) # gg
             }
           }
    )
    
    mse[i1] <- mean(tmp_mse) 
    rm(tmp_mse)
  }
  mse
}

calcMseMsv <- function(nam, ret) {
  # calculate MSE for MSV model
  
  # nam - ticker names
  
  if(length(nam) != length(ret)) {
    warning('Lengths of names and logret are not the same')
  }
  
  mse <- rep(NA, length(nam))
  names(mse) <- nam
  
  e1 <- new.env()
  for (i1 in 1:length(nam)) {
    fn <- paste0('f:/1work/stoch_vol_local/', 'res_separ1_', nam[i1], '.RData')
    load(fn, envir = e1)
    tmp_msv <- get(ls(e1), envir = e1)
    tmp_y <- get_posterior_mean(tmp_msv)
    tmp_y <- tmp_y[grep('y_gen', rownames(tmp_y), value = T),3]
    tmp_y <- matrix(tmp_y, ncol=2, byrow = T)
    tmp_ret <- tail(ret[[i1]], nrow(tmp_y))
    tmp_mse <- rep(NA,nrow(tmp_y))
    for (i2 in 1:nrow(tmp_ret)) {
      tmp_mse[i2] <- sum((tmp_ret[i2,] - tmp_y[i2,])^2) 
    }
    mse[i1] <- mean(tmp_mse)
  }
  mse
}

calcMadMsv <- function(nam, ret) {
  # calculate MAD for MSV model
  
  # nam - ticker names
  
  if(length(na) != length(ret)) {
    warning('Lengths of names and logret are not the same')
  }
  
  mse <- rep(NA, length(na))
  names(mse) <- na
  
  e1 <- new.env()
  for (i1 in 1:length(nam)) {
    fn <- paste0('f:/1work/stoch_vol_local/', 'res_separ1_', na[i1], '.RData')
    load(fn, envir = e1)
    tmp_msv <- get(ls(e1), envir = e1)
    tmp_y <- get_posterior_mean(tmp_msv)
    tmp_y <- tmp_y[grep('y_gen', rownames(tmp_y), value = T),3]
    tmp_y <- matrix(tmp_y, ncol=2, byrow = T)
    tmp_ret <- tail(ret[[i1]], nrow(tmp_y))
    tmp_mse <- rep(NA,nrow(tmp_y))
    for (i2 in 1:nrow(tmp_ret)) {
      tmp_mse[i2] <- sum(abs(tmp_ret[i2,] - tmp_y[i2,])) 
    }
    mse[i1] <- mean(tmp_mse)
  }
  mse
}

CountLeader <- function(rra, model) {
  # make summary table for one model by 3 hedge efficiency criteria
  
  # rra - vector with risk aversion coefs
  # model - # of model 
  hr_gg_forc <- calcHedgeRat(est_gg[['forc']], forecast = T, rra=rra)
  hr_adcc_forc <- calcHedgeRat(est_adcc[['forc']], forecast = T, rra=rra)
  hr_cgarch_forc <- calcHedgeRatCgarch(est_cgarch[['forc']], rra=rra)
  # hr_const_forc <- calcHedgeRatConst(lr[c(-1,-n)], est_gg[['forc']], forecast = T, rra = rra)
  
  hret_gg_forc <- CalcHedgedRet(lr[c(-1,-n)], hr_gg_forc)
  hret_adcc_forc <- CalcHedgedRet(lr[c(-1,-n)], hr_adcc_forc)
  hret_cgarch_forc <- CalcHedgedRet(lr[c(-1,-n)], hr_cgarch_forc)
  # hret_const_forc <- CalcHedgedRetConst(lr[c(-1,-n)], est_gg[[2]],forecast = T, hr_const_forc)
  
  eff_gg_forc <- calcEff(lr[c(-1,-n)], hret_gg_forc)
  eff_adcc_forc <- calcEff(lr[c(-1,-n)], hret_adcc_forc)
  eff_cgarch_forc <- calcEff(lr[c(-1,-n)], hret_cgarch_forc)
  # eff_const_forc <- calcEff(lr[c(-1,-n)], hret_const_forc)
  
  hr_msv_forc <- calcHedgeRatMsv(na, rra = rra)
  hret_msv_forc <- CalcHedgedRet(lr[c(-1,-n)], hr_msv_forc)
  eff_msv_forc <- calcEff(lr[c(-1,-n)], hret_msv_forc)
  
  eff_all <- matrix(NA, nrow=n-2, ncol=4)
  dimnames(eff_all) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
  eff_all[,1] <- unlist(eff_adcc_forc)
  eff_all[,2] <- unlist(eff_gg_forc)
  eff_all[,3] <- unlist(eff_cgarch_forc)
  eff_all[,4] <- unlist(eff_msv_forc)
  # eff_all[,5] <- unlist(eff_const_forc)

  profit <- matrix(NA, nrow=n-2, ncol=4)
  dimnames(profit) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
  profit[,1] <- sapply(hret_adcc_forc, sum)
  profit[,2] <- sapply(hret_gg_forc, sum)
  profit[,3] <- sapply(hret_cgarch_forc, sum)
  profit[,4] <- sapply(hret_msv_forc, sum)

  # utility
  util <- matrix(NA, nrow=n-2, ncol=4)
  dimnames(util) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
  util[,1] <- sapply(hret_adcc_forc, function(x) {mean(x) - var(x)*rra/2})
  util[,2] <- sapply(hret_gg_forc, function(x) {mean(x) - var(x)*rra/2})
  util[,3] <- sapply(hret_cgarch_forc, function(x) {mean(x) - var(x)*rra/2})
  util[,4] <- sapply(hret_msv_forc, function(x) {mean(x) - var(x)*rra/2})

  res <- c(sum(apply(eff_all, 1, function(x) {x[model] == max(x)})),
           sum(apply(profit, 1, function(x) {x[model] == max(x)})),
           sum(apply(util, 1, function(x) {x[model] == max(x)})))
  res
}





