readTxt <- function(ticker) {
  # read data from *.txt got from Finam.ru (FEFO)
  fileName <- paste0(getwd(), '/', ticker, "_790101_141231.txt")
  readedData <- read.table(fileName, check.names = FALSE, header=F,
                           skip=1, sep=",", dec=".", 
                           colClasses = c("character", rep("NULL", 4), "double", "NULL"))
  names(readedData) <- c("date", "close")
  readedData$date <- as.Date(readedData$date, format="%Y%m%d")
  return(readedData)
}

readAll <- function(names) {
  # read data from n files in cycle
  l <- length(names)
  readyData <- vector(mode = "list", length = l)  # vector for the result
  tickerAll <- NULL # vector of ticker names
  for (i1 in 1:l) {
    ticker <- toupper(names[[i1]])
    readedData <- readTxt(ticker)    
    readyData[[i1]] <- readedData
  }
  names(readyData) <- names
  readyData # explicit return() is a little bit slower
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

cara2 <- function(pm1, pm2, l) {
  # l - absolute risk aversion
  -exp(-l*pm1) * (1 + (l^2)/2*pm2)
}

cara4 <- function(pm1, pm2, pm3, l) {
  # l - absolute risk aversion
  l2 = (l^2)/2
  l3 = (l^3)/factorial(3)
  -exp(-l*pm1) * (1 + l2*pm2 - l3*pm3)  
}

riskm <- function(r, w, oos) {
  # evaluate risk measures
  
  # r - returns
  # w - weigths without tech-part
  # oos - out-of-sample length
  
  tst <- nrow(w)
  tech <- nrow(r)-tst
  riskcomp <- matrix(NA, tst, 4) # risk measures for ts-tech in-saple and oos out-of-sample points
  dimnames(riskcomp) <- list(c(rep('isa', tst-oos), rep('oos', oos)),
                             c('mad', 'sd', "cvar", "cdar"))
  
  for (i1 in 1:tst) {
    #       browser()
    
    riskcomp[i1,'mad'] <- riskfun(weights = w[i1,], Data = r[1:(tech+i1),], 
                                  risk = 'mad')
    riskcomp[i1,'sd'] <- sqrt(riskfun(weights = w[i1,], Data = r[1:(tech+i1),], 
                                      risk = 'ev'))
    riskcomp[i1,'cvar'] <- riskfun(weights = w[i1,], Data = r[1:(tech+i1),], 
                                   risk = 'cvar', alpha = 0.05)
    #     print(i1)
    riskcomp[i1,'cdar'] <- riskfun(weights = w[i1,], Data = r[1:(tech+i1),], 
                                   risk = 'cdar', alpha = 0.05)
  }  
  riskcomp  
}

plotComp <- function(x, d, ts, tech, oos, ylim, main, ylab, pch) {
  # make comparison plots
  
  # x - what to compare, vector
  # d - dates
  # ts - # of all points
  # tech - points for risk measures calculation
  # oos - out-of-sample
  
  plot(head(d, ts-tech-oos), x[(tech+1):(ts-oos)], 
       pch=pch, col='blue', type='b',
       xlim=c(d[1], tail(d,1)), ylim=ylim,
       xlab='', ylab='')
  par(new=T)
  plot(tail(d, oos), tail(x, oos), 
       pch=pch, col='magenta', type='b',
       xlim=c(d[1], tail(d,1)), ylim=ylim,
       xlab='Time', ylab=ylab, main=main)
  par(new=F)
}

varm <- function(x) {
  # skewness of matrix
  
  # x - matrix with time series (T x n)
  
  n <- ncol(x)
  ts <- nrow(x)
  z <- matrix(NA, nrow=n, ncol=n)
  for (i1 in 1:nrow(x)) {
    for (i2 in 1:ncol(x)) {
      if(i1<=i2) {
        y1 <- x[,i1]
        y2 <- x[,i2]
        z[i1,i2] <- (t(y1-mean(y1)) %*% (y2-mean(y2)))/(ts-1)
      }
    }
  }
  # browser()
  z[lower.tri(z,diag = F)] <- z[upper.tri(z, diag = F)]
  z
}

skewm <- function(x) {
  # skewness of matrix, w
  
  # x - matrix or xts with portfolio returns (T x n)

  n <- ncol(x)
  ts <- nrow(x)
  z <- array(0, c(n,n,n))
  for (i1 in 1:n) {
    for (i2 in 1:n) {
      for (i3 in 1:n)
      if(i1<=i2) {
        y1 <- x[,i1]
        y2 <- x[,i2]
        y3 <- x[,i3]
        z[i1,i2,i3] <- sum((y1-mean(y1))*(y2-mean(y2))*(y3-mean(y3)))/(ts-1)
      }
    }
  }
  for (i4 in 1:n) {
    z[,,i4] <- z[,,i4] + t(z[,,i4]) - diag(diag(z[,,i4]))
  }
  z
}

skewp <- function(s, w, x) {
  # http://quant.stackexchange.com/questions/1557/how-do-i-calculate-the-skewness-of-a-portfolio-of-assets
  # skewness of portfolio
  
  # s - coskewness matrix (tensor n x n x n)
  # w - vector of weights (n x 1)
  # x - returns
  
  m3 <- 0
  for (i1 in 1:n) {
    for (i2 in 1:n) {
      for (i3 in 1:n)
        m3 <- m3 + w[i1]*w[i2]*w[i3]*s[i1,i2,i3]
    }
  }
  # sp <- 1/(sqrt(t(w) %*% var(x) %*% w)^3)*m3
  # sp
  m3
}

skewpt <- function(r, w) {
  # portfolio skewness in time
  
  # r - xts of returns
  # w - matrix of weights (n x ts)

  n <- ncol(r)
  ts <- nrow(r)
  
  s <- rep(0, ts) # vector of portfolio skewnesses
  for (i1 in 1:ts) {
    # browser()
    s[i1] <- skewp(s = skewm(x = r[1:i1, ]), w = w[i1, ], x = r[1:i1, ])
  }
  s
}