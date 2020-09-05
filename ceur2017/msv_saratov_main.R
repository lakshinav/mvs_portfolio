setwd('f:/Google Диск/Nauka/stoch_vol/r_stoch_vol/saratov2017')
localpath <- 'd:/stoch_vol_homestation/'

require(xts)
require(rmgarch)
require(beepr) # beep(5)

# load('msv_saratov_start.RData') # OCHL + volume
load('msv_saratov_2.RData')
load('d:/stoch_vol_homestation/msv_saratov_3.RData')
load('f:/1work/stoch_vol_local/msv_saratov_3.RData')

source('msv_saratov_func.r')

# prepare data ####
files <- dir()
tickers <- sapply(files, function(x) {
  regmatches(x, regexpr('(?<=SPFB\\.)[A-Z]{2,4}', x, perl = T))
})
names(tickers) <- NULL
tickers <- unlist(tickers[-1])
tickers <- tolower(tickers)
rm(files)

tickers <- c("alrs", "chmf", "fees", "gazp", "gmkn", "hydr", "lkoh", "mgnt", "nlmk",
             "nvtk", "rosn", "rtkm", "sber", "sngs", "tatn", "trnf", "urka", "vtbr", "yndx")
n <- length(tickers)

dataFromTxt <- readAll(tickers)

sapply(1:length(dataFromTxt), function(x) 
  {plot(dataFromTxt[[x]], main = names(dataFromTxt)[x])}) # plot all

mergedData <- mergeData(dataFromTxt)

prices <- accountLotSize(mergedData)
plot(prices$fees)
plot(prices$sber)

lr <- calcLogRet(prices) # log returns
plot(lr$gazp)
plot(lr$sber)

sapply(lr, length)/2

rm("dataFromTxt", "mergedData", "prices")

# garch ####
est_gg <- estimGoGarch(lr[c(-1,-n)]) # alrs and yndx have too few observations
est_gg[['fit']][[3]]@mfit$ufit@fit[[1]]@fit$robust.matcoef
plot(fitted(est_gg[['fit']][[15]])[,2])

ar_spec <- lapply(est_gg[['fit']], function(x) {rownames(x@model$arcoef)})
# => all models have ARMA(1,0) mean model

est_adcc <- estimAdcc(lr[c(-1,-n)])
est_adcc[['fit']][[3]]@mfit$matcoef
plot(fitted(est_adcc[['fit']][[15]])[,2])

est_cgarch <- estimCopGarch(lr[c(-1,-n)]) 

############ aux ####
# arrange est_msvt
na
files_dir <- grep('res_separ1', dir(path = mydir), value = T)
files_dir <- gsub('res_separ1_', replacement = "", files_dir)
files_dir <- gsub('\\.RData', replacement = "", files_dir)
length(files_dir)

files_env <- grep('est_msvt_', ls(), value = T)
files_env <- gsub('est_msvt_', replacement = "", files_env)
files_env <- files_env[!files_env == 'lkoh_2']
length(files_env)

needed <- sapply(files_dir, function(x) {x %in% files_env})
which(needed)

sapply(files_dir, function(x) {
  if(!(x %in% files_env)) {
    load(paste0(mydir, 'res_separ1_', x, '.RData'), envir = .GlobalEnv)
  }
})
########## end aux

# decsritpive stat ####
require(moments)
ds <- matrix(NA, nrow=n-2, ncol=1+2*4) # N,mean,sd,skew,kurt for stock and fut
dimnames(ds) <- list(names(lr[c(-1,-n)]), 
                     c('N', rep(c('Mean', 'St.dev.', 'Skew.', 'Kurt.'),2)))
for (i1 in 2:(n-1)) {
  # stocks
  ds[i1-1,1] <- nrow(lr[[i1]])
  ds[i1-1,2] <- mean(lr[[i1]][,1])
  ds[i1-1,3] <- sd(lr[[i1]][,1])
  ds[i1-1,4] <- skewness(lr[[i1]][,1])
  ds[i1-1,5] <- kurtosis(lr[[i1]][,1])
  # futures
  ds[i1-1,2+4] <- mean(lr[[i1]][,2])
  ds[i1-1,3+4] <- sd(lr[[i1]][,2])
  ds[i1-1,4+4] <- skewness(lr[[i1]][,2])
  ds[i1-1,5+4] <- kurtosis(lr[[i1]][,2])
}
rm(i1)

# msv ####
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

file.edit('./msv_t.stan')

# d <- as.matrix(lr[[1]]*1e2)
d <- as.matrix(residuals(est_gg[['fit']][[1]]))
TT <- nrow(d)
nn <- ncol(d)

msv_t1<- stan( # no autoregression in volatility (no M matrix)
  file = "./msv_t.stan", 
  data = list(TT = TT, n=nn, y=d, oos=est_gg[['forc']][[1]]@model$n.roll+1, mu=rep(0, 2)),
  chains = 2,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 3000,            # total number of iterations per chain
  # cores = 2,              # number of cores (using 1 just for the vignette)
  refresh = 2000          # show progress every 'refresh' iterations
)
beep(1)
rm(vc,L,d,TT,n)

# get pars mean
pars <- get_posterior_mean(msv_t1)

# diagnostics
# rstan
options(digits=7)
summary(msv_t1, pars = c('H', 'y_gen[100,1]', 'y_gen[100,2]'), probs = c(0.025, 0.975), digits=3)
stan_plot(msv_t1, pars = c('H[1,1]', 'y_gen[100,1]', 'y_gen[100,2]'))
stan_trace(msv_t1, nrow=3, pars = c('H[1,1]', 'H[2,1]', 'y_gen[100,1]', 'y_gen[100,2]'))
stan_dens(msv_t1,separate_chains = T, pars = c('H[1,1]', 'H[2,1]', 'y_gen[100,1]', 'y_gen[100,2]'))
stan_scat(msv_t1, pars = c('H[1,1]', 'H[1,2]'))
stan_diag(msv_t1)
stan_par(msv_t1, par = c('H[1,1]'))
stan_rhat(msv_t1, bins=35)
stan_ess(msv_t1)
stan_mcse(msv_t1)
stan_ac(msv_t1, pars = c('H[1,1]','H[2,1]', 'y_gen[100,1]', 'y_gen[100,2]'))
stan_diag(msv_t1,
          information = c("sample","stepsize", "treedepth","divergence")[3])

# coda diagnostics
require(coda)
coda.options(digits = 3) 
coda.options(default=TRUE)

msv_t1_coda <- As.mcmc.list(msv_t1, pars = c('H', 'y_gen[100,1]', 'y_gen[100,2]')) # convert to coda object 'y_gen[100,1]', 'y_gen[100,2]'
msv_t1_mcmc <- as.mcmc(msv_t1)

plot(msv_t1_coda)
summary(msv_t1_coda)
crosscorr.plot(msv_t1_coda)
betterPairs(as.matrix(msv_t1_coda[[1]])[, c(1,2,5,6)])
autocorr.plot(msv_t1_coda) # just AC plot
autocorr.diag(msv_t1_coda) # just AC
# cumuplot(msv_t1_coda) # ???
effectiveSize(msv_t1_coda) # Sample size adjusted for autocorrelation (needed for se calculation)
lapply(msv_t1_coda,effectiveSize) # for each chain

# errors (take overdispersed inits)
gelman.diag(msv_t1_coda, multivariate=T) # PSRF
gelman.plot(msv_t1_coda)

geweke.diag(msv_t1_coda) # test for equality of the means of the first and last part of a Markov chain
geweke.plot(msv_t1_coda)

heidel.diag(msv_t1_coda)

codamenu()
# 1: Geweke’s Z -scores: convergence: if < 1.96 -> converged
# 2: Gelman and Rubin: PSRF (potential scale reduction factor, see lecture_on_MCMC.pdf), R_hat?: convergence; > 1 is bad
# 3: Raftery and Lewis
# 4: Heidelberger and Welch: convergence (like Geweke, but using Cramer-von-Mises statistic instead of t)
# + DIC
# + credible intervals: http://people.stat.sc.edu/Hitchcock/stat535slidesday3.pdf
# ESS (aktivnyi razmer byborki): if ESS << S (number of samples), it's bad (see Tsyplakov, 2007); n with respect to correlation
# ? inefficiency factors
# IACT (autocorrelation time, see Yu, Meyer), rho_k: high AC = slow mixing
# traceplots: should be as WN for good convergence + mixing (WN is good)
# if acceptance rate is high, it's bad (too small scale in tuning parameters); shoud be 20-70%
# density: not porcupine; if multi-model, then non-identifiability; plot each chain separately to see convergence
# ASDSF, AWTY (no in coda, see rwty package): convergence

# hedge ####
# rra <- c(0.01, 0.05, 0.2, 1)
# rra <- 0.2
rra <- 0.05
hr_gg_fit <- calcHedgeRat(est_gg[['fit']], forecast = F, rra=rra)
hr_gg_forc <- calcHedgeRat(est_gg[['forc']], forecast = T, rra=rra)
# plot(hr_gg_fit[[5]])
# max(hr_gg_forc[[5]])
# sapply(hr_gg_fit, mean)

hr_adcc_fit <- calcHedgeRat(est_adcc[['fit']], forecast = F, rra=rra)
hr_adcc_forc <- calcHedgeRat(est_adcc[['forc']], forecast = T, rra=rra)
plot(hr_adcc_fit[[5]])
max(hr_adcc_forc[[5]])
sapply(hr_adcc_fit, mean)

hr_cgarch_fit <- calcHedgeRat(est_cgarch[['fit']], forecast = F, rra=rra)
hr_cgarch_forc <- calcHedgeRatCgarch(est_cgarch[['forc']], rra=rra)
plot(hr_cgarch_fit[[5]])
max(hr_cgarch_forc[[5]])
sapply(hr_cgarch_fit, mean)

# hedge ratios with infinite risk aversion
hr_gg_fit_infrra <- calcHedgeRat(est_gg[['fit']], forecast = F, rra=Inf)
hr_gg_forc_infrra <- calcHedgeRat(est_gg[['forc']], forecast = T, rra=Inf)
plot(hr_gg_fit_infrra[[5]])
max(hr_gg_forc_infrra[[5]])
sapply(hr_gg_fit_infrra, mean)

hr_adcc_fit_infrra <- calcHedgeRat(est_adcc[['fit']], forecast = F, rra=Inf)
hr_adcc_forc_infrra <- calcHedgeRat(est_adcc[['forc']], forecast = T, rra=Inf)
plot(hr_adcc_fit_infrra[[5]])
max(hr_adcc_forc_infrra[[5]])

hr_cgarch_fit_infrra <- calcHedgeRat(est_cgarch[['fit']], forecast = F, rra=Inf)
hr_cgarch_forc_infrra <- calcHedgeRatCgarch(est_cgarch[['forc']], rra=Inf)
plot(hr_cgarch_fit_infrra[[5]])
max(hr_cgarch_forc_infrra[[5]])

# constant hedge ratios (OLS)
hr_const_fit <- calcHedgeRatConst(lr[c(-1,-n)], est_gg[['fit']], forecast = F, rra = rra)
hr_const_forc <- calcHedgeRatConst(lr[c(-1,-n)], est_gg[['forc']], forecast = T, rra = rra)
unlist(hr_const_fit)

## calculate hedged returns
hret_gg_fit <- CalcHedgedRet(lr[c(-1,-n)], hr_gg_fit)
hret_gg_forc <- CalcHedgedRet(lr[c(-1,-n)], hr_gg_forc)
plot(hret_gg_fit[[10]])
plot(hret_gg_forc[[5]])
hist(hret_gg_forc[[5]])

hret_adcc_fit <- CalcHedgedRet(lr[c(-1,-n)], hr_adcc_fit)
hret_adcc_forc <- CalcHedgedRet(lr[c(-1,-n)], hr_adcc_forc)
plot(hret_adcc_fit[[5]])
max(hret_adcc_forc[[5]])

hret_cgarch_fit <- CalcHedgedRet(lr[c(-1,-n)], hr_cgarch_fit)
hret_cgarch_forc <- CalcHedgedRet(lr[c(-1,-n)], hr_cgarch_forc)
plot(hret_cgarch_fit[[5]])
max(hret_cgarch_forc[[5]])

# hedge returns with infinite risk aversion
hret_gg_fit_infrra <- CalcHedgedRet(lr[c(-1,-n)], hr_gg_fit_infrra)
hret_gg_forc_infrra <- CalcHedgedRet(lr[c(-1,-n)], hr_gg_fit_infrra)
plot(hret_gg_fit_infrra[[5]])
max(hret_gg_forc_infrra[[5]])

hret_adcc_fit_infrra <- CalcHedgedRet(lr[c(-1,-n)], hr_adcc_fit_infrra)
hret_adcc_forc_infrra <- CalcHedgedRet(lr[c(-1,-n)], hr_adcc_forc_infrra)
plot(hret_adcc_fit_infrra[[5]])
max(hret_adcc_forc_infrra[[5]])

hret_cgarch_fit_infrra <- CalcHedgedRet(lr[c(-1,-n)], hr_cgarch_fit_infrra)
hret_cgarch_forc_infrra <- CalcHedgedRet(lr[c(-1,-n)], hr_cgarch_forc_infrra)
plot(hret_cgarch_fit_infrra[[5]])
max(hret_cgarch_forc_infrra[[5]])

# constant hedge returns (OLS)
hret_const_fit <- CalcHedgedRetConst(lr[c(-1,-n)], est_gg[[1]],forecast = F, hr_const_fit)
hret_const_forc <- CalcHedgedRetConst(lr[c(-1,-n)], est_gg[[2]],forecast = T, hr_const_forc)
plot(hret_const_fit[[5]])
max(hret_const_forc[[5]])

## calculate efficiency
eff_gg_fit <- calcEff(lr[c(-1,-n)], hret_gg_fit)
eff_gg_forc <- calcEff(lr[c(-1,-n)], hret_gg_forc)
plot(eff_gg_fit[[5]])
max(eff_gg_forc[[5]])
unlist(eff_gg_forc)

eff_adcc_fit <- calcEff(lr[c(-1,-n)], hret_adcc_fit)
eff_adcc_forc <- calcEff(lr[c(-1,-n)], hret_adcc_forc)
plot(eff_adcc_fit[[5]])
max(eff_adcc_forc[[5]])
unlist(eff_adcc_forc)

eff_cgarch_fit <- calcEff(lr[c(-1,-n)], hret_cgarch_fit)
eff_cgarch_forc <- calcEff(lr[c(-1,-n)], hret_cgarch_forc)
plot(eff_cgarch_fit[[5]])
max(eff_cgarch_forc[[5]])
unlist(eff_cgarch_forc)

# hedge ratios with infinite risk aversion
eff_gg_fit_infrra <- calcEff(lr[c(-1,-n)], hret_gg_fit_infrra)
eff_gg_forc_infrra <- calcEff(lr[c(-1,-n)], hret_gg_fit_infrra)
plot(eff_gg_fit_infrra[[5]])
max(eff_gg_forc_infrra[[5]])

eff_adcc_fit_infrra <- calcEff(lr[c(-1,-n)], hret_adcc_fit_infrra)
eff_adcc_forc_infrra <- calcEff(lr[c(-1,-n)], hret_adcc_forc_infrra)
plot(eff_adcc_fit_infrra[[5]])
max(eff_adcc_forc_infrra[[5]])

eff_cgarch_fit_infrra <- calcEff(lr[c(-1,-n)], hret_cgarch_fit_infrra)
eff_cgarch_forc_infrra <- calcEff(lr[c(-1,-n)], hret_cgarch_forc_infrra)
plot(eff_cgarch_fit_infrra[[5]])
max(eff_cgarch_forc_infrra[[5]])

# constant hedge ratios (OLS)
eff_const_fit <- calcEff(lr[c(-1,-n)], hret_const_fit)
eff_const_forc <- calcEff(lr[c(-1,-n)], hret_const_forc)
max(unlist(eff_const_fit))
unlist(eff_const_forc)

# hedge ratios for msv
# here forecast for volatility only
hr_msv_forc <- calcHedgeRatMsv(na, rra = rra)
plot(hr_msv_forc[[1]], type='l')
mean(hr_msv_forc[[1]])

hret_msv_forc <- CalcHedgedRet(lr[c(-1,-n)], hr_msv_forc)
plot(hret_msv_forc[[15]], type='l')

eff_msv_forc <- calcEff(lr[c(-1,-n)], hret_msv_forc)

# short version ####
# full version see above
# rra <- c(0.05, 0.2, 4)
rra <- 4

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

# comparison ####
eff_all <- matrix(NA, nrow=n-2, ncol=4)
dimnames(eff_all) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
eff_all[,1] <- unlist(eff_adcc_forc)
eff_all[,2] <- unlist(eff_gg_forc)
eff_all[,3] <- unlist(eff_cgarch_forc)
eff_all[,4] <- unlist(eff_msv_forc)
# eff_all[,5] <- unlist(eff_const_forc)
sum(apply(eff_all, 1, function(x) {x[4] == max(x)}))

profit <- matrix(NA, nrow=n-2, ncol=4)
dimnames(profit) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
profit[,1] <- sapply(hret_adcc_forc, sum)
profit[,2] <- sapply(hret_gg_forc, sum)
profit[,3] <- sapply(hret_cgarch_forc, sum)
profit[,4] <- sapply(hret_msv_forc, sum)
sum(apply(profit, 1, function(x) {x[4] == max(x)}))

# utility
util <- matrix(NA, nrow=n-2, ncol=4)
dimnames(util) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
util[,1] <- sapply(hret_adcc_forc, function(x) {mean(x) - var(x)*rra/2})
util[,2] <- sapply(hret_gg_forc, function(x) {mean(x) - var(x)*rra/2})
util[,3] <- sapply(hret_cgarch_forc, function(x) {mean(x) - var(x)*rra/2})
util[,4] <- sapply(hret_msv_forc, function(x) {mean(x) - var(x)*rra/2})
sum(apply(util, 1, function(x) {x[4] == max(x)}))

# smth like VaR
tipa_valatrisk <- matrix(NA, nrow=n-2, ncol=4)
dimnames(tipa_valatrisk) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
tipa_valatrisk[,1] <- sapply(hret_adcc_forc, function(x) {-quantile(x, prob=0.05)})
tipa_valatrisk[,2] <- sapply(hret_gg_forc, function(x) {-quantile(x, prob=0.05)})
tipa_valatrisk[,3] <- sapply(hret_cgarch_forc, function(x) {-quantile(x, prob=0.05)})
tipa_valatrisk[,4] <- sapply(hret_msv_forc, function(x) {-quantile(x, prob=0.05)})

# summary
rra <- c(0.05, 0.2, 4)
rra <- seq(1e-4, 10, length.out = 10)

summ_rra <- vector('list', 4)
for (i2 in 1:length(summ_rra)) {
  summ_rra_tmp <- matrix(NA, nrow=length(rra), ncol=3)
  dimnames(summ_rra_tmp) <- list(rra, c('eff', 'profit', 'util')) 
  for (i1 in 1:length(rra)) {
    summ_rra_tmp[i1,] <- CountLeader(rra[i1], i2) 
    # 1 - ADCC, 2 - GO-GARCH, 3 - cop-GARCH, 4 - MSV
  }
  summ_rra[[i2]] <- summ_rra_tmp
}
summ_rra
rm(i1, i2, summ_rra_tmp)

# likelihood
ll <- matrix(NA, nrow=n-2, ncol=4)
dimnames(ll) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
for (i1 in 1:length(na)) {
  tmp <- get(paste0('est_msvt_', na[i1]), envir = .GlobalEnv)
  tmp <- get_posterior_mean(tmp)['lp__',3]
  ll[i1,] <- c(likelihood(est_adcc[[1]][[i1]]),
               likelihood(est_gg[[1]][[i1]]),
               likelihood(est_cgarch[[1]][[i1]]),
               tmp)
}
rm(i1, tmp)
sum(apply(ll,1,function(x) {x[4] == max(x)}))

bic <- matrix(NA, nrow=n-2, ncol=4)
dimnames(bic) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
for (i1 in 1:length(na)) {
  bic[i1,] <- c('todo')
}
rm(i1, tmp)
sum(apply(ll,1,function(x) {x[4] == max(x)}))

# mse
mse_adcc <- calcMse(est_adcc$forc, lr, 'adcc')
mse_cgarch <- calcMse(est_cgarch$forc, lr, 'cgarch')
mse_gg <- calcMse(est_gg$forc, lr, 'gg')
mse_msv <- calcMseMsv(na, lr)

plot(mse_adcc, type='b', col='red', ylim=c(4,25))
par(new=T)
plot(mse_cgarch, type='b', col='yellow', ylim=c(4,25))
par(new=T)
plot(mse_gg, type='b', col='green', ylim=c(4,25))
par(new=F)

# mad
mad_adcc <- calcMad(est_adcc$forc, lr, 'adcc')
mad_cgarch <- calcMad(est_cgarch$forc, lr, 'cgarch')
mad_gg <- calcMad(est_gg$forc, lr, 'gg')
mad_msv <- calcMadMsv(na, lr)


# mean hr
hr_mean <-  matrix(NA, nrow=n-2, ncol=4)
dimnames(hr_mean) <- list(na, c('ADCC', 'GO-GARCH', 'cop-GARCH', 'MSV'))
for (i1 in 1:length(na)) {
  hr_mean[i1,] <- c(mean(hr_adcc_forc[[i1]]),
                    mean(hr_gg_forc[[i1]]),
                    mean(hr_cgarch_forc[[i1]]),
                    mean(hr_msv_forc[[i1]]))
}
hr_mean
rm(i1)

# Geweke test for stanfit
gew_test <- matrix(NA, nrow=n-2, ncol=3)
dimnames(gew_test) <- list(na, c('Sigma[2,1]', 'Sigma[2,2]', 'y_fut'))
for (i1 in 1:length(na)) {
  msv_t1_coda <- As.mcmc.list(get(paste0('est_msvt_', na[i1])), 
                              pars = c('Sigma[10,2,1]', 'Sigma[10,2,2]', 'y_gen[100,2]')) 
  tmp <- geweke.diag(msv_t1_coda)
  gew_test[i1,] <- pnorm(abs(tmp[[1]]$z))
}
gew_test
rm(i1,tmp,msv_t1_coda)



# plot
barplot(apply(t(eff_all), 2, function(x) {x-min(x)}), 
        beside = T, col = heat.colors(5)[5:1])
barplot(apply(t(profit), 2, function(x) {x-min(x)}), 
        beside = T, col = heat.colors(4)[4:1])
barplot(apply(t(tipa_valatrisk), 2, function(x) {x-max(x)}), 
        beside = T, col = heat.colors(4)[4:1])

# output ####
# excel
prnt <- tickers
write.table(sprintf('%s', prnt), #%0.3f
            file='clipboard',
            quote = F, sep = '\t', row.names = F, col.names = F) #eol = '\t'
rm(prnt)

write(names(lr), file='clipboard') # easy clipboard

m1 <- paste('&', sprintf('%0.3f', ds))
m1 <- matrix(m1, nrow = nrow(ds), ncol = ncol(ds), byrow = F)
m2 <- paste0('\\textsc{', names(lr[c(-1,-n)]), '}')
write.table(cbind(m2, m1), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' \\\\ \n')
rm(m1,m2)

write.table(paste('&', colnames(ds)), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' ')


m1 <- paste('&', sprintf('%0.3f', eff_all))
m1 <- matrix(m1, nrow = nrow(eff_all), ncol = ncol(eff_all), byrow = F)
# m2 <- paste0('\\textsc{', names(lr[c(-1,-n)]), '}')
write.table(cbind(toupper(na), m1), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' \\\\ \n')
rm(m1,m2)
write.table(paste('&', colnames(eff_all)), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' ')

# profit
m1 <- paste('&', sprintf('%0.3f', profit))
m1 <- matrix(m1, nrow = nrow(profit), ncol = ncol(profit), byrow = F)
# m2 <- paste0('\\textsc{', names(lr[c(-1,-n)]), '}')
write.table(cbind(toupper(na), m1), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' \\\\ \n')
rm(m1,m2)
write.table(paste('&', colnames(profit)), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' ')
# summ_rra
m1 <- paste('&', sprintf('%0.0f', unlist(summ_rra)))
m1 <- matrix(m1, nrow=10, ncol=12, byrow = F)
write.table(cbind(sprintf('%0.2f',rra), m1), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' \\\\ \n')
write.table(paste('&', rep(colnames(summ_rra[[1]]),4)), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' ')
# ll
m1 <- paste('&', sprintf('%0.3f', ll))
m1 <- matrix(m1, nrow=n-2, ncol=4, byrow = F)
m2 <- paste0('\\textsc{', names(lr[c(-1,-n)]), '}')
write.table(cbind(m2, m1), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' \\\\ \n')
write.table(paste('&', colnames(profit)), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' ')
# mean hr
m1 <- paste('&', sprintf('%0.3f', hr_mean))
m1 <- matrix(m1, nrow=n-2, ncol=4, byrow = F)
m2 <- paste0('\\textsc{', names(lr[c(-1,-n)]), '}')
write.table(cbind(m2, m1), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' \\\\ \n')
write.table(paste('&', colnames(profit)), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' ')
rm(m1,m2)
# gew_test
m1 <- paste('&', sprintf('%0.3f', gew_test))
m1 <- matrix(m1, nrow=n-2, ncol=3, byrow = F)
m2 <- paste0('\\textsc{', names(lr[c(-1,-n)]), '}')
write.table(cbind(m2, m1), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' \\\\ \n')
write.table(paste('&', colnames(profit)), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' ')
rm(m1,m2)






# in the end of the day ####
save.image('d:/stoch_vol_homestation/msv_saratov_3.RData')
save.image('f:/1work/stoch_vol_local/msv_saratov_4_mvnorm.RData')

