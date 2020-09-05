# msv ####
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

file.edit('./msv_t.stan')


# mydir <- 'd:/stoch_vol_homestation/'
mydir <- 'f:/1work/stoch_vol_local/'
# est_msvt <- vector('list', n-2)
# names(est_msvt) <- names(est_gg[[1]])
na <- names(est_gg[[1]])
for (i1 in 6) { 
  d <- as.matrix(residuals(est_gg[['fit']][[i1]]))
  TT <- nrow(d)
  nn <- ncol(d)
  
  msv_t1 <- tryCatch(
    {
      msv_t1 <- stan( # no autoregression in volatility (no M matrix)
        file = "./msv_t.stan",
        data = list(TT = TT, n=nn, y=d, oos=est_gg[['forc']][[i1]]@model$n.roll+1, mu=rep(0, 2)),
        chains = 2,             # number of Markov chains
        warmup = 1000,          # number of warmup iterations per chain
        iter = 2e4,            # total number of iterations per chain
        # cores = 2,              # number of cores (using 1 just for the vignette)
        refresh = -2000          # show progress every 'refresh' iterations
      )
      # msv_t1 <- errTest(i1, d)
      # est_msvt[[i1]] <- msv_t1
      assign(paste0('est_msvt_', na[i1], '_2'), value = msv_t1)
      save(list=paste0('est_msvt_', na[i1]), file = paste0(mydir, 'res_separ2_', na[i1], '.RData'))
      print(i1)
      rm(msv_t1)
      gc()
    },
    warning = function(err) {msv_t1 <- NULL; print(err)},
    error = function(err) {msv_t1 <- NULL; print(err)}
  )
  # beep(1)
}
rm(vc,L,d,TT,n)

errTest <- function(j1, d) {
  if(j1 == 3) {
    stop(err = 'here i\'ll stop')
    # print('here i\'ll stop')
  } else {
    tmp <- length(d)
  }
  tmp
}

# intermediate saving ####
# mydir <- 'd:/stoch_vol_homestation/'
# li <- vector('list', 5)
for (i1 in 1:5) {
  tmp <- i1+rnorm(1)
  # li[[i1]] <- tmp
  assign(paste0('test_save_', i1), tmp)
  save(list=grep('test_save_', ls(), value = T), file = paste0(mydir, 'test_save_', i1, '.RData'))
  # if(i1==3) {
  #   stop('here i\'ll stop')
  # }
}


e1 <- new.env()
load(file = paste0(mydir, 'res.RData'), envir = e1)
ls(e1)
li <- get(x = 'est_msvt', envir = e1)
li

rm(tmp, i1, e1, li)

# diagnostics ####
stan_dens(est_msvt$sngs,separate_chains = T, 
          pars = c('H[1,1]', 'H[2,1]', 'y_gen[100,1]', 'y_gen[100,2]'))
stan_rhat(est_msvt$sngs, bins=35)

tmp <- get_posterior_mean(est_msvt_lkoh)[,3]
Sigma21 <- tmp[grep('Sigma\\[\\d{1,3},2,1\\]', names(tmp))] # covariance
Sigma22 <- tmp[grep('Sigma\\[\\d{1,3},2,2\\]', names(tmp))] # futures variance
futMean <- tmp[grep('y_gen\\[\\d{1,4},2\\]', names(tmp))] # futures mean

plot(Sigma21, type='l')
plot(Sigma22, type='l')
plot(futMean, type='l')
hist(futMean, 25)

rm(tmp, Sigma21, Sigma22, futMean)

# diagnostics
# rstan
msv_t1 <- est_msvt_chmf
options(digits=7)
summary(msv_t1, pars = c('H', 'y_gen[100,1]', 'y_gen[100,2]'), probs = c(0.025, 0.975), digits=3)
stan_plot(msv_t1, pars = c('H[1,1]', 'y_gen[100,1]', 'y_gen[100,2]'))
stan_trace(msv_t1, nrow=3, pars = c('Sigma[1,1,1]', 'H[2,1]', 'y_gen[100,1]', 'y_gen[100,2]'))
stan_dens(msv_t1,separate_chains = T, pars = c('H[1,1]', 'H[2,1]', 'y_gen[100,1]', 'y_gen[100,2]'))
stan_scat(msv_t1, pars = c('H[1,1]', 'H[1,2]'))
stan_diag(msv_t1)
stan_par(msv_t1, par = c('H[1,1]'))
stan_rhat(msv_t1, bins=35)
stan_ess(msv_t1)
stan_mcse(msv_t1)
stan_ac(msv_t1, pars = c('Sigma[10,2,1]', 'Sigma[10,2,2]', 'y_gen[100,2]'))
stan_diag(msv_t1,
          information = c("sample","stepsize", "treedepth","divergence")[3])

# coda diagnostics
require(coda)
coda.options(digits = 3) 
coda.options(default=TRUE)

msv_t1_coda <- As.mcmc.list(msv_t1, pars = c('Sigma[50,2,1]', 'Sigma[50,2,2]', 'y_gen[50,1]', 'y_gen[50,2]')) # convert to coda object 'y_gen[100,1]', 'y_gen[100,2]'
msv_t1_mcmc <- as.mcmc(msv_t1)

plot(msv_t1_coda)
summary(msv_t1_coda)
crosscorr.plot(msv_t1_coda)
betterPairs(as.matrix(msv_t1_coda[[1]]))
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
# ESS: if ESS << S (number of samples), it's bad (see Tsyplakov, 2007); n with respect to correlation
# ? inefficiency factors
# IACT (autocorrelation time, see Yu, Meyer), rho_k: high AC = slow mixing
# traceplots: should be as WN for good convergence + mixing (WN is good)
# if acceptance rate is high, it's bad (too small scale in tuning parameters); shoud be 20-70%
# density: not porcupine; if multi-model, then non-identifiability; plot each chain separately to see convergence
# ASDSF, AWTY (no in coda, see rwty package): convergence

