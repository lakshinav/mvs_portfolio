# latex ####
require(xtable)
names <- colnames(ret_sc)


# write part of latex table to clipboard
write.table(paste('&', sprintf('%0.3f', tail(mv.go.sk[,"util", "5"],6))), 
            file='clipboard',
            quote = F, sep = ' ', row.names = F, col.names = F, eol = ' ')

# coeffs
gfit.garchcoef.p <- sapply(gfit@mfit$ufit@fit, function(x) 
{ x@fit$matcoef[,"Pr(>|t|)"]} ) #  "Pr(>|t|)"
gfit.sk.garchcoef.p <- sapply(gfit.sk@mfit$ufit@fit, function(x) 
{ x@fit$matcoef[,"Pr(>|t|)"]} ) # "Pr(>|t|)"
sum(gfit.garchcoef.p<0.05)
sum(gfit.sk.garchcoef.p<0.05)

gfit.garchcoef.se <- sapply(gfit@mfit$ufit@fit, function(x) 
{ x@fit$matcoef[," Std. Error"]} ) #  "Pr(>|t|)"
gfit.sk.garchcoef.se <- sapply(gfit.sk@mfit$ufit@fit, function(x) 
{ x@fit$matcoef[," Std. Error"]} ) # "Pr(>|t|)"

tab <- matrix(NA, n*2, 5*2)
colnames(tab) <- c(rownames(gfit.garchcoef.se), rownames(gfit.sk.garchcoef.se))

tab[2*(1:n)-1,1:nrow(gfit@mfit$garchcoef)] <- paste0(sprintf('%0.3f', t(gfit@mfit$garchcoef)),t(ifelse(gfit.garchcoef.p < 0.05, '*','')))
tab[2*(1:n),1:nrow(gfit@mfit$garchcoef)] <- paste0('(', sprintf('%0.3f', t(gfit.garchcoef.se)), ')')

tab[2*(1:n)-1,(nrow(gfit@mfit$garchcoef)+1):ncol(tab)] <- paste0(sprintf('%0.3f', t(gfit.sk@mfit$garchcoef)), t(ifelse(gfit.sk.garchcoef.p < 0.05, '*','')))
tab[2*(1:n),(nrow(gfit@mfit$garchcoef)+1):ncol(tab)] <- paste0('(', sprintf('%0.3f', t(gfit.sk.garchcoef.se)), ')')

# tab[2*((n+1):(2*n))-1,1:nrow(gfit.sk@mfit$garchcoef)] <- sprintf('%0.3f', t(gfit.sk@mfit$garchcoef))
# tab[2*((n+1):(2*n)),1:(nrow(gfit.sk@mfit$garchcoef)-2)] <- paste0('(', sprintf('%0.3f', t(gfit.sk.garchcoef.se)), ')')

na <- matrix(NA, n*2, 1)
na[2*(1:(n))-1] <- rep(paste0('\\multirow{2}{*}{\\textsc{', tolower(names), '}}'), 2)
tab <- cbind(na, tab)

tabp <- xtable(tab, label='tab:gogarch_coeffs', 
               caption = 'coeffs')
tt <- print(tabp, include.colnames = T, include.rownames = F, booktabs = T, 
            floating = T, print.results = F, sanitize.text.function = identity)
t1 <- gsub('[[:graph:][:space:]]*2016\\n', '', tt)
cat(t1, file='mvs_latex_output.txt', append = F)
file.edit('mvs_latex_output.txt')
rm(tab, tabp, na, tt, t1)

# average weights
# tab <- matrix(NA, n*2, length(ra))
# for (i1 in 1:length(ra)) {
#   tab[1:n,i1] <- sprintf('%0.3f', colMeans(w.go[,,which(ra[i1]==ra)]))
#   tab[(n+1):(2*n),i1] <- sprintf('%0.3f', colMeans(w.go.sk[,,which(ra[i1]==ra)]))
# }
# tab <- cbind(matrix(sapply(names, function(x) {paste0('\\textsc{', x, '}')}), 2*n,1), tab)
# tabp <- xtable(tab, label='tab:aver_weights', 
#                caption = 'Average optimal weights obtained from maximization of 2 and 3 moment utility')
# tt <- print(tabp, include.colnames = T, include.rownames = F, booktabs = T, 
#             floating = T, print.results = F, sanitize.text.function = identity)
# cat(tt, file='mvs_latex_output.txt', append = F)
# # check
# sum(as.numeric(tab[1:n, 2]))
# rm(i1, tab, tabp, tt)

tab <- matrix(NA, n, length(ra)*2)
colnames(tab) <- c(ra, ra)

avw <- apply(w.go[(isa+1):ts,,], c(2,3), mean)
tab[1:n,1:length(ra)] <- paste0(sprintf('%0.3f', t(gfit@mfit$garchcoef)),t(ifelse(gfit.garchcoef.p < 0.05, '*','')))
tab[2*(1:n),1:nrow(gfit@mfit$garchcoef)] <- paste0('(', sprintf('%0.3f', t(gfit.garchcoef.se)), ')')

tab[2*(1:n)-1,(nrow(gfit@mfit$garchcoef)+1):ncol(tab)] <- paste0(sprintf('%0.3f', t(gfit.sk@mfit$garchcoef)), t(ifelse(gfit.sk.garchcoef.p < 0.05, '*','')))
tab[2*(1:n),(nrow(gfit@mfit$garchcoef)+1):ncol(tab)] <- paste0('(', sprintf('%0.3f', t(gfit.sk.garchcoef.se)), ')')

na <- matrix(NA, n*2, 1)
na[2*(1:(n))-1] <- rep(paste0('\\multirow{2}{*}{\\textsc{', tolower(names), '}}'), 2)
tab <- cbind(na, tab)

tabp <- xtable(tab, label='tab:gogarch_coeffs', 
               caption = 'coeffs')
tt <- print(tabp, include.colnames = T, include.rownames = F, booktabs = T, 
            floating = T, print.results = F, sanitize.text.function = identity)
t1 <- gsub('[[:graph:][:space:]]*2016\\n', '', tt)
cat(t1, file='mvs_latex_output.txt', append = F)


# average moments
tab <- matrix(NA, 4*2, length(ra))
for (i1 in 1:length(ra)) {
  tab[1:2,i1] <- sprintf('%0.3f', colMeans(mv.go[,c("portf.mean", "portf.sd"),which(ra[i1]==ra)]))
  tab[3,i1] <- sprintf('%0.3f', mean(sk2mom[,which(ra[i1]==ra)], na.rm = T)*100)
  tab[4,i1] <- sprintf('%0.3f', mean(mv.go[,"util",which(ra[i1]==ra)]))
  
  tab[5:8,i1] <- sprintf('%0.3f', colMeans(mv.go.sk[,c("portf.mean", "portf.sd", "portf.sk", "util"),which(ra[i1]==ra)]*
                                             cbind(matrix(1, ts, 2), rep(100, ts), rep(1, ts))))
}
tab <- cbind(rep(c('Mean', 'Sd', 'Skew*100', 'Utility'), 2), tab)
tabp <- xtable(tab[,1:6], label='tab:aver_moments', 
               caption = 'Average moments obtained from maximization of 2 and 3 moment utility')
tt <- print(tabp, include.colnames = T, include.rownames = F, booktabs = T, 
            floating = T, print.results = F, sanitize.text.function = identity)
cat(tt, file='mvs_latex_output.txt', append = F)
rm(i1, tab, tabp, tt)

# descriptive statistics
tabp <- xtable(t(ds), label='tab:descr_stat', digits = 3,
               caption = 'Descriptive statistics for log returns')
tt <- print(tabp, include.colnames = T, include.rownames = T, booktabs = T, 
            floating = T, print.results = F, sanitize.text.function = identity)
cat(tt, file='mvs_latex_output.txt', append = F)
rm(tabp, tt)
file.show('mvs_latex_output.txt')

# risk measures for forecast
tab <- matrix(NA, 3*2, length(ra))
mp <- matrix(rep(c(rep(0.5/(isa-tech), isa-tech), rep(0.5/oos, oos)), 3), nrow = ts-tech)
for (i1 in 1:length(ra)) {
  tab[1:3,i1] <- sprintf('%0.3f', -colSums(riskcomp[,c("sd", "cvar", "cdar"),which(ra[i1]==ra)]*mp))
  tab[4:6,i1] <- sprintf('%0.3f', -colSums(riskcomp.sk[ts-tech,c("sd", "cvar", "cdar"),which(ra[i1]==ra)]*mp))
}
tab <- cbind(rep(c('SD', 'CVaR', 'CDaR'), 2), tab)
tabp <- xtable(tab, label='tab:risk_measures',
               caption = 'Mean values of risk measures obtained from maximization of 2 and 3 moment utility')
tt <- print(tabp, include.colnames = T, include.rownames = F, booktabs = T, 
            floating = T, print.results = F, sanitize.text.function = identity)
cat(tt, file='mvs_latex_output.txt', append = F)
rm(i1, tab, tabp, tt, mp)
file.show('mvs_latex_output.txt')

# out-of-sample table
tab <- matrix(NA, 5*2, oos)

tab[1,] <- sprintf('%0.3f', -riskcomp[(dim(riskcomp)[1]-oos+1):dim(riskcomp)[1], 'cvar', "5"])
tab[2,] <- sprintf('%0.3f', -riskcomp[(dim(riskcomp)[1]-oos+1):dim(riskcomp)[1], 'cdar', "5"])
tab[3,] <- sprintf('%0.3f', mv.go[(dim(mv.go)[1]-oos+1):dim(mv.go)[1], 'portf.mean', "5"])
tab[4,] <- sprintf('%0.3f', mv.go[(dim(mv.go)[1]-oos+1):dim(mv.go)[1], 'portf.sd', "5"])
tab[5,] <- sprintf('%0.3f', 100*sk2mom[(dim(sk2mom)[1]-oos+1):dim(sk2mom)[1], 2])

tab[6,] <- sprintf('%0.3f', -riskcomp.sk[(dim(riskcomp.sk)[1]-oos+1):dim(riskcomp.sk)[1], 'cvar', "5"])
tab[7,] <- sprintf('%0.3f', -riskcomp.sk[(dim(riskcomp.sk)[1]-oos+1):dim(riskcomp.sk)[1], 'cdar', "5"])
tab[8,] <- sprintf('%0.3f', mv.go.sk[(dim(mv.go.sk)[1]-oos+1):dim(mv.go.sk)[1], 'portf.mean', "5"])
tab[9,] <- sprintf('%0.3f', mv.go.sk[(dim(mv.go.sk)[1]-oos+1):dim(mv.go.sk)[1], 'portf.sd', "5"])
tab[10,] <- sprintf('%0.3f', 100*mv.go.sk[(dim(mv.go.sk)[1]-oos+1):dim(mv.go.sk)[1], 'portf.sk', "5"])

tab <- cbind(rep(c('CVaR', 'CDaR', 'Mean', 'SD', 'Skew*100'), 2), tab)
tabp <- xtable(tab, label='tab:oos_measures',
               caption = 'Risk measures and moments obtained from maximization of 2 and 3 moment utility for out-of-sample period')
tt <- print(tabp, include.colnames = T, include.rownames = F, booktabs = T, 
            floating = T, print.results = F, sanitize.text.function = identity)
cat(tt, file='mvs_latex_output.txt', append = F)
rm(tab, tabp, tt)
file.show('mvs_latex_output.txt')


# plots ####
d <- tail(index(ret), ts-tech) # dates
# plot(mv.go.sk[,'portf.sk', '5'], type='l')
# plot(mv.go.sk[,'util', '15'], type='l')
# plot(mv.go[,'util', '5'], type='l')

# average weights for ra=3...20
require(plyr)
require(ggplot2)
require(reshape2)

tab <- matrix(NA, n*2, length(ra))
for (i1 in 1:length(ra)) {
  tab[1:n,i1] <- colMeans(w.go[,,which(ra[i1]==ra)])
  tab[(n+1):(2*n),i1] <- colMeans(w.go.sk[,,which(ra[i1]==ra)])
}
rm(i1)

gg.w.go <- adply(tab[1:n, 1:5], .margins = c(2,1), .id = c('ra', 'ticker'), expand = F)
gg.w.go$ra <- as.numeric(gg.w.go$ra)
gg.w.go$ticker <- factor(gg.w.go$ticker, levels = 1:n, labels = names)

p1 <- ggplot(data=gg.w.go, aes(x=ra, y=V1, color=ticker))
  p1 + geom_line(size=1.2) + geom_point(show_guide = FALSE, size=3.5) + 
    ylab('Weight') + xlab('Risk aversion') + 
    scale_x_continuous(labels=ra[1:5]) +
    theme_bw()
  

# gogarch coeffs (no skew)
gfit.garchcoef.prob <- sapply(gfit@mfit$ufit@fit, function(x) 
{ x@fit$matcoef[,ncol(x@fit$matcoef)]} )

for (i1 in 1:nrow(gfit@mfit$garchcoef)) {
pl <- ggplot(data.frame(x=names, y=gfit@mfit$garchcoef[i1,], f=gfit.garchcoef.prob[i1,]), 
                  aes(x = x, y = y, fill=f<0.05)) + 
  geom_bar(stat = 'identity', color='gray') + 
  theme_bw() + 
  scale_x_discrete(name='') + 
  scale_y_continuous(name=rownames(gfit@mfit$garchcoef)[i1])
eval(parse(text=paste0('barcoef', i1, '<- pl'))) 
}
rm(i1, pl)
# source("F:\\Google Диск\\Nauka\\r_useful_staff\\multiplot.R")
multiplot(barcoef1, barcoef2, barcoef3, cols = nrow(gfit@mfit$garchcoef))

# gogarch coeffs (with skew)
gfit.sk.garchcoef.prob <- sapply(gfit.sk@mfit$ufit@fit, function(x) 
{ x@fit$matcoef[,ncol(x@fit$matcoef)]} )

for (i1 in 1:nrow(gfit.sk@mfit$garchcoef)) {
  pl <- ggplot(data.frame(x=names, y=gfit.sk@mfit$garchcoef[i1,], f=gfit.sk.garchcoef.prob[i1,]), 
               aes(x = x, y = y, fill=f<0.05)) + 
    geom_bar(stat = 'identity', color='gray') + 
    theme_bw() + 
    scale_x_discrete(name='') + 
    scale_y_continuous(name=rownames(gfit.sk@mfit$garchcoef)[i1])
  eval(parse(text=paste0('barcoef', i1, '<- pl'))) 
}
rm(i1, pl)
# source("F:\\Google ????\\Nauka\\r_useful_staff\\multiplot.R")
multiplot(barcoef1, barcoef2, barcoef3, barcoef4, barcoef5, cols = nrow(gfit.sk@mfit$garchcoef))

# arma coeffs
for (i1 in 1:nrow(gfit@mfit$arcoef)) {
  pl <- ggplot(data.frame(x=names, y=gfit@mfit$arcoef[i1,]), 
               aes(x = x, y = y)) + 
    geom_bar(stat = 'identity', color='gray') + 
    theme_bw() + 
    scale_x_discrete(name='') + 
    scale_y_continuous(name=rownames(gfit@mfit$arcoef)[i1])
  eval(parse(text=paste0('barcoef', i1, '<- pl'))) 
}
rm(i1, pl)
# source("F:\\Google ????\\Nauka\\r_useful_staff\\multiplot.R")
multiplot(barcoef1, barcoef2, barcoef3, cols = 3)
    

# confidence intervals for weights
# huita...
require(plyr)
require(ggplot2)
gg.w.go <- adply(w.go[,3,], .margins = c(1,2), .id = c('date', 'ra')) # gmkn
gg.w.go$ra <- factor(gg.w.go$ra, levels = 1:7, labels = ra)
# head(gg.w.go)

se <- ddply(gg.w.go, .(ra), summarize, sd=sd(V1))
v1mean <- ddply(gg.w.go, .(ra), summarize, mean=mean(V1))
res <- cbind(v1mean$mean - 1.96*se$sd, v1mean$mean + 1.96*se$sd)

limits <- aes(ymax = res[,2], ymin=res[,1])
c <- ggplot(data.frame(), aes(x=1:length(ra), y=v1mean$mean))
c +  geom_point(colour = "red", size = 5) +
geom_errorbar(limits, width=.2) +
scale_x_continuous(breaks=(1:length(ra)), labels=ra, name='Degree of risk aversion') +
scale_y_continuous(name='Average weight') + 
# ggtitle('SNGS')

boxplot(x = w.go[,3,])

tt <- apply(w.go[,3,], MARGIN = 2, median)

# utility
ra. <- '5' # risk aversion to plot
postscript(paste0("utility_", "ra_", ra., ".eps"), width=8, height=4, horizontal = F, 
           paper="special") # colormodel = 'gray', 
plotComp(x = mv.go[,'util', ra.], d = d, ts = ts, tech = tech, oos = oos, 
         ylim = c(-6, 0), pch=20,
         main = '', ylab='')
par(new=T)
plotComp(x = mv.go.sk[,'util', ra.], d = d, ts = ts, tech = tech, oos = oos, 
         ylim = c(-6, 0), pch=17,
         main = '', ylab='Utility')
par(new=F)
legend('bottomright', 
       c("no skew","skew"),
       col = c("black", "black"),
       pch = c(20, 17), 
       pt.bg = c("black", "black"), bty='n',
       text.font=2, pt.cex=1.5, 
       y.intersp=1, text.width=200)
dev.off()

# portfolio return
ra. <- '5' # risk aversion to plot
postscript(paste0("returns_", "ra_", ra., ".eps"), width=8, height=4, horizontal = F, 
           colormodel = 'gray', paper="special")
plotComp(x = mv.go[,'portf.mean', ra.], d = d, ts = ts, tech = tech, oos = oos, 
         ylim = c(-.01, .07), pch=20,
         main = '', ylab='')
par(new=T)
plotComp(x = mv.go.sk[,'portf.mean', ra.], d = d, ts = ts, tech = tech, oos = oos, 
         ylim = c(-.01, .07), pch=17,
         main = '', ylab='Return')
# par(new=T) #  -.3, .2
# plotComp(x = ret %*% rep(1/n, n), d = d, ts = ts, tech = tech, oos = oos, 
#          ylim = c(-.3, .2), pch=15,
#          main = '', ylab='Return')
par(new=F)
legend('topleft', 
       c("no skew","skew"),
       col = c("black", "black"),
       pch = c(20, 17), 
       pt.bg = c("black", "black"), bty='n',
       text.font=2, pt.cex=1.5,
       y.intersp=1, text.width=200)
dev.off()

# portfolio sd
# ra. <- '15' # risk aversion to plot
# postscript(paste0("sd_", "ra_", ra., ".eps"), width=8, height=4, horizontal = F, 
#            colormodel = 'gray', paper="special")
# plotComp(x = mv.go[,'portf.sd', ra.], d = d, ts = ts, tech = tech, oos = oos, 
#          ylim = c(.05, .12), pch=20,
#          main = '', ylab='')
# par(new=T)
# plotComp(x = mv.go.sk[,'portf.sd', ra.], d = d, ts = ts, tech = tech, oos = oos, 
#          ylim = c(.05, .12), pch=17,
#          main = '', ylab='Standard deviation')
# par(new=F)
# legend('topleft', 
#        c("no skew","skew"),
#        col = c("black", "black"),
#        pch = c(20, 17), 
#        pt.bg = c("black", "black"), bty='n',
#        text.font=2, pt.cex=1.5,
#        y.intersp=1, text.width=200)
# dev.off()

# portfolio skewness
ra. <- '3'
empsk <- skewpt(ret, w.go[,,which(as.numeric(ra.)==ra)]) # portfolio skewness in 2 mom optimization

postscript(paste0("skewness_", "ra_", ra., ".eps"), width=8, height=4, horizontal = F, 
           colormodel = 'gray', paper="special")
plotComp(x = empsk, d = d, ts = ts, tech = tech, oos = oos, 
         ylim = c(-.08, .08), pch=20,
         main = '', ylab='Skewness')
par(new=T)
plotComp(x = mv.go.sk[,'portf.sk', ra.], d = d, ts = ts, tech = tech, oos = oos, 
         ylim = c(-.08, .08), pch=17,
         main = '', ylab='Skewness')
par(new=F)
legend('bottomleft', 
       c("no skew","skew"),
       col = c("black", "black"),
       pch = c(20, 17), 
       pt.bg = c("black", "black"), bty='n',
       text.font=2, pt.cex=1.5,
       y.intersp=1, text.width=200)
dev.off()
rm(empsk)

# cvar
ra. <- 10
postscript(paste0("cvar_", "ra_", ra., ".eps"), width=8, height=4, horizontal = F, 
            paper="special") # colormodel = 'gray',
plotComp(x = c(rep(NA, tech), -riskcomp[,'cvar', as.character(ra.)]), 
         d = d, ts = ts, tech = tech, oos = oos, 
         ylim = c(-.3, -.1), pch=20,
         main = '', ylab='CVaR')
par(new=T)
plotComp(x = c(rep(NA, tech), -riskcomp.sk[,'cvar', as.character(ra.)]), 
         d = d, ts = ts, tech = tech, oos = oos,
         ylim = c(-.3, -.1), pch=17,
         main = '', ylab='CVaR')
par(new=F)
legend('topleft', 
       c("no skew","skew"),
       col = c("black", "black"),
       pch = c(20, 17), 
       pt.bg = c("black", "black"), bty='n',
       text.font=2, pt.cex=1.5,
       y.intersp=1, text.width=200)
dev.off()

# cdar
ra. <- '5'
postscript(paste0("cdar_", "ra_", ra., ".eps"), width=8, height=4, horizontal = F, 
           colormodel = 'gray', paper="special")
plotComp(x = c(rep(NA, tech), -riskcomp[1,,'cdar']), 
         d = d, ts = ts, tech = tech, oos = oos, 
         ylim = c(-1.5, -.5), pch=20,
         main = '', ylab='CDaR')
par(new=T)
plotComp(x = c(rep(NA, tech), -riskcomp[2,,'cdar']), 
         d = d, ts = ts, tech = tech, oos = oos,
         ylim = c(-1.5, -.5), pch=17,
         main = '', ylab='CDaR')
par(new=F)
legend('topright', 
       c("no skew","skew"),
       col = c("black", "black"),
       pch = c(20, 17), 
       pt.bg = c("black", "black"), bty='n',
       text.font=2, pt.cex=1.5,
       y.intersp=1, text.width=200)
dev.off()

# mad
# ra. <- '3'
# postscript(paste0("mad_", "ra_", ra., ".eps"), width=8, height=4, horizontal = F, 
#            colormodel = 'gray', paper="special")
# plotComp(x = c(rep(NA, tech), -riskcomp[1,,'mad']), 
#          d = d, ts = ts, tech = tech, oos = oos, 
#          ylim = c(-.1, -.05), pch=20,
#          main = '', ylab='MAD')
# par(new=T)
# plotComp(x = c(rep(NA, tech), -riskcomp[2,,'mad']), 
#          d = d, ts = ts, tech = tech, oos = oos,
#          ylim = c(-.1, -.05), pch=17,
#          main = '', ylab='MAD')
# par(new=F)
# legend('bottomright', 
#        c("no skew","skew"),
#        col = c("black", "black"),
#        pch = c(20, 17), 
#        pt.bg = c("black", "black"), bty='n',
#        text.font=2, pt.cex=1.5,
#        y.intersp=1, text.width=200)
# dev.off()

# risk measures
# riskcomp structure changed!!!!!
# no asymmetry
ra. <- '5'
pm <- ret %*% w.go[ts,,which(as.numeric(ra.)==ra)]
pp <- ecdf(sort(pm))
postscript(paste0("riskmeasure_", "ra_", ra., ".eps"), width=8, height=4, horizontal = F, 
           colormodel = 'gray', paper="special")
plot(sort(pm), pp(sort(pm)), type = "l", lty=1,
     col = "black", ylab = "Probability", xlab="Portfolio return", main = "",
     xlim=c(-.22,-.21), ylim=c(0, .06))

points(-riskcomp['no skew',ts-tech,'cvar'], pp(-riskcomp[1,ts-tech,'cvar']), col = "black", pch=20, cex=2) # 
par(new=T)
rm(pm)
# with asymmetry
pm <- ret %*% w.go.sk[ts,,which(as.numeric(ra.)==ra)]
pp <- ecdf(sort(pm))
plot(sort(pm), pp(sort(pm)), type = "l", lty=2,
     col = "black", ylab = "Probability", xlab="Portfolio return", main = "CVaR",
     xlim=c(-.22,-.21), ylim=c(0, .06)) 
points(-riskcomp['skew',ts-tech,'cvar'], pp(-riskcomp['skew',ts-tech,'cvar']), col = "black", pch=17, cex=1.8)
legend('topleft', 
       c("no skew","skew"),
       col = c("black", "black"),
       pch = c(20, 17), 
       pt.bg = c("black", "black"), bty='n',
       text.font=2, pt.cex=1.5,
       y.intersp=1, text.width=200)
par(new=F)
dev.off()
file.show(paste0("riskmeasure_", "ra_", ra., ".eps"))
rm(pp,pm)


rm(tech,d,ra.)
