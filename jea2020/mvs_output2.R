require(ggplot2)
require(reshape2)
require(gridExtra)
require(cowplot) # one legend for grid plot

plot.box <- function(m,r) {
  # r - ra vector
  # m - risk measure
  # browser()
  p <- vector(mode = 'list', length(r))
  for (i1 in 1:length(r)) {
  rm_df <- data.frame(nosk=tail(riskcomp[,m,i1],oos), 
                      sk=tail(riskcomp.sk[,m,i1],oos),time=1:oos)
  rm_df <- melt(rm_df, measure.vars = 1:2)
  p[[i1]] <- ggplot(data = rm_df, aes(y=value, x=variable, col=variable)) + geom_boxplot() +
    ggtitle(paste('measure', m, 'ra', r[i1])) + theme_bw()
  }
  grid.arrange(grobs = p)
}
plot.box(2,ra)

plot.lines <- function(m,r) {
  # r - ra vector
  # m - risk measure
  # browser()
  p <- vector(mode = 'list', length(r))
  for (i1 in 1:length(r)) {
  rm_df <- data.frame(nosk=tail(riskcomp[,m,i1],oos), 
                      sk=tail(riskcomp.sk[,m,i1],oos),
                      time=as.Date(tail(rownames(ret6),oos)))
  rm_df <- melt(rm_df, measure.vars = 1:2)
  p[[i1]] <- ggplot(data = rm_df, aes(x=time, y=value, col=variable)) + geom_line() +
    ggtitle(paste('Risk aversion =', r[i1])) + 
    xlab('') + ylab('') + 
    theme_bw() + theme(legend.position="none") #+
    # scale_x_date(limits = as.Date(range(tail(rownames(ret6),oos))))
  }
  legend <- cowplot::get_legend(p[[1]] + 
                        scale_color_discrete(name='Utility:', labels = c("2 moment", "3 moment")) +
                        theme(legend.position = "bottom"))
  p_grid <- cowplot::plot_grid(plotlist = p, ncol = 2)
  cowplot::plot_grid(p_grid, legend, ncol = 1, rel_heights = c(1.1, 0.1))
  # grid.arrange(grobs = p)
  # print(p)
}

setEPS()
postscript(file='sd_grid.eps')
plot.lines(2,ra) # ok
dev.off()

setEPS()
postscript(file='mad_grid.eps')
plot.lines(1,ra) # ok
dev.off()

setEPS()
postscript(file='cvar_grid.eps')
plot.lines(3,ra) # ok
dev.off()

plot.lines(4,ra)
plot.lines(5,ra)
dimnames(riskcomp)[[2]]

rm_df <- melt(data.frame(noskew=tail(unlist(lapply(putil.go, parmareward)),oos),
                         skew=tail(unlist(lapply(putil.go.sk, parmareward)),oos),
                         time=1:oos), measure.vars = 1:2)
ggplot(data = rm_df, aes(x=time, y=value, col=variable)) + 
  geom_line() + theme_bw()


plot.lines_perf <- function(m,r) {
  # r - ra vector
  # m - perf measure
  # browser()
  p <- vector(mode = 'list', length(r))
  for (i1 in 1:length(r)) {
    rm_df <- data.frame(nosk=tail(mv.go[,m,i1],oos), 
                        sk=tail(mv.go.sk[,m,i1],oos),
                        time=as.Date(tail(rownames(ret6),oos)))
    rm_df <- melt(rm_df, measure.vars = 1:2)
    p[[i1]] <- ggplot(data = rm_df, aes(x=time, y=value, col=variable)) + geom_line() +
      ggtitle(paste('Risk aversion =', r[i1])) + # dimnames(mv.go)[[2]][m], ',',
      xlab('') + ylab('') + 
      theme_bw() + theme(legend.position="none")
  }
  legend <- cowplot::get_legend(p[[1]] + 
                                  scale_color_discrete(name='Utility:', labels = c("2 moment", "3 moment")) +
                                  theme(legend.position = "bottom"))
  p_grid <- cowplot::plot_grid(plotlist = p, ncol = 2)
  cowplot::plot_grid(p_grid, legend, ncol = 1, rel_heights = c(1.1, 0.1))
  # grid.arrange(grobs = p)
  # print(p)
}
dimnames(mv.go)[[2]]
dimnames(mv.go.sk)[[2]]

setEPS()
postscript(file='reward_grid.eps')
plot.lines_perf(1,ra)
dev.off()

plot.lines_perf('util',ra)
plot.lines_perf('status',ra)

plot.boxplot_perf <- function(m,r) {
  # r - ra vector
  # m - perf measure
  p <- vector(mode = 'list', length(r))
  for (i1 in 1:length(r)) {
    rm_df <- data.frame(nosk=tail(mv.go[,m,i1],oos), 
                        sk=tail(mv.go.sk[,m,i1],oos),
                        time=as.Date(tail(rownames(ret6),oos)))
    rm_df <- melt(rm_df, measure.vars = 1:2)
    p[[i1]] <- ggplot(data = rm_df, aes(x=time, y=value, col=variable)) + geom_boxplot() +
      ggtitle(paste(dimnames(mv.go)[[2]][m], ',', 'Risk aversion =', r[i1])) + 
      xlab('') + ylab('') + 
      theme_bw() + theme(legend.position="none")
  }
  legend <- cowplot::get_legend(p[[1]] + 
                                  scale_color_discrete(name='Utility:', labels = c("2 moment", "3 moment")) +
                                  theme(legend.position = "bottom"))
  p_grid <- cowplot::plot_grid(plotlist = p, ncol = 2)
  cowplot::plot_grid(p_grid, legend, ncol = 1, rel_heights = c(1.1, 0.1))
}
plot.boxplot_perf('util',ra)

avw <- apply(w.go[(isa+1):ts,,], c(2,3), mean)
avw.sk <- apply(w.go.sk[(isa+1):ts,,], c(2,3), mean)
colSums(avw)
colSums(avw.sk)
t1 <- ds[8,apply(avw.sk, 1, function(z) {z[1] < z[6]})] # skewness for increasing weights
t2 <- ds[8,apply(avw.sk, 1, function(z) {z[1] > z[6]})] # skewness for decreasing weights
t1 <- sort(t1)
t2 <- sort(t2)
plot(1:length(t1), t1, type='l', xlim=c(1,11), ylim=c(-2,3))
par(new=T)
plot(1:length(t2), t2, type='l', col='red', xlim=c(1,11), ylim=c(-2,3))
par(new=F)
rm(t1,t2)

# variance is almost the same across assets
t1 <- apply(w.go[(isa+1):ts,,], c(2), sd)
t2 <- apply(w.go.sk[(isa+1):ts,,], c(2), sd)
plot(1:length(t1), t1, type='l', xlim=c(1,20), ylim=c(0,1))
par(new=T)
plot(1:length(t2), t2, type='l', col='red', xlim=c(1,20), ylim=c(0,1))
par(new=F)
which.min(t2)
rm(t1,t2)

t3 <- sapply(gfit.sk@mfit$ufit@fit, function(x)
{ x@fit$matcoef[,' Estimate']} ) # skew parameter in go-garch
t3 <- t3['skew',]
names(t3) <- nam
avw <- avw[order(t3, decreasing = T),]
avw.sk <- avw.sk[order(t3, decreasing = T),]
t3 <- sort(t3, decreasing = T)

plot.lines_wei <- function() {
  # 1 plot for each asset
  # r - ra vector
  # m - risk measure
  # browser()
  p <- vector(mode = 'list', n)
  for (i1 in 1:n) {
    rm_df <- data.frame(noskew=avw[i1,], 
                        skew=avw.sk[i1,],
                        ra=ra)
    rm_df <- melt(rm_df, measure.vars = 1:2, value.name = 'weight')
    p[[i1]] <- ggplot(data = rm_df, aes(x=ra, y=weight, col=variable)) + geom_line() +
      labs(title=gsub('\\.ME', '', t3)[i1], x='', y='', color = '') +
      theme_bw() + theme(legend.position="none")
  }
  # p_leg <- ggplot(data = rm_df, aes(x=ra, y=weight, col=variable)) + geom_line() +
  #   ggtitle(gsub('\\.ME', '', nam)[i1]) + xlab('') + ylab('') + theme_bw() +
    
  legend <- cowplot::get_legend(p[[1]] + 
                                  scale_color_discrete(name='Utility:', labels = c("2 moment", "3 moment")) +
                                  theme(legend.position = "bottom"))
  p_grid <- cowplot::plot_grid(plotlist = p, ncol = 4)
  cowplot::plot_grid(p_grid, legend, ncol = 1, rel_heights = c(1.2, 0.1))
    # grid.arrange(grobs = p)
  # print(p)
}
setEPS()
postscript(file = 'weights_grid.eps')
plot.lines_wei()
dev.off()

plot.boxplot_rm <- function(m) {
  # 1 plot for all risk measures
  # m - risk measure

    rm_df.nsk <- data.frame(riskcomp[(isa+1):ts,m,]) 
    colnames(rm_df.nsk) <- ra
    rm_df.nsk <- melt(rm_df.nsk, measure.vars = 1:6, value.name = 'rm')

    rm_df.sk <- data.frame(riskcomp.sk[(isa+1):ts,m,]) 
    colnames(rm_df.sk) <- ra
    rm_df.sk <- melt(rm_df.sk, measure.vars = 1:6, value.name = 'rm')

    rm_df <- merge(rm_df.nsk, rm_df.sk, by = 'variable')
    colnames(rm_df) <- c('ra', 'noskew', 'skew')
    rm_df <- melt(rm_df, measure.vars = 2:3, value.name = 'rm')
    
    p <- ggplot(data = rm_df, aes(x=ra, y=rm, fill=variable)) + 
      geom_boxplot() +
      # labs(title=paste('measure', m)) +
      labs(title='', x='', y='') +
      theme_bw() + theme(legend.position="bottom") + 
      scale_fill_discrete(name='Utility:', labels = c("2 moment", "3 moment"))
    print(p)
}
plot.boxplot_rm(1)

plot.boxplot_perf <- function(m) {
  # 1 plot for all risk measures
  # m - risk measure
  browser()
  rm_df.nsk <- data.frame(mv.go[(isa+1):ts,m,]) 
  colnames(rm_df.nsk) <- ra
  rm_df.nsk <- melt(rm_df.nsk, measure.vars = 1:6, value.name = 'perf')
  
  rm_df.sk <- data.frame(mv.go.sk[(isa+1):ts,m,]) 
  colnames(rm_df.sk) <- ra
  rm_df.sk <- melt(rm_df.sk, measure.vars = 1:6, value.name = 'perf')
  
  rm_df <- merge(rm_df.nsk, rm_df.sk, by = 'variable')
  colnames(rm_df) <- c('ra', 'noskew', 'skew')
  rm_df <- melt(rm_df, measure.vars = 2:3, value.name = 'perf')
  
  p <- ggplot(data = rm_df, aes(x=ra, y=perf, color=variable)) + 
    geom_boxplot() +
    labs(title=paste('measure', m)) +
    # labs(title='', x='', y='') +
    theme_bw() + theme(legend.position="bottom") + 
    scale_fill_discrete(name='Utility:', labels = c("2 moment", "3 moment"))
  print(p)
}
plot.boxplot_perf('util')

setEPS()
postscript(file = 'sd_con_box.eps')
plot.boxplot_rm(2)
dev.off()

setEPS()
postscript(file = 'mad_con_box.eps')
plot.boxplot_rm(1)
dev.off()

setEPS()
postscript(file = 'cvar_con_box.eps')
plot.boxplot_rm(3)
dev.off()

plot.boxplot_rm(3)
plot.boxplot_rm(4) # no
plot.boxplot_rm(5) # no
dimnames(riskcomp)[[2]]

postscript(file = 'rm.ps')
plot.lines_rm()
dev.off()

plot.lines_wei2 <- function() {
  # plot 2 panels for sk and nosk (not very informative)
  # r - ra vector
  # m - risk measure
  # browser()
  p <- vector(mode = 'list', 2)
    rm_df <- melt(avw) %>% 
      transmute(value = value, ra = Var2, asset = as.factor(Var1))
    p[[1]] <- ggplot(data = rm_df, aes(x=ra, y=value, col=asset)) + geom_line(size=0.8)  + theme_bw()
    rm_df <- melt(avw.sk) %>% 
      transmute(value = value, ra = Var2, asset = as.factor(Var1))
    p[[2]] <- ggplot(data = rm_df, aes(x=ra, y=value, col=asset)) + geom_line(size=0.8)  + theme_bw()
  grid.arrange(grobs = p, ncol=2)
  # print(p)
}
plot.lines_wei2()


plot.hist <- function(m,r) {
  # r - ra
  # m - risk measure
  # browser()
  lwd <- 0.9
  rm_df <- data.frame(nosk=tail(riskcomp[,m,r],oos), 
                      sk=tail(riskcomp.sk[,m,r],oos),time=1:oos)
  rm_df <- melt(rm_df, measure.vars = 1:2)
  mu <- group_by(.data = rm_df, variable) %>% 
    summarise(mean_rm = mean(value))
  p <- ggplot(data = rm_df, aes(x=value, col=variable)) + geom_density(size=lwd) +
    geom_vline(data=mu, aes(xintercept=mean_rm, color=variable),
               linetype="dashed", size=lwd) +
    ggtitle(paste('measure', m, 'ra', r)) +
    theme_bw()
  print(p)
}

plot.weights <- function(a, r) {
  # a - assets' vector
  # r - ra
  # browser()
  lwd <- 0.5
  # rm_df <- data.frame(nosk=tail(riskcomp[,m,r],oos), 
  #                     sk=tail(riskcomp.sk[,m,r],oos),time=1:oos)
  p <- vector(mode = 'list', length(a))
  for (i1 in 1:length(a)) {
  wei_df <- data.frame(nosk=w.go[(isa+1):ts,a[i1],r],
                      sk=w.go.sk[(isa+1):ts,a[i1],r],
                      time = 1:oos)
  wei_df <- melt(wei_df, measure.vars = 1:2)
  # mu <- group_by(.data = rm_df, variable) %>% 
  #   summarise(mean_rm = mean(value))
  p[[i1]] <- ggplot(data = wei_df, aes(x=time, y=value, col=variable)) + geom_line(size=lwd) +
    ggtitle(paste('asset', a[i1], 'ra', r)) +
    theme_bw()
  }
  grid.arrange(grobs = p, ncol=4)
  # print(p)
}

# latex ####
# cumulative returns
tab <- matrix(NA,2,length(ra))
tab[1,] <- colSums(tail(mv.go[,1,],oos))
tab[2,] <- colSums(tail(mv.go.sk[,1,],oos))

tabp <- xtable(tab, label='tab:cum_ret', 
               caption = 'Cumulative portfolio returns obtained from maximization of 2 and 3 moment utility ')
tt <- print(tabp, include.colnames = T, include.rownames = F, booktabs = T, 
            floating = T, print.results = F, sanitize.text.function = identity)
t1 <- gsub('[[:graph:][:space:]]*2016\\n', '', tt)
cat(t1, file='mvs_latex_output.txt', append = F)

# utility
tab <- matrix(NA,2,length(ra))
# tab[1,] <- colMeans(tail(mv.go[,'util',],oos))
# tab[2,] <- colMeans(tail(mv.go.sk[,'util',],oos))
tab[1,] <- apply(mv.go[(isa+1):ts,'util',],2,function(z) {quantile(z, 0.25)}) # mean(z) 
tab[2,] <- apply(mv.go.sk[(isa+1):ts,'util',],2,function(z) {quantile(z, 0.25)})

tabp <- xtable(tab, label='tab:util', digits = 4,
               caption = 'Investors utility obtained from maximization of 2 and 3 moment utility')
tt <- print(tabp, include.colnames = T, include.rownames = F, booktabs = T, 
            floating = T, print.results = F, sanitize.text.function = identity)
t1 <- gsub('[[:graph:][:space:]]*2016\\n', '', tt)
cat(t1, file='mvs_latex_output.txt', append = F)

# weights analysis
nam <- colnames(ret6)
t1 <- apply(w.go.sk[(isa+1):ts,,],2,function(z) {colMeans(z)})[6,] # weights for ra=10
names(t1) <- nam
sort(t1)
t2 <- sapply(gfit.sk@mfit$ufit@fit, function(x)
{ x@fit$matcoef[,"Pr(>|t|)"]} ) # signif of skew in go-garch
t2 <- t2['skew',]
names(t2) <- nam
sort(t2, decreasing = T)
t3 <- sapply(gfit.sk@mfit$ufit@fit, function(x)
{ x@fit$matcoef[,' Estimate']} ) # skew in go-garch
t3 <- t3['skew',]
names(t3) <- nam
sort(t3)

sum(t1[which(t2<0.05)])
t1[names(t1) %in% names(t3)[1:10]]
# => nothing here

nam <- gsub('\\.ME', '', nam)
avw['BSPB' == nam,]
avw.sk['BSPB' == nam,]

avw['CNTL' == nam,]
avw.sk['CNTL' == nam,]

ass <- 'CNTL' #'URKA'  # 'BSPB' # 
gg_df <- melt(w.go[(isa+1):ts, ass == nam,])
gg_df$Var2 <- as.factor(gg_df$Var2)
ggplot(data=gg_df, aes(x=Var1, y=value, col=Var2)) + geom_line() + 
  ggtitle(ass) + scale_color_discrete(labels=ra)
rm(gg_df, ass)
