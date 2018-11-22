##
## File name: plot_best_models.R
## Author: Fabio Silva
## Date created: 22/10/2018
## R version 3.5.0
## Description: This code takes teh results of model_selection.R to plot the 
## regression fits of the best breakpoint values for each model
##


# Invoke Support Functions and Packages ---------------------------------------
source('./src.R')
require(minpack.lm)


# Load SPD Data ---------------------------------------------------------------
data <- read.csv('./data/SPD_All_ConfidenceIntervals.csv', sep=';', header=F)
aux <- t(as.matrix(data))
data <- data.frame(x=-t(as.matrix(data))[,1], y=t(as.matrix(data))[,4])
data$x <- -data$x


# Choose or Import Parameters -------------------------------------------------
bp0 <- -16600
ll <- seq(18000,8000,-1000); for (i in seq(2,NROW(ll),2)) { ll[i] <- "" } ## labels for x-axis of plots
load('modSel_results.RData')


# Initialize Plots ------------------------------------------------------------
par(mfrow=c(2,3), mar=c(5,1,2,1))

## Model A: Single Exponential
plot(data, type='l', xlim=c(-18000,-8000), xlab='Calendar Age (BP)', ylab='', main='Model A', axes=F, cex.main=1.5); axis(1, at=seq(-18000,-8000,1000), labels=ll)
param <- modBestParam['model A',]
start <- param[1]; end <- -8000
aux <- subset(data, x>=start & x<=end)
modA <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.001, a=1), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modA, newdata = list(x=seq(start,end,1))), col='red', lwd=2)

legend('topleft', legend=c('Population Proxy', 'Model Fits with best-fitting thresholds'), col=c('black', 'red'), lty=c(1,1), lwd=c(1,2), bty='n') ## Legend


## Model B: Single Logistic
plot(data, type='l', xlim=c(-18000,-8000), xlab='Calendar Age (BP)', ylab='', main='Model B', axes=F, cex.main=1.5); axis(1, at=seq(-18000,-8000,1000), labels=ll)
param <- modBestParam['model B',]
start <- param[1]; end <- -8000
aux <- subset(data, x>=start & x<=end)
modB <- nlsLM(y~ SSlogisX(x, A, xmid, s, ysc), data=aux, start=list(A=0.00012, xmid=-10000, s=100, ysc=0.00001), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modB, newdata = list(x=seq(start,end,1))), col='red', lwd=2)


## Model C: Dual Exponential
plot(data, type='l', xlim=c(-18000,-8000), xlab='Calendar Age (BP)', ylab='', main='Model C', axes=F, cex.main=1.5); axis(1, at=seq(-18000,-8000,1000), labels=ll)
param <- modBestParam['model C',]
start <- param[1]; end <- param[2]
aux <- subset(data, x>=start & x<=end)
modC1 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.00001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modC1, newdata = list(x=seq(start,end,1))), col='red', lwd=2)

abline(v=param[2], col='grey', lty=2); text(param[2], 0, paste0(-param[2]," cal BP"), col='grey', pos=4, offset=0.5)
start <- param[2]+1; end <- -8000
aux <- subset(data, x>=start & x<=end)
modC2 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.001, a=1), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modC2, newdata = list(x=seq(start,end,1))), col='red', lwd=2)


## Model D: Mixed Exponential + Logistic
plot(data, type='l', xlim=c(-18000,-8000), xlab='Calendar Age (BP)', ylab='', main='Model D', axes=F, cex.main=1.5); axis(1, at=seq(-18000,-8000,1000), labels=ll)
param <- modBestParam['model D',]
start <- param[1]; end <- param[2]
aux <- subset(data, x>=start & x<=end)
modD1 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.00001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modD1, newdata = list(x=seq(start,end,1))), col='red', lwd=2)

abline(v=param[2], col='grey', lty=2); text(param[2], 0, paste0(-param[2]," cal BP"), col='grey', pos=4, offset=0.5)
start <- param[2]+1; end <- -8000
aux <- subset(data, x>=start & x<=end)
modD2 <- nlsLM(y~ SSlogisX(x, A, xmid, s, ysc), data=aux, start=list(A=0.00012, xmid=-10000, s=100, ysc=0.00001), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modD2, newdata = list(x=seq(start,end,1))), col='red', lwd=2)


## Model E: Mixed Exp + Exp w/baseline + Exp
plot(data, type='l', xlim=c(-18000,-8000), xlab='Calendar Age (BP)', ylab='', main='Model E', axes=F, cex.main=1.5); axis(1, at=seq(-18000,-8000,1000), labels=ll)
param <- modBestParam['model E',]
start <- param[1]; end <- param[2]
aux <- subset(data, x>=start & x<=end)
modE1 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.00001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modE1, newdata = list(x=seq(start,end,1))), col='red', lwd=2)

abline(v=param[2], col='grey', lty=2); text(param[2], 0, paste0(-param[2]," cal BP"), col='grey', pos=4, offset=0.5)
start <- param[2]+1; end <- param[3]
aux <- subset(data, x>=start & x<=end)
modE2 <- nlsLM(y~fexp(x, a, b, y0), data=aux, start=list(b=-0.0001, a=0.0001, y0=0), control=nls.lm.control(maxiter=1000, maxfev=5000), upper=c(0,Inf,Inf))
lines(seq(start,end,1), predict(modE2, newdata = list(x=seq(start,end,1))), col='red', lwd=2)

abline(v=param[3], col='grey', lty=2); text(param[3], 0, paste0(-param[3]," cal BP"), col='grey', pos=4, offset=0.5)
start <- param[3]+1; end <- -8000
aux <- subset(data, x>=start & x<=end)
modE3 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.0001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modE3, newdata = list(x=seq(start,end,1))), col='red', lwd=2)


## Model F: Mixed Exp + Exp w/baseline + Logistic
plot(data, type='l', xlim=c(-18000,-8000), xlab='Calendar Age (BP)', ylab='', main='Model F', axes=F, cex.main=1.5); axis(1, at=seq(-18000,-8000,1000), labels=ll)
param <- modBestParam['model F',]
start <- param[1]; end <- param[2]
aux <- subset(data, x>=start & x<=end)
modF1 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.00001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modF1, newdata = list(x=seq(start,end,1))), col='red', lwd=2)

abline(v=param[2], col='grey', lty=2); text(param[2], 0, paste0(-param[2]," cal BP"), col='grey', pos=4, offset=0.5)
start <- param[2]+1; end <- param[3]
aux <- subset(data, x>=start & x<=end)
modF2 <- nlsLM(y~fexp(x, a, b, y0), data=aux, start=list(b=-0.0001, a=0.0001, y0=0), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modF2, newdata = list(x=seq(start,end,1))), col='red', lwd=2)

abline(v=param[3], col='grey', lty=2); text(param[3], 0, paste0(-param[3]," cal BP"), col='grey', pos=4, offset=0.5)
start <- param[3]+1; end <- -8000
aux <- subset(data, x>=start & x<=end)
modF3 <- nlsLM(y~ SSlogisX(x, A, xmid, s, ysc), data=aux, start=list(A=0.00012, xmid=-10000, s=100, ysc=0.00001), control=nls.lm.control(maxiter=1000, maxfev=5000))
lines(seq(start,end,1), predict(modF3, newdata = list(x=seq(start,end,1))), col='red', lwd=2)




# Save Plot -------------------------------------------------------------------
dev.print(device=png, "plot_best_models.png", width=8/2*3, height=6, units='in', res=300)


