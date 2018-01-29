##
## File name: growth_rates.R
## Author: Fabio Silva
## Date created: 25/10/2017
## R version 3.4.2
## Description: This code computes the mean growth rate and associated 
## confidence envelope for the model that best-fits the SPD, over the 
## entire range of breakpoint values.
##

# Invoke Support Functions and Packages ---------------------------------------
source('./src.R')
require(minpack.lm)
require(doParallel)
require(foreach)


# Initialize Parallel Cluster -------------------------------------------------
nstreams <- detectCores()-1
cl <- makeCluster(nstreams, type = "FORK")
registerDoParallel(cl); getDoParWorkers()


# Load SPD Data ---------------------------------------------------------------
data <- read.csv('./Data/SPD_All_ConfidenceIntervals.csv', sep=';', header=F)
aux <- t(as.matrix(data))
data <- data.frame(x=-t(as.matrix(data))[,1], y=t(as.matrix(data))[,4])


# Choose Breakpoint Values ----------------------------------------------------
bp0 <- -14800
bp1 <- -12900
bp2 <- -10200

# Creates a range of 'extra' years on either side of breakpoint, every 'res' years
extra <- 150
res <- 10
br1 <- seq(bp1 - extra, bp1 + extra, res)
br2 <- seq(bp2 - extra, bp2 + extra, res)


# Re-Run Best Model and Save Growth Rates -------------------------------------
xxx <<- seq(bp0,-8000)

gr <- foreach (j = 1:NROW(br1), .combine=cbind) %dopar% {
  gr <- matrix(NA, 7+length(xxx), NROW(br2))
  for (k in 1:NROW(br2)) {
    aux2 <- matrix(NA, length(xxx), 1)
    
    ## Model E: Mixed Exp + Exp w/baseline + Exp
    start <- bp0; end <- br1[j]
    aux <- subset(data, x>=start & x<=end)
    modE1 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.00001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000))
    ind <- which(xxx >= start & xxx <= end); aux2[ind] <- (exp(coefficients(modE1)[1])-1)*100

    start <- br1[j]+1; end <- br2[k]
    aux <- subset(data, x>=start & x<=end)
    modE2 <- nlsLM(y~fexp(x, a, b, y0), data=aux, start=list(b=-0.000001, a=0.0001, y0=0), control=nls.lm.control(maxiter=1000, maxfev=5000))
    ind <- which(xxx >= start & xxx <= end)
    x <- seq(2,NROW(aux),1); yp <- predict(modE2, newdata = list(x=seq(start,end,1))); aux2[ind] <- c(NA,extractGrowthRate(x,yp))
    
    start <- br2[k]+1; end <- -8000
    aux <- subset(data, x>=start & x<=end)
    modE3 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.0001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000))
    ind <- which(xxx >= start & xxx <= end); aux2[ind] <- (exp(coefficients(modE3)[1])-1)*100
    
    ## Save regression coefficients
    gr[,k] <- c(coefficients(modE1), coefficients(modE2), coefficients(modE3), aux2)
  }
  gr
}
stopCluster(cl)
cc <- gr[1:7,]
gr <- gr[8:NROW(gr),]


# Extract Growth Rate and Pop Proxy Uncertainties -----------------------------
gr.phase1 <- matrix(NA, length(xxx),NCOL(gr))
gr.phase2 <- matrix(NA, length(xxx),NCOL(gr))
gr.phase3 <- matrix(NA, length(xxx),NCOL(gr))
pop.phase1 <- matrix(NA, length(xxx),NCOL(cc))
pop.phase2 <- matrix(NA, length(xxx),NCOL(cc))
pop.phase3 <- matrix(NA, length(xxx),NCOL(cc))

k <- 1
for (i in 1:NROW(br1)) {
  for (j in 1:NROW(br2)) {
    ind1 <- 1:which(xxx==br1[i])
    ind2 <- (which(xxx==br1[i]) + 1):which(xxx==br2[j])
    ind3 <- (which(xxx==br2[j]) + 1):NROW(xxx)
    
    ## Separate growth rates for the three phases
    gr.phase1[ind1,k] <- gr[ind1,k]
    gr.phase2[ind2,k] <- gr[ind2,k]
    gr.phase3[ind3,k] <- gr[ind3,k]
    
    ## Re-calculate fit predictions for different breakpoint values
    start <- bp0; end <- br1[i]; aux <- rev(subset(data, x>=start & x<=end)$x)
    pp <- cbind(cc[1,k], cc[2,k])
    pop.phase1[ind1,k] <- fexp(aux, pp[2], pp[1], 0)
    start <- br1[i]+1; end <- br2[j]; aux <- rev(subset(data, x>=start & x<=end)$x)
    pp <- cbind(cc[3,k], cc[4,k], cc[5,k])
    pop.phase2[ind2,k] <- fexp(aux, pp[2], pp[1], pp[3])
    start <- br2[j]+1; end <- -8000; aux <- rev(subset(data, x>=start & x<=end)$x)
    pp <- cbind(cc[6,k], cc[7,k])
    pop.phase3[ind3,k] <- fexp(aux, pp[2], pp[1], 0)
    
    k <- k + 1
  }
}



# Plot Means and Confidence Envelopes -----------------------------------------
## Regressions Fits
par(mfrow=c(2,1)); par(mar=c(4,4,2,2))
plot(data, type='l', xlim=c(-15000,-8000), xlab='Calendar Age (BP)', ylab='Density', main='', axes=F)
ll <- seq(18000,8000,-1000); for (i in seq(2,NROW(ll),2)) { ll[i] <- "" }
axis(1, at=seq(-18000,-8000,1000), labels=ll); axis(2); box()
polygon(c(min(br1),max(br1),max(br1),min(br1)), c(-1,-1,1,1), border=NA, col=MESS::col.alpha('grey',0.3))
polygon(c(min(br2),max(br2),max(br2),min(br2)), c(-1,-1,1,1), border=NA, col=MESS::col.alpha('grey',0.3))

plotCI(xxx, pop.phase1, 'blue', level=.95); text(mean(xxx[which(!is.na(pop.phase1))], na.rm=T), 0.000025, labels='Phase 1', col='blue')
plotCI(xxx, pop.phase2, 'red', level=.95); text(mean(xxx[which(!is.na(pop.phase2))], na.rm=T), 0.000025, labels='Phase 2', col='red')
plotCI(xxx, pop.phase3, 'darkgreen', level=.95); text(mean(xxx[which(!is.na(pop.phase3))], na.rm=T), 0.000025, labels='Phase 3', col='darkgreen')

## Growth Rates
plot(-100,-100,  xlim=c(-15000,-8000), xlab='Calendar Age (BP)', ylab='Annual Growth Rate', main='', ylim=c(-0.05,0.13), axes=F)
ll <- seq(18000,8000,-1000); for (i in seq(2,NROW(ll),2)) { ll[i] <- "" }
axis(1, at=seq(-18000,-8000,1000), labels=ll); axis(2, at=seq(-0.06,0.16,0.02)); box()
polygon(c(min(br1),max(br1),max(br1),min(br1)), c(-1,-1,1,1), border=NA, col=MESS::col.alpha('grey',0.3))
polygon(c(min(br2),max(br2),max(br2),min(br2)), c(-1,-1,1,1), border=NA, col=MESS::col.alpha('grey',0.3))
abline(h=0)

litGrowthRates() ## values from the literature

plotCI(xxx, gr.phase1, 'blue', level=.95); text(mean(xxx[which(!is.na(gr.phase1))], na.rm=T), 0.01, labels='Phase 1', col='blue')
plotCI(xxx, gr.phase2, 'red', level=.95); text(mean(xxx[which(!is.na(gr.phase2))], na.rm=T), 0.01, labels='Phase 2', col='red')
plotCI(xxx, gr.phase3, 'darkgreen', level=.95); text(mean(xxx[which(!is.na(gr.phase3))], na.rm=T), 0.01, labels='Phase 3', col='darkgreen')

dev.print(device=png, filename="Growth_Rates_uncertainty.png", width=12, height=14, units='in', res=300)

## Min and Max Growth Rates
grRates <- rbind( c(min(gr.phase1, na.rm=T), max(gr.phase1, na.rm=T)), c(min(gr.phase2, na.rm=T), max(gr.phase2, na.rm=T)) , c(min(gr.phase3, na.rm=T), max(gr.phase3, na.rm=T)))
rownames(grRates) <- c('Phase 1', 'Phase 2', 'Phase 3')
colnames(grRates) <- c('min','max')
grRates


## Save results for later use
save.image('growthRate_results.RData')


