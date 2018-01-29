##
## File name: model_selection.R
## Author: Fabio Silva
## Date created: 25/10/2017
## R version 3.4.2
## Description: This code does regression fitting for the six models considered
## over a range of values of the breakpoints and uses the information-theoretic 
## approach of Burnham and Anderson (2002) to identify the model that best-fits
## the SPD data
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


# Plot SDP and Breakpoints ----------------------------------------------------
plot(data, type='l', xlim=c(-18000,-8000), xlab='Calendar Age (BP)', ylab='Density', main='')
abline(v = bp0, col = 'blue', lty=2)
abline(v = bp1, col = 'blue', lty=2)
abline(v = bp2, col = 'blue', lty=2)
abline(v = min(br1), col = 'grey', lty=2); abline(v = max(br1), col = 'grey', lty=2)
abline(v = min(br2), col = 'grey', lty=2); abline(v = max(br2), col = 'grey', lty=2)
dev.print(device=png, filename="SPD_breakpoints.png", width=12, height=6, units='in', res=300)


# Parallel Runs of Model Regression Fitting -----------------------------------
res <- foreach (j = 1:NROW(br1), .combine=rbind) %dopar% {
  out <- rep(NA, 17)
  for (k in 1:NROW(br2)) {
    
    ## Model A: Single Exponential
    start <- bp0; end <- -8000
    aux <- subset(data, x>=start & x<=end)
    modA <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.001, a=1), control=nls.lm.control(maxiter=1000, maxfev=5000))
    K <- 3; n <- NROW(aux); logL<- as.numeric(logLik(modA))
    AICc.A <- -2*logL + 2*K + 2*K*(n/(n-K-1)); BIC.A <- -2*logL + K * log(n)
    
    ## Model B: Dual Exponential
    start <- bp0; end <- br1[j]
    aux <- subset(data, x>=start & x<=end)
    modB1 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.00001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000))
    K <- 3; n <- NROW(aux); logL <- as.numeric(logLik(modB1))
    AICc.B1 <- -2*logL + 2*K + 2*K*(n/(n-K-1)); BIC.B1 <- -2*logL + K * log(n)
    
    start <- br1[j]+1; end <- -8000
    aux <- subset(data, x>=start & x<=end)
    modB2 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.001, a=1), control=nls.lm.control(maxiter=1000, maxfev=5000))
    K <- 3; n <- NROW(aux); logL <- as.numeric(logLik(modB2))
    AICc.B <- AICc.B1 - 2*logL + 2*K + 2*K*(n/(n-K-1)); BIC.B <- BIC.B1 - 2*logL + K * log(n)
    
    ## Model C: Mixed Exponential + Logistic
    modC1 <- modB1
    AICc.C <- AICc.B1; BIC.C <- BIC.B1
    
    start <- br1[j]+1; end <- -8000
    aux <- subset(data, x>=start & x<=end)
    modC2 <- nlsLM(y~ SSlogisX(x, A, xmid, s, ysc), data=aux, start=list(A=0.00012, xmid=-10000, s=100, ysc=0.00001), control=nls.lm.control(maxiter=1000, maxfev=5000))
    K <- 5; n <- NROW(aux); logL <- as.numeric(logLik(modC2))
    AICc.C <- AICc.C - 2*logL + 2*K + 2*K*(n/(n-K-1)); BIC.C <- BIC.C - 2*logL + K * log(n)
    
    ## Model D: Triple Exponential
    modD1 <- modB1
    AICc.D <- AICc.B1; BIC.D <- BIC.B1
    
    start <- br1[j]+1; end <- br2[k]
    aux <- subset(data, x>=start & x<=end)
    # tries different starting values
    modD2.a <- try(nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=-0.0001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000)))
    modD2.b <- try(nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=-0.001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000)))
    if (class(modD2.a)!='try-error') { modD2 <- modD2.a } else if (class(modD2.b)!='try-error') { modD2 <- modD2.b } else { stop() }
    K <- 3; n <- NROW(aux); logL <- as.numeric(logLik(modD2))
    AICc.D <- AICc.D - 2*logL + 2*K + 2*K*(n/(n-K-1)); BIC.D <- BIC.D - 2*logL + K * log(n)
    
    start <- br2[k]+1; end <- -8000
    aux <- subset(data, x>=start & x<=end)
    modD3 <- nlsLM(y~fexp(x, a, b, 0), data=aux, start=list(b=0.0001, a=0.001), control=nls.lm.control(maxiter=1000, maxfev=5000))
    K <- 3; n <- NROW(aux); logL <- as.numeric(logLik(modD3))
    AICc.D3 <- - 2*logL + 2*K + 2*K*(n/(n-K-1)); BIC.D3 <- - 2*logL + K * log(n)
    AICc.D <- AICc.D + AICc.D3; BIC.D <- BIC.D + BIC.D3
    
    ## Model E: Mixed Exp + Exp w/baseline + Exp
    modE1 <- modB1
    AICc.E <- AICc.B1; BIC.E <- BIC.B1
    
    start <- br1[j]+1; end <- br2[k]
    aux <- subset(data, x>=start & x<=end)
    modE2 <- nlsLM(y~fexp(x, a, b, y0), data=aux, start=list(b=-0.000001, a=0.0001, y0=0), control=nls.lm.control(maxiter=1000, maxfev=5000))
    K <- 4; n <- NROW(aux); logL <- as.numeric(logLik(modE2))
    AICc.E <- AICc.E - 2*logL + 2*K + 2*K*(n/(n-K-1)); BIC.E <- BIC.E - 2*logL + K * log(n)
    
    modE3 <- modD3
    AICc.E <- AICc.E + AICc.D3; BIC.E <- BIC.E + BIC.D3
    
    ## Model F: Mixed Exp + Log + Exp
    modF1 <- modB1
    AICc.F <- AICc.B1; BIC.F <- BIC.B1
    
    start <- br1[j]+1; end <- br2[k]
    aux <- subset(data, x>=start & x<=end)
    modF2 <- nlsLM(y~flog(x, a, b, y0), data=aux, start=list(b=-br1[j], a=-0.00001, y0=0), control=nls.lm.control(maxiter=1000, maxfev=5000))
    K <- 4; n <- NROW(aux); logL <- as.numeric(logLik(modF2))
    AICc.F <- AICc.F - 2*logL + 2*K + 2*K*(n/(n-K-1)); BIC.F <- BIC.F - 2*logL + K * log(n)
    
    modF3 <- modD3
    AICc.F <- AICc.F + AICc.D3; BIC.F <- BIC.F + BIC.D3
    
    
    ## Save Model Selection Indices
    mm <- c(BIC.A, BIC.B, BIC.C, BIC.D, BIC.E, BIC.F, AICc.A, AICc.B, AICc.C, AICc.D, AICc.E, AICc.F)
    modSel <- c(which.min(mm[1:6]), which.min(mm[7:12]))
    out <- rbind(out,c(bp0,br1[j],br2[k], mm, modSel))
  }
  out[-1,]
}
stopCluster(cl)
colnames(res) <- c('bp0', 'bp1', 'bp2', 'BIC.A', 'BIC.B', 'BIC.C', 'BIC.D', 'BIC.E', 'BIC.F', 'AICc.A', 'AICc.B', 'AICc.C', 'AICc.D', 'AICc.E', 'AICc.F', 'selBIC', 'selAIC')
res <- data.frame(res); res$selBIC <- c("A","B","C","D","E")[res$selBIC]; res$selAIC <- c("A","B","C","D","E")[res$selAIC]


# Quick Analysis of Results ---------------------------------------------------
pracma::strcmp(res$selBIC, res$selAIC) ## Check if AIC always picks same model as BIC

table(res$selAIC)/NROW(res)*100 ## Proportion of times different models are picked

## Best breakpoint values for each model
modBestParam <- rbind(as.numeric(res[which.min(as.matrix(res[,9])),c(1:3,4,10)]), as.numeric(res[which.min(as.matrix(res[,10])),c(1:3,5,11)]), 
                      as.numeric(res[which.min(as.matrix(res[,11])),c(1:3,6,12)]), as.numeric(res[which.min(as.matrix(res[,12])),c(1:3,7,13)]),
                      as.numeric(res[which.min(as.matrix(res[,13])),c(1:3,8,14)]), as.numeric(res[which.min(as.matrix(res[,14])),c(1:3,9,15)]))
rownames(modBestParam) <- c('model A', 'model B', 'model C', 'model D', 'model E', 'model F')
colnames(modBestParam) <- c('bp0','bp1','bp2','BIC','AIC')
modBestParam

## Very best  model and breakpoint values
ind <- arrayInd(which.min(as.matrix(res[,10:15])), dim(res[,10:15]))
res[ind[1],]

## Save results for later use
save.image('modSel_results.RData')


