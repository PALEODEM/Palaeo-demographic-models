##
## File name: src.R
## Author: Fabio Silva
## Date created: 25/10/2017
## R version 3.4.2
## Decription: Support Functions for Model Selection
##


# Exponential for regression fitting --------------------------------------
fexp = function(x, a, b, y0){
  return(a*exp(b*x)+y0)
}


# Logarithm for regression fitting --------------------------------------
flog <- function(x, a, b, y0){
  aux <- a*log(x+b)+y0
  return(aux)
}


# Self-Starter for Logistic Regression --------------------------------------
SSlogisX <- selfStart(~ Asym/(1 + exp((xmid1 - x)/scal1)) + yscal,
                      
                      function(mCall, data, LHS)
                      {
                        xy <- sortedXyData(mCall[["x"]], LHS, data)
                        if(nrow(xy) < 4) {
                          stop("Too few distinct x values to fit a logistic")
                        }
                        z <- xy[["y"]]
                        if (min(z) <= 0) { z <- z + 0.05 * max(z) } # avoid zeroes
                        z <- z/(1.05 * max(z))                      # scale to within unit height
                        xy[["z"]] <- log(z/(1 - z))                 # logit transformation
                        aux <- coef(lm(x ~ z, xy))
                        parameters(xy) <- list(xmid1 = aux[1], scal1 = aux[2], yscal = aux[3])
                        pars <- as.vector(coef(nls(y ~ 1/(1 + exp((xmid1 - x)/scal1)) + yscal,
                                                   data = xy, algorithm = "plinear")))
                        setNames(c(pars[4], pars[1], pars[2], pars[3]),
                                 mCall[c("Asym", "xmid1", "scal1", "yscal")])
                      }, c("Asym", "xmid1", "scal1", "yscal"))


# Function to extract growth rate from any curve --------------------------------------
extractGrowthRate = function(i, x) {
  return(100*(x[i]/x[i-1] - 1))
}


# Function to calcualte and plot a Confidence Envelope/Interval as a shaded polygon --------------------------------------
plotCI = function(x, y, col, level=.682, plot=T, output=F) {
  # Conf Envelope
  y.ci <- apply(t(y), 2, quantile, probs=c(.5-level[1]/2, .5+level[1]/2), na.rm=T)
  if (plot) {
    ind <- which(is.na(y.ci[1,]))
    if (length(ind)>0) { aux <- x[-ind]; y.ci2 <- y.ci[,-ind] } else { aux <- x; y.ci2 <- y.ci }
    x.p <- c(aux,rev(aux)); y.p <- c(y.ci2[1,], rev(y.ci2[2,]))
    polygon(x.p, y.p, border=NA, col=MESS::col.alpha(col, .3))
  }
  
  # Mean
  y.m <- apply(t(y), 2, mean, na.rm=T)
  if (plot) { lines(x, y.m, col=col,lwd=2) }
  
  if (output) {
    df <- data.frame(x=x, mean = y.m, CI.bot = y.ci[1,], CI.top = y.ci[2,])
    return(df)
  }
}


# Function to plot Annual Growth Rates from the Literature ----------------
litGrowthRates = function(){
  ## Silva and Vander Linden 2017
  Silva_WMed <- list(x=c(-12000,-5000), y=c(0.04599-0.00051, 0.04599+0.00051))  
  xx <- c(Silva_WMed$x[1], Silva_WMed$x[1], Silva_WMed$x[2], Silva_WMed$x[2]); yy <- c(Silva_WMed$y[1], Silva_WMed$y[2], Silva_WMed$y[2], Silva_WMed$y[1])
  polygon(xx, yy, col='grey50', border=NA)
  text(min(Silva_WMed$x),max(Silva_WMed$y)+0.002,labels="Western Med Mesolithic", cex=0.8, col="black", pos=4)
  
  ## Goldberg et al 2016
  Goldberg <- list(x=c(-14000,-5500), y=c(0.131, 0.131))  
  lines(Goldberg$x, Goldberg$y, lwd=1, lty=3)
  text(-8000,min(Goldberg$y)+0.002,labels="South America", cex=0.8, col="black", pos=2)
  
  ## Zahid et al 2015
  Zahid <- list(x=c(-13000,-6000), y=c(0.041-0.003,0.041+0.003))  
  xx <- c(Zahid$x[1], Zahid$x[1], Zahid$x[2], Zahid$x[2]); yy <- c(Zahid$y[1], Zahid$y[2], Zahid$y[2], Zahid$y[1])
  polygon(xx, yy, col='black', border=1, density = 10)
  text(-12700,mean(Zahid$y),labels="Wyoming & Colorado", cex=0.8, col="black", pos=4)
}

