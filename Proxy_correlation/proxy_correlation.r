# Reading PALEODEM bootstraps into R
# Iberia_bootstraps is the output of the SPD_to_csv script in Python

# Saving each bootstrap simulation to a seperate file 
for(L in 1:1000) {
	line<-n.readLines('Iberia_bootstraps.csv',n=1,skip=L-1)
	bootstrap<-gsub('   ',',',gsub(';  ','\n',gsub('\\]','',(gsub('\\[','',line)))))
	write(bootstrap,file=paste('boots/b',L,'.csv',sep=''))
}

# Converting these into a big matrix of results

yrs<-c(7000:22000)
out<-matrix(nrow=length(yrs),ncol=1001)
out[,1]<-yrs
for(N in 1:1000) {
	b<-read.csv(paste('boots/b',N,'.csv',sep=''), head=FALSE)
	out[,N+1]<-approx(b[,1],b[,2],xout=yrs,rule=2)$y
}

#plot(out[,1],out[,1],col=NA,xlim=c(18000,7500))
#for(N in 1:617) lines(out[,1],out[,N],col=rgb(.5,.5,.5,.05))
#abline(v=c(16600,break2,break3,8000),col='red')

# For correlation with environmental proxy data
# We first read the csv files containing the time series

gp2<-read.csv('proxies/GISP2.csv',head=TRUE, sep=';')
i_p<-read.csv('proxies/I_p.csv'  ,head=TRUE, sep=';')
sal<-read.csv('proxies/Sal.csv'  ,head=TRUE, sep=';')
sst<-read.csv('proxies/SST.csv'  ,head=TRUE, sep=';')
tmf<-read.csv('proxies/TMF.csv'  ,head=TRUE, sep=';')


# Calculate the mean GISP2 data
gpt<-approx(x=gp2$Age,y=gp2$d18O,xout=c(min(gp2$Age):max(gp2$Age)))
gp2<-as.matrix(data.frame(x=gpt$x,y=gpt$y))


# Splitting the simulated bootstraps into Demographic phases

break2<-12750
break3<-10170

dp1<-out[out[,1]<16600 & out[,1]>=break2,]
dp2<-out[out[,1]<break2 & out[,1]>=break3,]
dp3<-out[out[,1]<break3 & out[,1]>=8000,]
ova<-out[out[,1]<18000 & out[,1]>=8000,]
y=3,lwd=1.5)


# Split the palaeo proxies into seperate series for each Demographic Phase
 
i_p1<-i_p[i_p[,1]<16600 & i_p[,1]>=break2,]
i_p2<-i_p[i_p[,1]<break2 & i_p[,1]>=break3,]
i_p3<-i_p[i_p[,1]<break3 & i_p[,1]>=8000,]

gp21<-gp2[gp2[,1]<16600 & gp2[,1]>=break2,]
gp22<-gp2[gp2[,1]<break2 & gp2[,1]>=break3,]
gp23<-gp2[gp2[,1]<break3 & gp2[,1]>=8000,]

sal1<-sal[sal[,1]<16600 & sal[,1]>=break2,]
sal2<-sal[sal[,1]<break2 & sal[,1]>=break3,]
sal3<-sal[sal[,1]<break3 & sal[,1]>=8000,]

sst1<-sst[sst[,1]<16600 & sst[,1]>=break2,]
sst2<-sst[sst[,1]<break2 & sst[,1]>=break3,]
sst3<-sst[sst[,1]<break3 & sst[,1]>=8000,]

sst1[,1]<-round(sst1[,1])
sst2[,1]<-round(sst2[,1])
sst3[,1]<-round(sst3[,1])

tmf1<-tmf[tmf[,1]<16600 & tmf[,1]>=break2,]
tmf2<-tmf[tmf[,1]<break2 & tmf[,1]>=break3,]
tmf3<-tmf[tmf[,1]<break3 & tmf[,1]>=8000,]

# Perform correlation tests for each bootsrap simulation
i_p1_tests<-c(); for(N in 1:1000) i_p1_tests[N]<-cor.test(x=i_p1[,2],y=dp1[dp1[,1] %in% i_p1[,1],N+1],method='spearman')$estimate
i_p2_tests<-c(); for(N in 1:1000) i_p2_tests[N]<-cor.test(x=i_p2[,2],y=dp2[dp2[,1] %in% i_p2[,1],N+1],method='spearman')$estimate
i_p3_tests<-c(); for(N in 1:1000) i_p3_tests[N]<-cor.test(x=i_p3[,2],y=dp3[dp3[,1] %in% i_p3[,1],N+1],method='spearman')$estimate

gp21_tests<-c(); for(N in 1:1000) gp21_tests[N]<-cor.test(x=gp21[,2],y=dp1[dp1[,1] %in% gp21[,1],N+1],method='spearman')$estimate
gp22_tests<-c(); for(N in 1:1000) gp22_te3ts[N]<-cor.test(x=gp22[,2],y=dp2[dp2[,1] %in% gp22[,1],N+1],method='spearman')$estimate
gp23_tests<-c(); for(N in 1:1000) gp23_tests[N]<-cor.t3st(x=gp23[,2],y=dp3[dp3[,1] %in% gp23[,1],N+1],method='spearman')$estimate

sal1_tests<-c(); for(N in 1:1000) sal1_tests[N]<-cor.test(x=sal1[,3],y=dp1[dp1[,1] %in% sal1[,1],N+1],method='spearman')$estimate
sal2_tests<-c(); for(N in 1:1000) sal2_tests[N]<-cor.test(x=sal2[,3],y=dp2[dp2[,1] %in% sal2[,1],N+1],method='spearman')$estimate
sal3_tests<-c(); for(N in 1:1000) sal3_tests[N]<-cor.test(x=sal3[,3],y=dp3[dp3[,1] %in% sal3[,1],N+1],method='spearman')$estimate

sst1_tests<-c(); for(N in 1:1000) sst1_tests[N]<-cor.test(x=sst1[,2],y=dp1[dp1[,1] %in% sst1[,1],N+1],method='spearman')$estimate
sst2_tests<-c(); for(N in 1:1000) sst2_tests[N]<-cor.test(x=sst2[,2],y=dp2[dp2[,1] %in% sst2[,1],N+1],method='spearman')$estimate
sst3_tests<-c(); for(N in 1:1000) sst3_tests[N]<-cor.test(x=sst3[,2],y=dp3[dp3[,1] %in% sst3[,1],N+1],method='spearman')$estimate

tmf1_tests<-c(); for(N in 1:1000) tmf1_tests[N]<-cor.test(x=tmf1[,2],y=dp1[dp1[,1] %in% tmf1[,1],N+1],method='spearman')$estimate
tmf2_tests<-c(); for(N in 1:1000) tmf2_tests[N]<-cor.test(x=tmf2[,2],y=dp2[dp2[,1] %in% tmf2[,1],N+1],method='spearman')$estimate
tmf3_tests<-c(); for(N in 1:1000) tmf3_tests[N]<-cor.test(x=tmf3[,2],y=dp3[dp3[,1] %in% tmf3[,1],N+1],method='spearman')$estimate

# Plotting the results
#pdf("Proxy_correlation_boxplots.pdf",width=12,height=13)
#par(mfrow=c(4,1))
#boxplot(list("Precipitation"=i_p_tests,"GISP2 d18O"=gp2_tests,"Aridity"=sal_tests,"Sea Surface Temp."=sst_tests,"TMF"=tmf_tests),ylim=c(-1,1)); title(main='18,000 to 8000 cal BP')
#boxplot(list("Precipitation"=i_p1_tests,"GISP2 d18O"=gp21_tests,"Aridity"=NA        ,"Sea Surface Temp."=sst1_tests,"TMF"=tmf1_tests),ylim=c(-1,1)); title(main='Demographic phase 1')
#boxplot(list("Precipitation"=i_p2_tests,"GISP2 d18O"=gp22_tests,"Aridity"=sal2_tests,"Sea Surface Temp."=sst2_tests,"TMF"=tmf2_tests),ylim=c(-1,1)); title(main='Demographic phase 2')
#boxplot(list("Precipitation"=i_p3_tests,"GISP2 d18O"=gp23_tests,"Aridity"=sal3_tests,"Sea Surface Temp."=sst3_tests,"TMF"=tmf3_tests),ylim=c(-1,1)); title(main='Demographic phase 3')
#dev.off()

# Printing the results
paste(
paste(round(mean(gp21_tests),2),round(sd(gp21_tests),2),sep='±'),
paste(round(mean(gp22_tests),2),round(sd(gp22_tests),2),sep='±'),
paste(round(mean(gp23_tests),2),round(sd(gp23_tests),2),sep='±'),sep='  ')

paste(
paste(round(mean(i_p1_tests),2),round(sd(i_p1_tests),2),sep='±'),
paste(round(mean(i_p2_tests),2),round(sd(i_p2_tests),2),sep='±'),
paste(round(mean(i_p3_tests),2),round(sd(i_p3_tests),2),sep='±'),sep='  ')

paste(
paste(round(mean(tmf1_tests),2),round(sd(tmf1_tests),2),sep='±'),
paste(round(mean(tmf2_tests),2),round(sd(tmf2_tests),2),sep='±'),
paste(round(mean(tmf3_tests),2),round(sd(tmf3_tests),2),sep='±'),sep='  ')

paste(
paste(round(mean(sst1_tests),2),round(sd(sst1_tests),2),sep='±'),
paste(round(mean(sst2_tests),2),round(sd(sst2_tests),2),sep='±'),
paste(round(mean(sst3_tests),2),round(sd(sst3_tests),2),sep='±'),sep='  ')

paste(
paste(round(mean(sal2_tests),2),round(sd(sal2_tests),2),sep='±'),
paste(round(mean(sal3_tests),2),round(sd(sal3_tests),2),sep='±'),sep='  ')


# Do the same for the full range of time, i.e. 18000 to 8000 cal BP

tmfo<-tmf[tmf[,1]<18000 & tmf[,1]>=8000,]
salo<-sal[sal[,1]<13000 & sal[,1]>=8000,]
gp2o<-gp2[gp2[,1]<18000 & gp2[,1]>=8000,]
i_po<-i_p[i_p[,1]<18000 & i_p[,1]>=8000,]
ssto<-sst[sst[,1]<18000 & sst[,1]>=8000,]
ssto[,1]<-round(ssto[,1])

i_p_tests<-c(); for(N in 1:616) i_p_tests[N]<-cor.test(x=i_po[,2],y=out[out[,1] %in% i_po[,1],N+1],method='spearman')$estimate
gp2_tests<-c(); for(N in 1:616) gp2_tests[N]<-cor.test(x=gp2o[,2],y=out[out[,1] %in% gp2o[,1],N+1],method='spearman')$estimate
sal_tests<-c(); for(N in 1:616) sal_tests[N]<-cor.test(x=salo[,3],y=out[out[,1] %in% salo[,1],N+1],method='spearman')$estimate
sst_tests<-c(); for(N in 1:616) sst_tests[N]<-cor.test(x=ssto[,2],y=out[out[,1] %in% ssto[,1],N+1],method='spearman')$estimate
tmf_tests<-c(); for(N in 1:616) tmf_tests[N]<-cor.test(x=tmfo[,2],y=out[out[,1] %in% tmfo[,1],N+1],method='spearman')$estimate

paste(round(mean(gp2_tests),2),round(sd(gp2_tests),2),sep='±')
paste(round(mean(sst_tests),2),round(sd(sst_tests),2),sep='±')
paste(round(mean(sal_tests),2),round(sd(sal_tests),2),sep='±')
paste(round(mean(sst_tests),2),round(sd(sst_tests),2),sep='±')
paste(round(mean(tmf_tests),2),round(sd(tmf_tests),2),sep='±')


# Investigating p-value scores
i_p_pv<-c(); for(N in 1:1000) i_p_pv[N]<-cor.test(x=i_po[,2],y=out[out[,1] %in% i_po[,1],N+1],method='spearman')$p.value
gp2_pv<-c(); for(N in 1:1000) gp2_pv[N]<-cor.test(x=gp2o[,2],y=out[out[,1] %in% gp2o[,1],N+1],method='spearman')$p.value
sal_pv<-c(); for(N in 1:1000) sal_pv[N]<-cor.test(x=salo[,3],y=out[out[,1] %in% salo[,1],N+1],method='spearman')$p.value
sst_pv<-c(); for(N in 1:1000) sst_pv[N]<-cor.test(x=ssto[,2],y=out[out[,1] %in% ssto[,1],N+1],method='spearman')$p.value
tmf_pv<-c(); for(N in 1:1000) tmf_pv[N]<-cor.test(x=tmfo[,2],y=out[out[,1] %in% tmfo[,1],N+1],method='spearman')$p.value

i_p1_pv<-c(); for(N in 1:1000) i_p1_pv[N]<-cor.test(x=i_p1[,2],y=dp1[dp1[,1] %in% i_p1[,1],N+1],method='spearman')$p.value
i_p2_pv<-c(); for(N in 1:1000) i_p2_pv[N]<-cor.test(x=i_p2[,2],y=dp2[dp2[,1] %in% i_p2[,1],N+1],method='spearman')$p.value
i_p3_pv<-c(); for(N in 1:1000) i_p3_pv[N]<-cor.test(x=i_p3[,2],y=dp3[dp3[,1] %in% i_p3[,1],N+1],method='spearman')$p.value

gp21_pv<-c(); for(N in 1:1000) gp21_pv[N]<-cor.test(x=gp21[,2],y=dp1[dp1[,1] %in% gp21[,1],N+1],method='spearman')$p.value
gp22_pv<-c(); for(N in 1:1000) gp22_pv[N]<-cor.test(x=gp22[,2],y=dp2[dp2[,1] %in% gp22[,1],N+1],method='spearman')$p.value
gp23_pv<-c(); for(N in 1:1000) gp23_pv[N]<-cor.test(x=gp23[,2],y=dp3[dp3[,1] %in% gp23[,1],N+1],method='spearman')$p.value

sal1_pv<-c(); for(N in 1:1000) sal1_pv[N]<-cor.test(x=sal1[,3],y=dp1[dp1[,1] %in% sal1[,1],N+1],method='spearman')$p.value
sal2_pv<-c(); for(N in 1:1000) sal2_pv[N]<-cor.test(x=sal2[,3],y=dp2[dp2[,1] %in% sal2[,1],N+1],method='spearman')$p.value
sal3_pv<-c(); for(N in 1:1000) sal3_pv[N]<-cor.test(x=sal3[,3],y=dp3[dp3[,1] %in% sal3[,1],N+1],method='spearman')$p.value

sst1_pv<-c(); for(N in 1:1000) sst1_pv[N]<-cor.test(x=sst1[,2],y=dp1[dp1[,1] %in% sst1[,1],N+1],method='spearman')$p.value
sst2_pv<-c(); for(N in 1:1000) sst2_pv[N]<-cor.test(x=sst2[,2],y=dp2[dp2[,1] %in% sst2[,1],N+1],method='spearman')$p.value
sst3_pv<-c(); for(N in 1:1000) sst3_pv[N]<-cor.test(x=sst3[,2],y=dp3[dp3[,1] %in% sst3[,1],N+1],method='spearman')$p.value

tmf1_pv<-c(); for(N in 1:1000) tmf1_pv[N]<-cor.test(x=tmf1[,2],y=dp1[dp1[,1] %in% tmf1[,1],N+1],method='spearman')$p.value
tmf2_pv<-c(); for(N in 1:1000) tmf2_pv[N]<-cor.test(x=tmf2[,2],y=dp2[dp2[,1] %in% tmf2[,1],N+1],method='spearman')$p.value
tmf3_pv<-c(); for(N in 1:1000) tmf3_pv[N]<-cor.test(x=tmf3[,2],y=dp3[dp3[,1] %in% tmf3[,1],N+1],method='spearman')$p.value


paste(
round(median(gp21_pv),2),
round(median(gp22_pv),2),
round(median(gp23_pv),2),sep='  ')

paste(
round(median(i_p1_pv),2),
round(median(i_p2_pv),2),
round(median(i_p3_pv),2),sep='  ')

paste(
round(median(tmf1_pv),2),
round(median(tmf2_pv),2),
round(median(tmf3_pv),2),sep='  ')

paste(
round(median(sst1_pv),2),
round(median(sst2_pv),2),
round(median(sst3_pv),2),sep='  ')

paste(
round(median(sal2_pv),2),
round(median(sal3_pv),2),sep='  ')


signif(median(i_p_pv),2)
signif(median(gp2_pv),2)
signif(median(sal_pv),2)
signif(median(sst_pv),2)
signif(median(tmf_pv),2)


# Plotting the timeseries

pdf('Bootstraps_vs_proxydata.pdf',height=10.5,width=10.5)
par(mar=c(6,4,2,4)+0.1)
plot(out[,1],out[,2],col=NA,xlim=c(18000,7500),yaxt='n',xlab='Calendar Age (BP)',ylim=c(0,7e-04),ylab='')
rug(seq(18500,7000,-250), ticksize = -0.01, side = 1)
for(N in 1:1000) lines(dp1[,1],dp1[,N+1],col=rgb(.8,0,0,.05))
for(N in 1:1000) lines(dp2[,1],dp2[,N+1],col=rgb(0,.5,0,.05))
for(N in 1:1000) lines(dp3[,1],dp3[,N+1],col=rgb(0,0,.5,.05))
text(15619,0.0001,'SPD bootsrap\nsimulations',col=rgb(.8,0,0,.8))
axis(2,at=c(0,0.00005,0.0001))

abline(v=c(16600,break2,break3,8000),col=rgb(.25,.25,.25,.5),lty=3,lwd=1.5)

lines(sst[,1],(sst[,2]/20)*3e-04,xlim=c(18000,8000),type='l')
tp<-c(12,14,16,18,20)
axis(4,at=(tp/20)*3e-04,lab=tp)
text(15619,0.00023,"Sea surface\ntemperature")

lines(gp2[,1],(gp2[,2]+55)/33*6.5e-04,type='l',xlim=c(18000,8000),col='indianred')
tp<-seq(-42,-34,2)
axis(2,at=(tp+55)/33*6.5e-04,lab=tp)
text(15619,0.00035,'GISP2',col='indianred')

lines(i_p[,1],(i_p[,2]/7000+0.00043),col='skyblue')
text(15619,0.00049,'Index of\nprecipitation',col='skyblue')
tp<-c(.2,.4,.6,.8)
axis(4,at=(tp/7000+0.00043),lab=tp)
text(15619,0.00035,'GISP2',col='indianred')

lines(tmf[,1],(tmf[,2]/400000+0.00055),col='darkgreen')
text(15619,0.00061,'TMF',col='darkgreen')
tp<-c(0,20,40,60)
axis(2,at=(tp/400000+0.00055),lab=tp)
text(15619,0.00035,'GISP2',col='indianred')
dev.off()


# Investigating 8.2

pdf('8k2event.pdf',width=5,height=5)

par(mar=c(6,4,2,4)+0.1)
plot(out[,1],out[,2],col=NA,xlim=c(8800,8000),yaxt='n',xlab='Calendar Age (BP)',ylim=c(5e-05,4.5e-04),ylab='')
for(N in 1:1000) lines(dp3[,1],dp3[,N+1],col=rgb(0,0,.5,.05))
axis(2,at=seq(0.00006,0.00022,0.00004))
lines(gp2[,1],(gp2[,2]+44)/15*6e-04,type='l',xlim=c(18000,8000),col='indianred')
tp<-c(-36,-35,-34)
axis(4,at=(tp+44)/15*6e-04,lab=tp)
text(8700,0.000395,'GISP2',col='indianred')
text(8700,0.0002,'SPD bootsrap\nsimulations',col=rgb(0,0,.5,.8))

dev.off()

# Detrending proxy data and calculating variance
library(pracma)

a<-paste("gp2 I",signif(var(detrend(gp21[,2])),2),"gp2 II", signif(var(detrend(gp22[,2])),2),"gp2 III", signif(var(detrend(gp23[,2])),2))
b<-paste("i_p I",signif(var(detrend(i_p1[,2])),2),"i_p II", signif(var(detrend(i_p2[,2])),2),"i_p III", signif(var(detrend(i_p3[,2])),2))
c<-paste("sst I",signif(var(detrend(sst1[,2])),2),"sst II", signif(var(detrend(sst2[,2])),2),"sst III", signif(var(detrend(sst3[,2])),2))
d<-paste("tmf I",signif(var(detrend(tmf1[,2])),2),"tmf II", signif(var(detrend(tmf2[,2])),2),"tmf III", signif(var(detrend(tmf3[,2])),2))
e<-paste("sal II",signif(var(detrend(sal2[,3])),2),"sal III",signif(var(detrend(sal3[,3])),2))

 
