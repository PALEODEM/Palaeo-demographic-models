# Convergenge check
# First sort columns sorted to reflect approximate similararity with the median SPD
# 	storing the results in a new matrix called 'test'
# The resulting convergence plot will therefore be of the 'worst case scenario' 

cor<-c()
for(N in 1:1000) cor[N]<-cor.test(T,out[,1+N],method='spearman')$estimate

test<-out[,c(1,1+order(cor))]

# Then compute average N itterations and the overall average and compute difference
# Do this for between 4 itterations and the number of bootstrap itterations contained in the input data

conv<-c(); nc<-ncol(test)-1
T<-apply(test[,2:nc+1],na.rm=T,FUN=median,MARGIN=1)
pb <- txtProgressBar(min=4,max=nc-1,initial=4)
for (N in 4:nc-1) {
 conv[N]<-sum(T,na.rm=T)-sum(apply(test[,2:N],na.rm=T,FUN=median,MARGIN=1))
 setTxtProgressBar(pb,N) }

plot(conv,pch=1,col=rgb(0.7,0.1,0,0.5),xlab="Run number", ylab="Convergence")
abline(h=0,col='grey')


# Alternative measure of convergence

conv2<-c(); nc<-ncol(test)-1
T<-apply(test[,2:nc+1],na.rm=T,FUN=median,MARGIN=1)
conv2<-c(); nc<-ncol(test)-1
pb <- txtProgressBar(min=4,max=nc-1,initial=4)
for (N in 4:nc-1) {
	conv2[N]<-sum(T-colMeans(test[,2:N],na.rm=T))
	setTxtProgressBar(pb,N) }
dev.set(3)
plot(conv2,pch=1,col=rgb(0.7,0.1,0,0.5),xlab="Run number", ylab="Convergence")
abline(h=0,col='grey')


