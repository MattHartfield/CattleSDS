# 18th December 2019
# Code to read in and analyse permutation test results

# Reading in data
setwd("/Users/hartfield/Documents/MilkSDS/HOL_Data_Analysis/")
Milk_HighN0 <- read.table("Permutation_Tests/Results/HOL_Permutation_Milk_HighN0.dat")
Milk_LowN0 <- read.table("Permutation_Tests/Results/HOL_Permutation_Milk_LowN0.dat")
pvalAMH0 <- Milk_HighN0[1,]
pvalAML0 <- Milk_LowN0[1,]
Milk_HighN0P <- Milk_HighN0[-1,]
Milk_LowN0P <- Milk_LowN0[-1,]

# How many permuted deviance values are greater than actual (P-value)?
pvalMH0 <- sum(pvalAMH0 < Milk_HighN0P)/length(Milk_HighN0P)
pvalML0 <- sum(pvalAML0 < Milk_LowN0P)/length(Milk_LowN0P)

# Histogram of permuted values, along with actual value
# Milk, High N0
png(paste0('OutFigures/SDS_Milk_Permutations_HighN0.png'),width=8,height=8,units = 'in',res=200)
par(mar=c(5,9.5,4,2) + 0.1)
maxx <- max(hist(Milk_HighN0P,breaks=30,plot=F)$breaks)
nt <- 5
while(pvalAMH0 > maxx){
	maxx <- maxx*1.2
	nt <- nt + 1
}
maxy <- max(hist(Milk_HighN0P,breaks=30,plot=F)$density)
hist(Milk_HighN0P, breaks=30, prob=T, col=rgb(0,0,0,0.25), xlab="", ylab="", main="", xaxt="n", yaxt="n", xlim=c(0,maxx))
axis(1, at=seq(0, maxx, maxx/nt), pos=0)
axis(2, at=seq(0,maxy,maxy/4), las=2)
title(bquote(bold("Randomised Values, Milk Protein Genes (High"~N[0]*")")),xlab="Deviance Values")
text(x=(-maxx/3.25),y=maxy/2,labels=paste("Density",sep="\n"),xpd=NA)
abline(v=pvalAMH0,lty=2,lwd=2)
text(x=(3*maxx/4),y=maxy/2,labels=bquote(italic('P')*"\u2013value ="~.(pvalMH0)),xpd=NA)
dev.off()

# Milk, Low N0
png(paste0('OutFigures/SDS_Milk_Permutations_LowN0.png'),width=8,height=8,units = 'in',res=200)
par(mar=c(5,9.5,4,2) + 0.1)
maxx <- max(hist(Milk_LowN0P,breaks=30,plot=F)$breaks)
nt <- 5
while(pvalAML0 > maxx){
	maxx <- maxx*1.2
	nt <- nt + 1
}
maxy <- max(hist(Milk_LowN0P,breaks=30,plot=F)$density)
hist(Milk_LowN0P, breaks=30, prob=T, col=rgb(0,0,0,0.25), xlab="", ylab="", main="", xaxt="n", yaxt="n", xlim=c(0,maxx))
axis(1, at=seq(0, maxx, maxx/nt), pos=0)
axis(2, at=seq(0,maxy,maxy/4), las=2)
title(bquote(bold("Randomised Values, Milk Protein Genes (Low"~N[0]*")")),xlab="Deviance Values")
text(x=(-maxx/3.25),y=maxy/2,labels=paste("Density",sep="\n"),xpd=NA)
abline(v=pvalAML0,lty=2,lwd=2)
text(x=(3*maxx/4),y=maxy/2,labels=bquote(italic('P')*"\u2013value ="~.(pvalML0)),xpd=NA)
dev.off()

# EOF