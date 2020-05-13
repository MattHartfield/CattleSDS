# 22nd May 2019
# R Script for analysing SDS results, both high and low N0 results

# First decide which dataset to use
# Set up so one can run from command line
# Otherwise one can set 'ind' manually
# 1 = High N0 case; 2 = Low N0 case
args <- commandArgs(trailingOnly = TRUE)
ind <- as.integer(args[1])
if(ind == 1){
	fname <- "High"
}else if(ind == 2){
	fname <- "Low"
}

# Half-Normal CDF (for calculating P-values)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
fncdf <- function(x){
	ifelse(x < 0, 0,  erf(x/sqrt(2)))
}
fnP <- function(x){
	1-fncdf(x)
}

library(qvalue)
setwd("/Users/hartfield/Documents/MilkSDS/HOL_Data_Analysis")
SDSres <- read.table(paste0("SDSAll_",fname,"N0.dat"),header=T)
cno <- c(1:24,26:29)
orderedChr=paste("Chr",cno, sep="")
SDSres$CHROMOSOME = factor(SDSres$CHROMOSOME,levels=orderedChr)

SDSres <- cbind(SDSres,as.numeric(cut(SDSres[,6],seq(0.05,0.95,0.05),include.lowest = T)))
row.names(SDSres) <- c(1:dim(SDSres)[1])
names(SDSres)[12] <- "Bin"

binM <- tapply(as.numeric(as.character(SDSres[,10])), SDSres[,12],function(x) (mean(x,na.rm=T)))
binSD <- tapply(as.numeric(as.character(SDSres[,10])), SDSres[,12],function(x) (sd(x,na.rm=T)))

## Normalising SDS scores per bin of DAF (taking absolute values)
SDSres <- cbind(SDSres,vector(mode="numeric",length=dim(SDSres)[1]),vector(mode="numeric",length=dim(SDSres)[1]))
names(SDSres)[13] <- "sSDS"
names(SDSres)[14] <- "pSDS"
for(i in 1:length(binM)){
	SDSres[SDSres[,12]==i,13] <- abs((SDSres[SDSres[,12]==i,10]-binM[i])/binSD[i])
	SDSres[SDSres[,12]==i,14] <- fnP(SDSres[SDSres[,12]==i,13])
}
# Calculating FDR values from p-values
qobj <- qvalue(p = SDSres[,14])
SDSres[,15] <- qobj$qvalues
names(SDSres)[15] <- "qSDS"

## Part 1: Finding regions with significantly high absolute SDS.
alpha <- (0.05/dim(SDSres)[1])
SigSDS <- intersect(which((SDSres[,14] < alpha)),which((as.numeric(as.character(SDSres[,9])) >= 1)))
SDSres[SigSDS,] -> SigSDSRes
SDS_CO <- min(SigSDSRes[,13])
cat(paste0(dim(SigSDSRes)[1]," SNPs with Bonferroni-corrected P < 0.05\n"))

# Extracting significant results considering q < 0.05 (FDR)
SigSDS2 <- intersect(which((SDSres[,15] < 0.05)),which((as.numeric(as.character(SDSres[,9])) >= 1)))
SDSres[SigSDS2,] -> SigSDSRes2
SDS_CO2 <- min(SigSDSRes2[,13])
cat(paste0(dim(SigSDSRes2)[1]," SNPs with FDR < 0.05\n"))

# Plotting absolute SDS, with outlier regions
colrs <- vector(mode="character",length=dim(SDSres)[1])
tickpos <- vector(mode="numeric",length=length(cno))

for(i in 1:24){
	if(which(cno==i)%%2 == 1){
		colrs[which(SDSres[,1]==paste0("Chr",i))] <- "gray30"
	}else if(which(cno==i)%%2 == 0){
		colrs[which(SDSres[,1]==paste0("Chr",i))] <- "gray60"
	}
	tickpos[which(cno==i)] <- floor(median(which(SDSres[,1]==paste0("Chr",i))))
}

for(i in 26:29){
	if(which(cno==i)%%2 == 0){
		colrs[which(SDSres[,1]==paste0("Chr",i))] <- "gray30"
	}else if(which(cno==i)%%2 == 1){
		colrs[which(SDSres[,1]==paste0("Chr",i))] <- "gray60"
	}
	tickpos[which(cno==i)] <- floor(median(which(SDSres[,1]==paste0("Chr",i))))
}

colrs[intersect(SigSDS,SigSDS2)] <- "red"
colrs[setdiff(SigSDS,SigSDS2)] <- "purple"
colrs[setdiff(SigSDS2,SigSDS)] <- "blue"

png(paste0('OutFigures/SDS_Chr_All_AbsVal_',fname,'N0.png'),width=4*(1+sqrt(5)),height=8,units = 'in',res=200)
par(mar=c(5,8,4,2) + 0.1)
plot(SDSres[,13],col=colrs,pch=16, xaxt="n", yaxt="n", xlab="Chromosome",ylab="",main=substitute(bold(paste("Absolute Standardised SDS for ",bolditalic('Bos taurus')," Autosomes"))),ylim=c(0,8))
axis(1, at=tickpos, labels=cno)
axis(2, at=seq(0,8,2), las=2)
text(x=-275000,y=4,labels=paste("Absolute","Standardised","SDS",sep="\n"),xpd=NA)
abline(SDS_CO,0,lty=2)
abline(SDS_CO2,0,lty=3)
legend("topleft", legend=c(expression("Bonferroni\u2013"*"Corrected "*italic("P")*" < 0.05"), "FDR < 0.05"),col=c("red", "blue"),pch=16)
dev.off()

## Part 2: Are any high SDS scores within milk-producing regions?
# Reading in and classifying milk-producing regions
milky <- read.table("Milk_Protein_Genes/milk_protein_names.bed",skip=1)[,1:3] # Only including gene locations and autosomal genes
names(milky) <- c("chrom","chromStart","chromEnd")
milky$chrom <- factor(milky$chrom, levels=unique(SDSres$CHROMOSOME))
sigmilk <- vector(mode="numeric",length=0)
for(i in unique(milky$chrom)){
	hpss <- as.matrix(subset(milky,chrom==i)[2:3])
	SDST <- subset(SigSDSRes,CHROMOSOME==i)
	SDST2 <- subset(SigSDSRes2,CHROMOSOME==i)	
	sigmilk <-c(sigmilk,unlist(apply(hpss,1,function(x) row.names(SDST[intersect(which(SDST[SDST$CHROMOSOME==i,5] >= x[1]), which(SDST[SDST$CHROMOSOME==i,5] <= x[2])),]))))
}
sigmilk		# Milk genes overlapping with highly-significant SDS regions
SigSDSRes2[row.names(SigSDSRes2)%in%sigmilk,]

# Print BED file of significant regions for bedtools analysis
space <- 10000		# How much space to add either side (to look for nearby genes as well)
sigbed <- cbind(as.numeric(unlist(strsplit(as.character(SigSDSRes$CHROMOSOME),"r"))[seq(2,2*dim(SigSDSRes)[1],2)]),SigSDSRes$POS-space,SigSDSRes$POS+space)
# Merging overlapping windows
for(i in cno){
	sbt <- sigbed[sigbed[,1]==i,]
	if(is.null(dim(sbt)) != T){
		if(dim(sbt)[1] > 0){
			sp <- 1
			ep <- 2
			done <- 0
			while(done == 0){
				if(sbt[sp,3] > sbt[ep,2]){
					sbt[sp,3] <- sbt[ep,3]
					sbt <- sbt[-ep,]
				}else if(sbt[sp,3] <= sbt[ep,2]){
					sp <- sp + 1
					ep <- ep + 1
				}
				if(is.null(dim(sbt))==T || sp == dim(sbt)[1]){
					sigbed <- sigbed[sigbed[,1]!=i,]
					sigbed <- rbind(sigbed,sbt)
					done <- 1
				}
			}
		}
	}else if(is.null(dim(sbt)) == T){
		sigbed <- sigbed[sigbed[,1]!=i,]
		sigbed <- rbind(sigbed,sbt)
	}
}
row.names(sigbed) <- NULL
cat("track name=sigSDSregions\n", file=paste0("Sigsds_",fname,"N0.bed"))
write.table(sigbed,file=paste0("Sigsds_",fname,"N0.bed"),col.names=F,row.names=F,quote=F,append=F,sep="\t")

# Removing significant SNPs
cno <- c(1:24,26:29)
orderedChr=paste("Chr",cno, sep="")
SDSres2 <- SDSres[!(row.names(SDSres)%in%row.names(SigSDSRes)),]
row.names(SDSres2) <- c(1:dim(SDSres2)[1])
SDSres2$CHROMOSOME = factor(SDSres2$CHROMOSOME,levels=orderedChr)

## Part 3: including milk-gene as a factor. Putting scores in bins, determining effect of milk on SDS.

# Classifying milk-producing regions
milkreg <- vector(mode="numeric",length=0)
for(i in unique(milky$chrom)){
	hpss <- as.matrix(subset(milky,chrom==i)[2:3])
	SDST <- subset(SDSres2,CHROMOSOME==i)
	milkreg<-c(milkreg,unlist(apply(hpss,1,function(x) row.names(SDST[intersect(which(SDST[SDST$CHROMOSOME==i,5] >= x[1]), which(SDST[SDST$CHROMOSOME==i,5] <= x[2])),]))))
}
names(milkreg) <- NULL

# Adding column stating if SNPs lie within a milk-producing region
SDSres2 <- cbind(SDSres2,vector(mode="numeric",length=dim(SDSres2)[1]))
names(SDSres2)[16] <- "Milk"
SDSres2[,16] <- 0
SDSres2[row.names(SDSres2)%in%milkreg,16] <- 1
SDSres2[!(row.names(SDSres2)%in%milkreg),16] <- 0
SDSres2$Milk <- factor(SDSres2$Milk,levels=unique(SDSres2$Milk))

# Plotting histogram of sSDS (red) with milk histogram on top (blue)
# Acknowledgement for overlapping histograms: https://www.r-bloggers.com/overlapping-histogram-in-r/
MilkSDS <- subset(SDSres2,Milk==1)$sSDS
maxx <- ceiling(max(hist(SDSres2$sSDS, breaks=30, plot=F)$breaks))
maxy <- max(max(hist(MilkSDS, breaks=30, plot=F)$density),max(hist(SDSres2$sSDS, breaks=30, plot=F)$density))
png(paste0('OutFigures/SDS_Abs_HistS_',fname,'N0.png'),width=8,height=8,units = 'in',res=200)
par(mar=c(5,6.5,4,1.5) + 0.1)
hist(SDSres2$sSDS, breaks=30, prob=T, col=rgb(1,0,0,0.5), xlim=c(0,maxx), ylim=c(0, maxy), xaxt="n", yaxt="n", xlab="",ylab="",main="")
axis(1,at=seq(0, maxx,1),pos=(0))
axis(2, at=seq(0,0.8,0.2), las=2)
title("Absolute Standardised SDS with Milk Protein Genes Values",xlab="Absolute Standardised SDS")
text(x=-1.25,y=0.4,labels=paste("Density",sep="\n"),xpd=NA)
hist(MilkSDS, breaks=30, prob=T, col=rgb(0,0,1,0.5), add=T)
legend("topright", legend=c("All regions", "Milk Protein Genes"),col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)),lwd=10)
dev.off()

# General Linear model fit, comparing milk and non-milk regions, Gamma link function.
SDSres2$CHROMOSOME <- factor(SDSres2$CHROMOSOME,levels=unique(SDSres2$CHROMOSOME))
SDSres2$Bin <- factor(SDSres2$Bin,levels=unique(SDSres2$Bin))
SDSres2$Bin <- relevel(SDSres2$Bin,ref="1")
SDSres2$Milk <- relevel(SDSres2$Milk,ref="0")
fitg <- glm(sSDS ~ CHROMOSOME + Bin + Milk,data=SDSres2,family=Gamma(link="inverse"));summary(fitg)
fitg0 <- glm(sSDS ~ CHROMOSOME + Bin,data=SDSres2,family=Gamma(link="inverse"));summary(fitg0)
anova(fitg0,fitg,test="Chisq")

# Print summary to file
sink(file=paste0("OutTables/Milk_GLM_Summary_",fname,"N0.txt"))
summary(fitg)
anova(fitg0,fitg,test="Chisq")
sink()

# Part 4: Reading in Stature QTL data. First where effect sizes estimated in 6 of 7 Holstein
# These QTLs have positions relative to old assembly (UMD)
milkQTL <- readRDS("QTLStatureDat/curated_stature_snps_6.rds")
QTLi <- noquote(matrix(data=unlist(strsplit(milkQTL$position,":")),nrow=dim(milkQTL)[1],ncol=2,byrow=T))
QTLi <- data.frame(CHROMOSOME=QTLi[,1],POS=QTLi[,2],milkQTL[,2])
names(QTLi)[3] <- "EFFECT"
QTLi <- QTLi[QTLi$CHROMOSOME!="Chr25",]
QTLi$CHROMOSOME <- factor(QTLi$CHROMOSOME,levels=orderedChr)
QTLi$POS <- as.numeric(QTLi$POS)
QTLi$EFFECT <- as.numeric(QTLi$EFFECT)

# Obtaining new QTL positions in ARS-UCD assembly
nQTL <- read.table("StatureQTLs/securenew.qtl",head=T)		 		# 'Secure' QTLs
nQTLn <- read.table("StatureQTLs/nonsecurenew.qtl",head=T)	 		# 'Nonsecure' QTLs
nQTL <- rbind(nQTL,nQTLn[,c(1:3,5)])
nQTL <- nQTL[order(nQTL$UMDChr, nQTL$UMDPos),]
nQTL$UMDChr <- paste("Chr",nQTL$UMDChr,sep="")
nQTL$ARSChr <- paste("Chr",nQTL$ARSChr,sep="")
nQTL <- nQTL[nQTL$UMDChr!="Chr25",]
nQTL$UMDChr <- factor(nQTL$UMDChr,levels=orderedChr)
nQTL$ARSChr <- factor(nQTL$ARSChr,levels=orderedChr)
row.names(nQTL) <- c(1:dim(nQTL)[1])

# Now matching old positions to new ones
nidx <- match(QTLi$POS,nQTL$UMDPos)								# Finding new positions in lookup table
QTLi <- cbind(QTLi,nQTL[nidx,c(3,4)])
QTLi <- QTLi[!is.na(QTLi$ARSPos),]
QTLi <- QTLi[QTLi$CHROMOSOME==QTLi$ARSChr,]
QTLi$POS <- QTLi$ARSPos
QTLi <- QTLi[,1:3]

# Initial number of QTL (excluding Chr 25)
nQ <- dim(QTLi)[1];nQ

# Finding nearest SNP to each QTL
QTLSDS <- matrix(data=NA,nrow=0,ncol=2)
QTLP <- vector(mode="numeric",length=0)
for(i in unique(SDSres2$CHROMOSOME)){
	SDST <- subset(SDSres2,CHROMOSOME==i)
	bounds <- c(min(SDST $POS),max(SDST$POS))
	QTLT <- subset(QTLi,CHROMOSOME==i)
	QTLT <- QTLT[((QTLT$POS>=bounds[1])+(QTLT$POS<=bounds[2]) == 2),] 	# Removing QTLs that lie in telomeric regions
	QTLT <- QTLT[!(apply(QTLT,1,function(x) SDST[which(abs(SDST$POS-as.numeric(x[2]))==min(abs(SDST$POS-as.numeric(x[2])))),5])%in%bounds),]	# Removing QTLs where nearest SNP is at edge of range
	PT <- apply(QTLT,1,function(x) SDST[which(abs(SDST$POS-as.numeric(x[2]))==min(abs(SDST$POS-as.numeric(x[2])))),5])
	QTLP <- c(QTLP,row.names(SDST[SDST$POS%in%PT,]))
	QTLSDS <- rbind(QTLSDS,cbind(SDST[SDST$POS%in%PT,13], abs(QTLT[,3])))
}
nQ2 <- dim(QTLSDS)[1];nQ2	# Number of QTL with SNPs assigned to them

# Setting up GLM, see if stature QTLs significantly differ from background
SDSres2 <- cbind(SDSres2,vector(mode="numeric",length=dim(SDSres2)[1]))
names(SDSres2)[17] <- "isnearQTL"
SDSres2[,17] <- 0
SDSres2[row.names(SDSres2)%in%QTLP,17] <- 1
SDSres2$isnearQTL <- factor(SDSres2$isnearQTL,levels=unique(SDSres2$isnearQTL))
SDSres2$isnearQTL <- relevel(SDSres2$isnearQTL,ref="0")

# Plotting histogram of sSDS (red) with stature QTL histogram on top (blue)
maxx <- ceiling(max(hist(SDSres2$sSDS, breaks=30, plot=F)$breaks))
maxy <- max(max(hist(QTLSDS[,1], breaks=30, plot=F)$density),max(hist(SDSres2$sSDS, breaks=30, plot=F)$density))
png(paste0('OutFigures/SDS_Abs_HistS_',fname,'N0_StatureQTL_6.png'),width=8,height=8,units = 'in',res=200)
par(mar=c(5,6.5,4,1.5) + 0.1)
hist(SDSres2$sSDS, breaks=30, prob=T, col=rgb(1,0,0,0.5), xlim=c(0,maxx), ylim=c(0, maxy+0.2), xaxt="n", yaxt="n", xlab="",ylab="",main="")
axis(1,at=seq(0, maxx,1),pos=(0))
axis(2, at=seq(0,maxy+0.2,0.2), las=2)
title("Absolute Standardised SDS with Stature QTL Values",xlab="Absolute Standardised SDS")
text(x=-1.25,y=(maxy/2),labels=paste("Density",sep="\n"),xpd=NA)
hist(QTLSDS[,1], breaks=30, prob=T, col=rgb(0,0,1,0.5), add=T)
legend("topright", legend=c("All regions", "SNPs near stature QTLs"),col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)),lwd=10)
dev.off()

# Running GLM analysis
fitgQ <- glm(sSDS ~ CHROMOSOME + Bin + isnearQTL,data=SDSres2,family=Gamma(link="inverse"));summary(fitgQ)
anova(fitg0,fitgQ,test="Chisq")

# Print summary to file
sink(file=paste0("OutTables/Stature_GLM_Summary_6_",fname,"N0.txt"))
cat(paste0("Number of QTL initially ",nQ,"; with SNP effects assigned, ",nQ2,"\n"))
summary(fitgQ)
anova(fitg0,fitgQ,test="Chisq")
sink()

# Part 4a: Repeating for QTLs with effect sizes in at least 5 of 7 Holstein
# Original co-ordinates in UMD assembly
milkQTL <- readRDS("QTLStatureDat/curated_stature_snps_5.rds")
QTLi <- noquote(matrix(data=unlist(strsplit(milkQTL$position,":")),nrow=dim(milkQTL)[1],ncol=2,byrow=T))
QTLi <- data.frame(CHROMOSOME=QTLi[,1],POS=QTLi[,2],milkQTL[,2])
names(QTLi)[3] <- "EFFECT"
QTLi <- QTLi[QTLi$CHROMOSOME!="Chr25",]
QTLi$CHROMOSOME <- factor(QTLi$CHROMOSOME,levels=orderedChr)
QTLi$POS <- as.numeric(QTLi$POS)
QTLi$EFFECT <- as.numeric(QTLi$EFFECT)

# Now matching old positions to new ones
nidx <- match(QTLi$POS,nQTL$UMDPos)								# Finding new positions in lookup table
QTLi <- cbind(QTLi,nQTL[nidx,c(3,4)])
QTLi <- QTLi[!is.na(QTLi$ARSPos),]
QTLi <- QTLi[QTLi$CHROMOSOME==QTLi$ARSChr,]
QTLi$POS <- QTLi$ARSPos
QTLi <- QTLi[,1:3]

nQ <- dim(QTLi)[1];nQ	# Initial number of QTL (excluding Chr 25)

# Finding nearest SNP to each QTL
QTLSDS <- matrix(data=NA,nrow=0,ncol=2)
QTLP <- vector(mode="numeric",length=0)
for(i in unique(SDSres2$CHROMOSOME)){
	SDST <- subset(SDSres2,CHROMOSOME==i)
	bounds <- c(min(SDST $POS),max(SDST$POS))
	QTLT <- subset(QTLi,CHROMOSOME==i)
	QTLT <- QTLT[((QTLT$POS>=bounds[1])+(QTLT$POS<=bounds[2]) == 2),] 	# Removing QTLs that lie in telomeric regions
	QTLT <- QTLT[!(apply(QTLT,1,function(x) SDST[which(abs(SDST$POS-as.numeric(x[2]))==min(abs(SDST$POS-as.numeric(x[2])))),5])%in%bounds),]	# Removing QTLs where nearest SNP is at edge of range
	PT <- apply(QTLT,1,function(x) SDST[which(abs(SDST$POS-as.numeric(x[2]))==min(abs(SDST$POS-as.numeric(x[2])))),5])
	QTLP <- c(QTLP,row.names(SDST[SDST$POS%in%PT,]))
	QTLSDS <- rbind(QTLSDS,cbind(SDST[SDST$POS%in%PT,13], abs(QTLT[,3])))
}
nQ2 <- dim(QTLSDS)[1];nQ2	# Number of QTL with SNPs assigned to them

# Setting up GLM, see if stature QTLs significantly differ from background
SDSres2[,17] <- 0
SDSres2[row.names(SDSres2)%in%QTLP,17] <- 1
SDSres2$isnearQTL <- factor(SDSres2$isnearQTL,levels=unique(SDSres2$isnearQTL))
SDSres2$isnearQTL <- relevel(SDSres2$isnearQTL,ref="0")

# Plotting histogram of sSDS (red) with stature QTL histogram on top (blue)
maxx <- ceiling(max(hist(SDSres2$sSDS, breaks=30, plot=F)$breaks))
maxy <- max(max(hist(QTLSDS[,1], breaks=30, plot=F)$density),max(hist(SDSres2$sSDS, breaks=30, plot=F)$density))
png(paste0('OutFigures/SDS_Abs_HistS_',fname,'N0_StatureQTL_5.png'),width=8,height=8,units = 'in',res=200)
par(mar=c(5,6.5,4,1.5) + 0.1)
hist(SDSres2$sSDS, breaks=30, prob=T, col=rgb(1,0,0,0.5), xlim=c(0,maxx), ylim=c(0, maxy+0.2), xaxt="n", yaxt="n", xlab="",ylab="",main="")
axis(1,at=seq(0, maxx,1),pos=(0))
axis(2, at=seq(0,maxy+0.2,0.2), las=2)
title("Absolute Standardised SDS with Stature QTL Values",xlab="Absolute Standardised SDS")
text(x=-1.25,y=(maxy/2),labels=paste("Density",sep="\n"),xpd=NA)
hist(QTLSDS[,1], breaks=30, prob=T, col=rgb(0,0,1,0.5), add=T)
legend("topright", legend=c("All regions", "SNPs near stature QTLs"),col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)),lwd=10)
dev.off()

# Running GLM analysis
fitgQ <- glm(sSDS ~ CHROMOSOME + Bin + isnearQTL,data=SDSres2,family=Gamma(link="inverse"));summary(fitgQ)
anova(fitg0,fitgQ,test="Chisq")

# Print summary to file
sink(file=paste0("OutTables/Stature_GLM_Summary_5_",fname,"N0.txt"))
cat(paste0("Number of QTL initially ",nQ,"; with SNP effects assigned, ",nQ2,"\n"))
summary(fitgQ)
anova(fitg0,fitgQ,test="Chisq")
sink()

# EOF