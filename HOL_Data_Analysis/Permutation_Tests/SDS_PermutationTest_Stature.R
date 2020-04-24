# 13th Dec 2019
# R Script for carrying out permutation test on stature QTLs

# First decide which dataset to use
# Set up so one can run from command line
# Otherwise one can set 'ind' manually
# 1 = High N0 case; 2 = Low N0 case
# This version only runs one permutation, so can be parallelised on cluster
args <- commandArgs(trailingOnly = TRUE)
ind <- as.integer(args[1])
idx <- as.integer(args[2])
minbr <- as.integer(args[3])
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
theseed <- sample(2147483647-1,1)
set.seed(theseed)
cat("Seed is ", theseed, "\n")
#setwd("/Users/hartfield/Documents/MilkSDS/HOL_Data_Analysis")
SDSres <- read.table(paste0("../SDSAll_",fname,"N0.dat"),header=T)
cno <- c(1:24,26:29)
orderedChr=paste("Chr",cno, sep="")
SDSres$CHROMOSOME = factor(SDSres$CHROMOSOME,levels=orderedChr)

# Obtaining (1) starting point of each chromosome (2) length of chromosome
# For permutation analyses, to determine where to place SNPs
minp <- vector(mode="numeric",length=length(cno))
maxp <- vector(mode="numeric",length=length(cno))
lengthp <- vector(mode="numeric",length=length(cno))
names(minp) <- factor(orderedChr)
names(maxp) <- factor(orderedChr)
names(lengthp) <- factor(orderedChr)
for(i in unique(SDSres$CHROMOSOME)){
	SDST <- subset(SDSres,CHROMOSOME==i)
	minp[i] <- min(SDST$POS)
	maxp[i] <- max(SDST$POS)
	lengthp[i] <- maxp[i] - minp[i]
}

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

## Finding regions with significantly high absolute SDS.
alpha <- (0.05/dim(SDSres)[1])
SigSDS <- intersect(which((SDSres[,14] < alpha)),which((as.numeric(as.character(SDSres[,9])) >= 1)))
SDSres[SigSDS,] -> SigSDSRes
SDS_CO <- min(SigSDSRes[,13])

# Extracting significant results considering q < 0.05 (FDR)
SigSDS2 <- intersect(which((SDSres[,15] < 0.05)),which((as.numeric(as.character(SDSres[,9])) >= 1)))
SDSres[SigSDS2,] -> SigSDSRes2
SDS_CO2 <- min(SigSDSRes2[,13])

# Removing significant SNPs
cno <- c(1:24,26:29)
orderedChr=paste("Chr",cno, sep="")
SDSres2 <- SDSres[!(row.names(SDSres)%in%row.names(SigSDSRes)),]
row.names(SDSres2) <- c(1:dim(SDSres2)[1])
SDSres2$CHROMOSOME = factor(SDSres2$CHROMOSOME,levels=orderedChr)

# Reading in Stature QTL data. These QTLs have positions relative to old assembly (UMD)
milkQTL <- readRDS(paste0("curated_stature_snps_",minbr,".rds"))
QTLi <- noquote(matrix(data=unlist(strsplit(milkQTL$position,":")),nrow=dim(milkQTL)[1],ncol=2,byrow=T))
QTLi <- data.frame(CHROMOSOME=QTLi[,1],POS=QTLi[,2],milkQTL[,2])
names(QTLi)[3] <- "EFFECT"
QTLi <- QTLi[QTLi$CHROMOSOME!="Chr25",]
QTLi$CHROMOSOME <- factor(QTLi$CHROMOSOME,levels=orderedChr)
QTLi$POS <- as.numeric(QTLi$POS)
QTLi$EFFECT <- as.numeric(QTLi$EFFECT)

# Obtaining new QTL positions in ARS-UCD assembly
nQTL <- read.table("securenew.qtl",head=T)		 		# 'Secure' QTLs
nQTLn <- read.table("nonsecurenew.qtl",head=T)	 		# 'Nonsecure' QTLs
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
nQ <- dim(QTLi)[1];

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
nQ2 <- dim(QTLSDS)[1];	# Number of QTL with SNPs assigned to them

# Setting up GLM, see if stature QTLs significantly differ from background
SDSres2 <- cbind(SDSres2,vector(mode="numeric",length=dim(SDSres2)[1]))
names(SDSres2)[16] <- "isnearQTL"
SDSres2[,16] <- 0
SDSres2[row.names(SDSres2)%in%QTLP,16] <- 1
SDSres2$isnearQTL <- factor(SDSres2$isnearQTL,levels=unique(SDSres2$isnearQTL))
SDSres2$isnearQTL <- relevel(SDSres2$isnearQTL,ref="0")
SDSres2$CHROMOSOME <- factor(SDSres2$CHROMOSOME,levels=unique(SDSres2$CHROMOSOME))
SDSres2$Bin <- factor(SDSres2$Bin,levels=unique(SDSres2$Bin))
SDSres2$Bin <- relevel(SDSres2$Bin,ref="1")

# Running GLM analysis
fitgQ <- glm(sSDS ~ CHROMOSOME + Bin + isnearQTL,data=SDSres2,family=Gamma(link="inverse"))
fitg0 <- glm(sSDS ~ CHROMOSOME + Bin,data=SDSres2,family=Gamma(link="inverse"))
#anova(fitg0,fitgQ,test="Chisq")
pvalA <- anova(fitg0,fitgQ)$"Deviance"[2]

# Print summary to file
if(idx == 1){
	sink(file=paste0("OutTables/Stature_GLM_Summary_",minbr,"_",fname,"N0_FromPermutationTest.txt"))
	cat(paste0("Number of QTL initially ",nQ,"; with SNP effects assigned, ",nQ2,"\n"))
	summary(fitgQ)
	anova(fitg0,fitgQ,test="Chisq")
	sink()
}

# Now permutation analyses

QPermPos <- vector(mode="numeric",length=0)
SDSres3 <- SDSres2
SDSres3[,16] <- 0

# Assigning QTLs to random sites
for(k in 1:nQ2){
	
	if(k%%10==0){
		print(paste0("Assigning QTL number ",k))
	}
	
	isdone <- 0
	while(isdone == 0){
		ctou <- which(rmultinom(1, 1, prob = lengthp)==1)	# Chromosome to place
		qpos <- runif(1,minp[ctou],minp[ctou]+lengthp[ctou])	# Start position
		
		# Defining new QTL position
		SDST <- subset(SDSres3,CHROMOSOME==orderedChr[ctou])
		qpt <- SDST[which(abs(SDST$POS-qpos)==min(abs(SDST$POS-qpos))),5]
		qprt <- row.names(SDST[SDST$POS%in%qpt,])
		
		# Check if region already exists, if not then proceed
		if(length(intersect(qpt,QPermPos)) == 0){
			SDSres3[row.names(SDSres3)%in%qprt,16] <- 1
			QPermPos <- c(QPermPos,qprt)
			isdone <- 1
		}
	}
}

# Performing model fit on permutated dataset
print("Performing model fit on permuted data")
SDSres2$isnearQTL <- factor(SDSres2$isnearQTL,levels=unique(SDSres2$isnearQTL))
SDSres2$isnearQTL <- relevel(SDSres2$isnearQTL,ref="0")
SDSres3$CHROMOSOME <- factor(SDSres3$CHROMOSOME,levels=unique(SDSres3$CHROMOSOME))
SDSres3$Bin <- factor(SDSres3$Bin,levels=unique(SDSres3$Bin))
SDSres2$Bin <- relevel(SDSres3$Bin,ref="1")
fit <- glm(sSDS ~ CHROMOSOME + Bin + isnearQTL,data=SDSres3,family=Gamma(link="inverse"))
pvals <- anova(fitg0,fit)$"Deviance"[2]	# P-value of comparison

# Outputting deviance values
write.table(c(pvalA,pvals),file=paste0("OutTables/PermTest_n",idx,"_",fname,"N0_Stature_",minbr,".dat"), quote = F, row.names = F, col.names = F)

# EOF