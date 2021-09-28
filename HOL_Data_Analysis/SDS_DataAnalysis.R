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
}else{
	cat("Input command needs to be 1 or 2. Exiting.\n")
	quit(status=1)
}
cat("Analysing",fname,"N0 data\n")

# Function for finding nearest position to QTL, polarising. Milk QTLs
QTLpol <- function(milkfQTL,SDSres2)
{
	# Replacing position, chromosome with ARS-UCD versions
	milkfQTL$pos <- unlist(strsplit(milkfQTL$`ARS-UCD2.1`,":"))[seq(2,dim(milkfQTL)[1]*2,2)]
	milkfQTL$chr <- chartr("c","C",unlist(strsplit(milkfQTL$`ARS-UCD2.1`,":"))[seq(1,(dim(milkfQTL)[1]*2-1),2)])
	milkfQTL$chr <- factor(milkfQTL$chr,levels=orderedChr)
	milkfQTL$pos <- as.numeric(milkfQTL$pos)

	# Reading in table of allele polarisation
	apol <- read.table(paste0("SDSAll_Polarisation.dat"),header=T)
	apol$CHROMOSOME = factor(apol$CHROMOSOME,levels=orderedChr)

	# Finding nearest SNP to each QTL, obtaining SDS score
	QTLmfs <- data.frame(CHROMOSOME=character(),POS=numeric(),sSDS=numeric(),Bin=numeric(),NotRef=numeric(),p=numeric(),German_Holstein=character(),NegEff=numeric())
	QTLmfpos <- vector(mode="numeric",length=0) # Distance to nearest SNP
	for(i in unique(SDSres2$CHROMOSOME)){
		SDST <- subset(SDSres2,CHROMOSOME==i)
		bounds <- c(min(SDST$POS),max(SDST$POS))
		QTLT <- subset(milkfQTL,chr==i)
		QTLT <- QTLT[((QTLT$pos>=bounds[1])+(QTLT$pos<=bounds[2]) == 2),] 	# Removing QTLs that lie in telomeric regions
		QTLT <- QTLT[!(apply(QTLT,1,function(x) SDST[which(abs(SDST$POS-as.numeric(x[3]))==min(abs(SDST$POS-as.numeric(x[3])))),5])%in%bounds),]		# Removing QTLs where nearest SNP is at edge of range
		PT <- apply(QTLT,1,function(x) SDST[which(abs(SDST$POS-as.numeric(x[3]))==min(abs(SDST$POS-as.numeric(x[3])))),5])
		QTLmfpos <- c(QTLmfpos,abs(PT-QTLT$pos))
		QTLmfsT <- SDST[SDST$POS%in%PT,c(1,5,13,12)]
		if(dim(QTLmfsT)[1]>0)
		{
			# Finding QTLs where REF!=ANC, switching signs
			QTLmfsT <- cbind(QTLmfsT,0)
			names(QTLmfsT)[5] <- "NotRef"
			apT <- filter(apol,CHROMOSOME==i,POS%in%PT,ANC!=REF)
			QTLmfsT[QTLmfsT$POS%in%apT$POS,5] <- 1
			QTLmfsT[QTLmfsT$POS%in%apT$POS,3] <- ifelse(QTLmfsT[QTLmfsT$POS%in%apT$POS,3]<0,abs(QTLmfsT[QTLmfsT$POS%in%apT$POS,3]),(-1)*(QTLmfsT[QTLmfsT$POS%in%apT$POS,3]))
			# Then polarising by QTL effect
			QTLmfsT <- cbind(QTLmfsT,abs(log10(QTLT[,5])),QTLT[,7],0)
			names(QTLmfsT)[8] <- "RePol"
			for(j in 1:dim(QTLmfsT)[1])
			{
				# Polarising SDS value with effect direction
				if(QTLmfsT[j,7]=="+" & QTLmfsT[j,3]<0)
				{
					QTLmfsT[j,3] <- abs(QTLmfsT[j,3])
					QTLmfsT[j,8] <- 1
				}else if(QTLmfsT[j,7]=="-" & QTLmfsT[j,3]>0)
				{
					QTLmfsT[j,3] <- (-1)*(QTLmfsT[j,3])
					QTLmfsT[j,8] <- 1
				}
			}
			QTLmfs <- rbind(QTLmfs,QTLmfsT)
		}
	}
	
	outv <- list("Distances"=QTLmfpos,"Values"=QTLmfs)
	return(outv)
	
}

# Function for finding nearest position to QTL, polarising. Stature QTLs
QTLspol <- function(statQTL6,SDSres2)
{
	QTLi <- noquote(matrix(data=unlist(strsplit(statQTL6$position,":")),nrow=dim(statQTL6)[1],ncol=2,byrow=T))
	QTLi <- data.frame(CHROMOSOME=QTLi[,1],POS=QTLi[,2],statQTL6[,c(2,3,8)])
	names(QTLi)[3] <- "EFFECT"
	names(QTLi)[4] <- "p"
	names(QTLi)[5] <- "QTLAncestral"
	QTLi <- QTLi[QTLi$CHROMOSOME!="Chr25",]
	QTLi$CHROMOSOME <- factor(QTLi$CHROMOSOME,levels=orderedChr)
	QTLi$POS <- as.numeric(QTLi$POS)
	QTLi$EFFECT <- as.numeric(QTLi$EFFECT)
	QTLi$p <- as.numeric(QTLi$p)
	QTLi$QTLAncestral <- as.character(QTLi$QTLAncestral)

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
	QTLi <- QTLi[,1:5]
	
	# Initial number of QTL (excluding Chr 25)
	nQ <- dim(QTLi)[1];nQ
	
	# Reading in table of allele polarisation
	apol <- read.table(paste0("SDSAll_Polarisation.dat"),header=T)
	apol$CHROMOSOME = factor(apol$CHROMOSOME,levels=orderedChr)
	
	# Finding nearest SNP to each QTL
	QTLsts <- data.frame(CHROMOSOME=character(),POS=numeric(),sSDS=numeric(),Bin=numeric(),QTLAncestral=character(),ANC=character(),NotRef=numeric(),EFFECT=numeric(),p=numeric(),NegEff=numeric())
	QTLsts2 <- data.frame(CHROMOSOME=character(),POS=numeric(),sSDS=numeric(),Bin=numeric(),QTLAncestral=character(),ANC=character(),NotRef=numeric(),EFFECT=numeric(),p=numeric(),NegEff=numeric())
	QTLspos <- vector(mode="numeric",length=0) # Distance to nearest SNP
	for(i in unique(SDSres2$CHROMOSOME)){
		SDST <- subset(SDSres2,CHROMOSOME==i)
		bounds <- c(min(SDST $POS),max(SDST$POS))
		QTLT <- subset(QTLi,CHROMOSOME==i)
		QTLT <- QTLT[((QTLT$POS>=bounds[1])+(QTLT$POS<=bounds[2]) == 2),] 	# Removing QTLs that lie in telomeric regions
		QTLT <- QTLT[!(apply(QTLT,1,function(x) SDST[which(abs(SDST$POS-as.numeric(x[2]))==min(abs(SDST$POS-as.numeric(x[2])))),5])%in%bounds),]	# Removing QTLs where nearest SNP is at edge of range
		PT <- apply(QTLT,1,function(x) SDST[which(abs(SDST$POS-as.numeric(x[2]))==min(abs(SDST$POS-as.numeric(x[2])))),5])
		if(length(PT)>0)
		{
			QTLspos <- c(QTLspos,abs(PT-QTLT$POS))
			# Now forming table of nearest SDS scores, along with ancestral allele calls
			PT2 <- apply(QTLT,1,function(x) c(SDST[which(abs(SDST$POS-as.numeric(x[2]))==min(abs(SDST$POS-as.numeric(x[2])))),5],x[5]))
			PT2 <- as.data.frame(t(PT2))
			names(PT2)[1] <- "POS"
			PT2$POS <- as.numeric(as.character(PT2$POS))
			PT3 <- SDST[SDST$POS%in%PT,c(1,5,13,12)]
			PT4 <- filter(apol,CHROMOSOME==i,POS%in%PT)[,2:3]
			PT5 <- inner_join(PT3,PT2,by="POS") %>% inner_join(PT4,by="POS")

			# Finding QTLs where REF!=ANC, switching signs
			apT <- filter(apol,CHROMOSOME==i,POS%in%PT,ANC!=REF)
			PT5 <- cbind(PT5,0)
			names(PT5)[7] <- "NotRef"
			PT5[PT5$POS%in%apT$POS,7] <- 1
			PT5[PT5$POS%in%apT$POS,3] <- ifelse(PT5[PT5$POS%in%apT$POS,3]<0,abs(PT5[PT5$POS%in%apT$POS,3]),(-1)*(PT5[PT5$POS%in%apT$POS,3]))
			# Then polarising by QTL effect
			PT5 <- cbind(PT5,QTLT[,3],abs(log10(QTLT[,4])),0)
			names(PT5)[8] <- "EFFECT"
			names(PT5)[9] <- "p"
			names(PT5)[10] <- "RePol"
			for(j in 1:dim(PT5)[1])
			{
				# Aligning effect size and SDS score
				if(PT5[j,8] < 0 & PT5[j,3] > 0)
				{
					PT5[j,3] <- (-1)*(PT5[j,3])
					PT5[j,10] <- 1
				}else if(PT5[j,8] > 0 & PT5[j,3] < 0)
				{
					PT5[j,3] <- abs(PT5[j,3])
					PT5[j,10] <- 1
				}
			}
			QTLsts <- rbind(QTLsts,PT5)
			QTLsts2 <- rbind(QTLsts2,filter(PT5,POS%in%PT[PT==QTLT$POS]))			
		}	
	}
	QTLsts2 <- filter(QTLsts2,QTLAncestral!="?")
	QTLsts2$QTLAncestral <- factor(QTLsts2$QTLAncestral,unique(QTLsts2$QTLAncestral))
	outv <- list("Distances"=QTLspos,"Values"=QTLsts,"Exact"=QTLsts2)
	return(outv)
}

# Loading libraries and reading in data
library(qvalue)
library(tidyverse)
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

## Normalising SDS scores per bin of DAF
SDSres <- cbind(SDSres,vector(mode="numeric",length=dim(SDSres)[1]),vector(mode="numeric",length=dim(SDSres)[1]))
names(SDSres)[13] <- "sSDS"
names(SDSres)[14] <- "pSDS"
for(i in 1:length(binM)){
	SDSres[SDSres[,12]==i,13] <- (SDSres[SDSres[,12]==i,10]-binM[i])/binSD[i]
	SDSres[SDSres[,12]==i,14] <- pnorm(SDSres[SDSres[,12]==i,13],lower.tail=F)
}

# Calculating FDR values from p-values
qobj <- qvalue(p = SDSres[,14])
SDSres[,15] <- qobj$qvalues
names(SDSres)[15] <- "qSDS"

# QQ values of normalised scores per bin
png(paste0('OutFigures/SDS_QQPlots_',fname,'N0.png'),width=4*(1+sqrt(5)),height=8,units = 'in',res=200)
par(mfrow=c(3,6))
for(i in 1:length(binM)){
	qqnorm(subset(SDSres,Bin==i)$sSDS, main=paste0("Normal Q-Q Plot, Bin ",i))
	qqline(subset(SDSres,Bin==i)$sSDS)
}
dev.off()

# Histogram of DAF
if(fname == "High")
{
	png(paste0('OutFigures/SDS_DAF_Hist.png'),width=8,height=8,units = 'in',res=200)
	hist(SDSres$DAF,breaks=seq(0.05,0.95,0.05),main=paste0("Histogram of Derived Allele Frequencies"),xlab="Derived Allele Frequency",xaxt="n")
	axis(1,at=seq(0.05,0.95,0.1),pos=(0))
	dev.off()
}

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
SDS_CO2_P <- SigSDSRes2[ which(SigSDSRes2[,13]==SDS_CO2),14]
cat(paste0(dim(SigSDSRes2)[1]," SNPs with FDR < 0.05\n"))
cat(paste0("Smallest SDS below FDR threshold: ",SDS_CO2,"\n"))

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

png(paste0('OutFigures/SDS_Chr_All_',fname,'N0.png'),width=4*(1+sqrt(5)),height=8,units = 'in',res=200)
par(mar=c(5,8.75,4,2) + 0.1)
ymax <- ceiling(max(-log10(SDSres[,14])))
plot(-log10(SDSres[,14]), col=colrs, pch=16, xaxt="n", yaxt="n", xlab="Chromosome",ylab="",main=substitute(bold(paste("sSDS Results for ",bolditalic('Bos taurus')," Autosomes"))),ylim=c(0,ymax))
axis(1, at=tickpos, labels=cno)
axis(2, at=seq(0,ymax,2), las=2)
text(x=-275000,y=(ymax/2),labels=paste("sSDS","\u2013Log10 P\u2013Value",sep="\n"),xpd=NA)
abline(-log10(alpha),0,lty=2)
abline(-log10(SDS_CO2_P),0,lty=3)
legend("topleft", legend=c(expression("Bonferroni\u2013"*"Corrected "*italic("P")*" < 0.05"), "FDR < 0.05"),col=c("red", "blue"),pch=16)
dev.off()

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

## Part 2: Are any high SDS scores within milk-producing regions?

# Start with fat percentage QTL
milkfQTL <- readRDS("MilkQTL/qtls_tbl_fat_curated.rds")
milkfQTL <- filter(milkfQTL,German_Holstein==French_Holstein)		# Filtering out those with conflicting effect sizes
milkfv <- QTLpol(milkfQTL,SDSres2)
nQmf <- dim(milkfv$Values)[1]
cat("For milk fat percentage:\n")
cat(nQmf,"QTLs have SDS scores assigned to them.\n")

# Next, protein content
milkpQTL <- readRDS("MilkQTL/qtls_tbl_prot_curated.rds")
milkpQTL <- filter(milkpQTL,German_Holstein==French_Holstein) %>% filter(chr!=25) %>% filter(`ARS-UCD2.1`!="deleted")		# Filtering out those with conflicting effect sizes; those on Chr 25; where not present in ARS background
milkpv <- QTLpol(milkpQTL,SDSres2)
nQmp <- dim(milkpv$Values)[1]
cat("For milk protein percentage:\n")
cat(nQmp,"QTLs have SDS scores assigned to them.\n")

# Histogram of distances to nearest SNP
png(paste0('OutFigures/QTL_Milk_dist_',fname,'N0.png'),width=4*(1+sqrt(5)),height=8,units = 'in',res=200)
par(mfrow=c(1,2))
hist(milkfv$Distances,main="Histogram of distances to nearest QTL\n(milk fat percentage)",xlab="Distance")
hist(milkpv$Distances,main="Histogram of distances to nearest QTL\n(milk protein content)",xlab="Distance")
dev.off()

milkfv$Values$Bin <- factor(milkfv$Values$Bin,levels=unique(milkfv$Values$Bin))
milkpv$Values$Bin <- factor(milkpv$Values$Bin,levels=unique(milkpv$Values$Bin))

cat("Spearman's rank test for milk traits, using all datapoints\n")
milkf_lm <- lm(sSDS ~ p,data=milkfv$Value)
milkp_lm <- lm(sSDS ~ p,data=milkpv$Value)
cor.test(milkfv$Value$p,milkfv$Value$sSDS,method="spearman",exact=F)
cor.test(milkpv$Value$p,milkpv$Value$sSDS,method="spearman",exact=F)

# Removing points with outlier p-value and repeating
milkfv2 <- subset(milkfv$Value,p<max(milkfv$Value$p))
milkpv2 <- subset(milkpv$Value,p<max(milkpv$Value$p))

cat("Spearman's rank test for milk traits, after removing high p-value\n")
milkf_lm2 <- lm(sSDS ~ p,data=milkfv2)
milkp_lm2 <- lm(sSDS ~ p,data=milkpv2)
cor.test(milkfv2$p,milkfv2$sSDS,method="spearman",exact=F)
cor.test(milkpv2$p,milkpv2$sSDS,method="spearman",exact=F)

# Part 3: Reading in Stature QTL data. First where effect sizes estimated in 6 of 7 Holstein
# These QTLs have positions relative to old assembly (UMD), the function transfers them to ARS-UCD assembly

statQTL6 <- readRDS("QTLStatureDat/curated_stature_snps_6_2020.rds")
stat6res <- QTLspol(statQTL6,SDSres2)
nQ26Q <- dim(stat6res$Values)[1]		# Number of QTL with SNPs assigned to them
s6Qmm <- filter(stat6res$Exact,QTLAncestral!=ANC)  # QTLs where listed ancestral state != new ancestral
cat("For stature QTLs with effects in 6 of 7 breeds:\n")
cat(nQ26Q,"QTLs have SDS scores assigned to them.\n")
cat(dim(s6Qmm)[1],"QTLs where anc SNP in database differs from current polarisation.\n")

# effect sizes in 5 of 7 Holstein
statQTL5 <- readRDS("QTLStatureDat/curated_stature_snps_5_2020.rds")
stat5res <- QTLspol(statQTL5,SDSres2)
nQ25Q <- dim(stat5res$Values)[1]		# Number of QTL with SNPs assigned to them
stat5res$Values$ANC <- factor(stat5res$Values$ANC,levels(stat5res$Values$QTLAncestral))
s5Qmm <- filter(stat5res$Exact,QTLAncestral!=ANC)  # QTLs where listed ancestral state != new ancestral
cat("For stature QTLs with effects in 5 of 7 breeds:\n")
cat(nQ25Q,"QTLs have SDS scores assigned to them.\n")
cat(dim(s5Qmm)[1],"QTLs where anc SNP in database differs from current polarisation.\n")

stat6res$Values$Bin <- factor(stat6res$Values$Bin,levels=unique(stat6res$Values$Bin))
stat5res$Values$Bin <- factor(stat5res$Values$Bin,levels=unique(stat5res$Values$Bin))

# Histogram of distances to nearest SNP
png(paste0('OutFigures/QTL_Stature_dist_',fname,'N0.png'),width=4*(1+sqrt(5)),height=8,units = 'in',res=200)
par(mfrow=c(1,2))
hist(stat6res$Distances,main="Histogram of distances to nearest QTL\n(stature, 6 of 7 breeds)",xlab="Distance")
hist(stat5res$Distances,main="Histogram of distances to nearest QTL\n(stature, 5 of 7 breeds)",xlab="Distance")
dev.off()

cat("Spearman's rank test for stature\n")
stat6_lm <- lm(sSDS ~ p,data=stat6res$Value)
stat5_lm <- lm(sSDS ~ p,data=stat5res$Value)
cor.test(stat6res$Value$p,stat6res$Value$sSDS,method="spearman",exact=F)
cor.test(stat5res$Value$p,stat5res$Value$sSDS,method="spearman",exact=F)

# Plots of correlations between sSDS, QTL P-values
# Seperate setups for high, low N0
pfunc <- function(data,mt,xp,ma,cin)
{
	par(mar=c(5,7.5,4,2) + 0.1)
	plot(sSDS~p,data=data,xlab="Absolute Log10 QTL P-Value",ylab="",main=ma,pch=16,las=1,cex.lab=cin,cex.axis=cin,cex.main=cin*1.25)
	ymp <- min(data$sSDS)+(max(data$sSDS)-min(data$sSDS))/2
	text(x=xp,y=ymp,labels=paste("tSDS",sep="\n"),xpd=NA,cex=cin)
	abline(mt)
}

if(ind==1){
	png(paste0('OutFigures/QTL_All_SDS_',fname,'N0.png'),width=12,height=12,units = 'in',res=200)
	par(mfrow=c(2,2))
	pfunc(milkfv$Values,milkf_lm,-10,"(a) Milk fat percentage",1)
	pfunc(milkpv$Values,milkp_lm,-45,"(b) Milk protein percentage",1)
	pfunc(stat6res$Values,stat6_lm,4.5,"(c) Stature, determined from 6 of 7 breeds",1)
	pfunc(stat5res$Values,stat5_lm,4,"(d) Stature, determined from 5 of 7 breeds",1)
	dev.off()
	
	# Plot of correlations, milk data without outlier p-value
	png(paste0('OutFigures/QTL_Milk_SDS_NoOutlier_',fname,'N0.png'),width=12,height=6,units = 'in',res=200)
	par(mfrow=c(1,2))
	pfunc(milkfv2,milkf_lm2,1.5,"(a) Milk fat percentage",1)
	pfunc(milkpv2,milkp_lm2,-15,"(b) Milk protein percentage",1)
	dev.off()
}else if(ind==2){
	png(paste0('OutFigures/QTL_All_SDS_',fname,'N0.png'),width=12,height=18,units = 'in',res=200)
	par(mfrow=c(3,2))
	pfunc(milkfv$Values,milkf_lm,-7.5,"(a) Milk fat percentage",1.5)
	pfunc(milkpv$Values,milkp_lm,-35,"(b) Milk protein percentage",1.5)
	pfunc(milkfv2,milkf_lm2,4.15,"(c) Milk fat percentage, outlier P-value removed",1.5)
	pfunc(milkpv2,milkp_lm2,-10,"(d) Milk protein percentage, outlier P-value removed",1.5)
	pfunc(stat6res$Values,stat6_lm,5,"(e) Stature, determined from 6 of 7 breeds",1.5)
	pfunc(stat5res$Values,stat5_lm,4.5,"(f) Stature, determined from 5 of 7 breeds",1.5)
	dev.off()
}

# Diagnostic plots for LM fits
png(paste0('OutFigures/QTL_SDS_LM_milkf_',fname,'N0.png'),width=12,height=12,units = 'in',res=200)
par(mfrow=c(2,2))
plot(milkf_lm)
dev.off()

png(paste0('OutFigures/QTL_SDS_LM_milkp_',fname,'N0.png'),width=12,height=12,units = 'in',res=200)
par(mfrow=c(2,2))
plot(milkp_lm)
dev.off()

# Diagnostic plots for LM fits
png(paste0('OutFigures/QTL_SDS_LM_milkf2_',fname,'N0.png'),width=12,height=12,units = 'in',res=200)
par(mfrow=c(2,2))
plot(milkf_lm2)
dev.off()

png(paste0('OutFigures/QTL_SDS_LM_milkp2_',fname,'N0.png'),width=12,height=12,units = 'in',res=200)
par(mfrow=c(2,2))
plot(milkp_lm2)
dev.off()

png(paste0('OutFigures/QTL_SDS_LM_stat6_',fname,'N0.png'),width=12,height=12,units = 'in',res=200)
par(mfrow=c(2,2))
plot(stat6_lm)
dev.off()

png(paste0('OutFigures/QTL_SDS_LM_stat5_',fname,'N0.png'),width=12,height=12,units = 'in',res=200)
par(mfrow=c(2,2))
plot(stat5_lm)
dev.off()

# EOF