# 22nd May 2019
# R Script for analysing SDS results, both high and low N0 results

# 7th Dec 2020
# Permutation test code

# First decide which dataset to use
# Set up so one can run from command line
# Otherwise one can set 'ind' manually
# 1 = High N0 case; 2 = Low N0 case

args <- commandArgs(trailingOnly = TRUE)
ind <- as.integer(args[1])
idx <- as.integer(args[2])
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
	apol <- read.table(paste0("../SDSAll_Polarisation.dat"),header=T)
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
	apol <- read.table(paste0("../SDSAll_Polarisation.dat"),header=T)
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
theseed <- sample(2147483647-1,1)
set.seed(theseed)
cat("Seed is ", theseed, "\n",sep="")
#setwd("/Users/hartfield/Documents/MilkSDS/HOL_Data_Analysis")
setwd("/usr/home/qgg/mhart/MilkSDS_Dec18/PermutationTest")
SDSres <- read.table(paste0("../SDSAll_",fname,"N0.dat"),header=T)
#SDSres <- read.table(paste0("SDSAll_",fname,"N0.dat"),header=T)
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

## Part 1: Finding regions with significantly high absolute SDS.
alpha <- (0.05/dim(SDSres)[1])
SigSDS <- intersect(which((SDSres[,14] < alpha)),which((as.numeric(as.character(SDSres[,9])) >= 1)))
SDSres[SigSDS,] -> SigSDSRes
SDS_CO <- min(SigSDSRes[,13])

# Extracting significant results considering q < 0.05 (FDR)
SigSDS2 <- intersect(which((SDSres[,15] < 0.05)),which((as.numeric(as.character(SDSres[,9])) >= 1)))
SDSres[SigSDS2,] -> SigSDSRes2
SDS_CO2 <- min(SigSDSRes2[,13])
SDS_CO2_P <- SigSDSRes2[ which(SigSDSRes2[,13]==SDS_CO2),14]

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

# Next, protein content
milkpQTL <- readRDS("MilkQTL/qtls_tbl_prot_curated.rds")
milkpQTL <- filter(milkpQTL,German_Holstein==French_Holstein) %>% filter(chr!=25) %>% filter(`ARS-UCD2.1`!="deleted")		# Filtering out those with conflicting effect sizes; those on Chr 25; where not present in ARS background
milkpv <- QTLpol(milkpQTL,SDSres2)

# Part 3: Reading in Stature QTL data. First where effect sizes estimated in 6 of 7 Holstein
# These QTLs have positions relative to old assembly (UMD), the function transfers them to ARS-UCD assembly
statQTL6 <- readRDS("QTLStatureDat/curated_stature_snps_6_2020.rds")
stat6res <- QTLspol(statQTL6,SDSres2)

# effect sizes in 5 of 7 Holstein
statQTL5 <- readRDS("QTLStatureDat/curated_stature_snps_5_2020.rds")
stat5res <- QTLspol(statQTL5,SDSres2)

# Permutation test code
list_dat <- list(milkfv$Values, milkpv$Values, stat6res$Values, stat5res$Values)
if(idx == 1){
	output_f2 <- vector(mode="numeric",length=length(list_dat))
	for(l in 1:length(list_dat))
	{
		tttab <- list_dat[[l]]
		lmtt <- summary(lm(sSDS ~ p,data=tttab))
		output_f2[l] <- lmtt$fstatistic[1]
	}
	write.table(t(output_f2),file=paste0("OutTables/PermTest_n",idx,"_",fname,"N0.dat"), quote = F, row.names = F, col.names = F)
}
output_f <- vector(mode="numeric",length=length(list_dat))
for(l in 1:length(list_dat))
{
	ttab <- list_dat[[l]]
	nSNPs <- dim(ttab)[1]
	permSDS <- data.frame(sSDS=numeric(),p=numeric())
	for(k in 1:nSNPs)
	{
		Chrk <- ttab[k,1]
		Bink <- ttab[k,4]
		sdsk <- as.numeric(SDSres2 %>% filter(CHROMOSOME==Chrk,Bin==Bink) %>% select(sSDS) %>% sample_n(1))
		# Now deciding whether to change polarisation of score 
		ispol <- rowSums(ttab[k,c("NotRef","RePol")])
		if(ispol==1){
			sdsk <- (-1)*sdsk
		}
		permSDS[k,1] <- sdsk
	}
	permSDS[,2] <- ttab$p
	lmt <- summary(lm(sSDS ~ p,data=permSDS))
	output_f[l] <- lmt$fstatistic[1]
}

if(idx == 1)
{
	write.table(t(output_f),file=paste0("OutTables/PermTest_n",idx,"_",fname,"N0.dat"), quote = F, row.names = F, col.names = F, append = T)
}else
{
	write.table(t(output_f),file=paste0("OutTables/PermTest_n",idx,"_",fname,"N0.dat"), quote = F, row.names = F, col.names = F)
}

# EOF