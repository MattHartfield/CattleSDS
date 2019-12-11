# 11th Dec 2019
# R Script for carrying out permutation test on milk-protein genes

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

# Obtaining (1) starting point of each chromosome (2) length of chromosome
# For permutation analyses, to determine where to place SNPs
minp <- vector(mode="numeric",length=length(cno))
lengthp <- vector(mode="numeric",length=length(cno))
names(minp) <- factor(orderedChr)
names(lengthp) <- factor(orderedChr)
for(i in unique(SDSres$CHROMOSOME)){
	SDST <- subset(SDSres,CHROMOSOME==i)
	minp[i] <- min(SDST$POS)
	lengthp[i] <- max(SDST$POS)-min(SDST$POS)
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

# Has masking worked? Should give integer(0)
intersect(which((SDSres2[,14] < alpha)),which((as.numeric(as.character(SDSres2[,9])) >= 1)))

# Reading in and classifying milk-producing regions
# From milk-protein gene set (additional data file 2 from Lemay et al 2009 Genome Biology)
# Then including milk-gene as a factor. Putting scores in bins, determining effect of milk on SDS.
milky <- read.table("Milk_Protein_Genes_Sept19/milk_protein_names_conv_Sept19.bed",skip=1)[,1:3] # Only including gene locations and autosomal genes
names(milky) <- c("chrom","chromStart","chromEnd")
milky$chrom <- factor(milky$chrom, levels=unique(SDSres$CHROMOSOME))
milkreg <- vector(mode="numeric",length=0)
muse <- vector(mode="numeric",length=length(sort(unique(milky$chrom))))
#mlen <- vector(mode="numeric",length=length(sort(unique(milky$chrom))))
glen <- as.data.frame(matrix(nrow=0,ncol=2))
names(muse) <- (unique(milky$chrom))
#names(mlen) <- (unique(milky$chrom))
for(i in unique(milky$chrom)){
	hpss <- as.matrix(subset(milky,chrom==i)[2:3])
	SDST <- subset(SDSres2,CHROMOSOME==i)
	listres <- apply(hpss,1,function(x) row.names(SDST[intersect(which(SDST[SDST$CHROMOSOME==i,5] >= x[1]), which(SDST[SDST$CHROMOSOME==i,5] <= x[2])),]))
	if(length(listres) > 0){
		if(dim(hpss)[1]!=1){
			muse[names(muse)==i] <- sum(lapply(listres,length)!=0)
#			mlen[names(mlen)==i] <- sum(hpss[which(lapply(listres,length)!=0),2]-hpss[which(lapply(listres,length)!=0),1])

			# To check and implement for hpss == 1
			glen <- rbind(glen,cbind(as.data.frame.factor(i),hpss[which(lapply(listres,length)!=0),2]-hpss[which(lapply(listres,length)!=0),1]))
		}else if(dim(hpss)[1]==1){
			muse[names(muse)==i] <- sum(length(listres)!=0)
			mlen[names(mlen)==i] <- hpss[which(length(listres)!=0),2]-hpss[which(length(listres)!=0),1]
		}
	}
	milkreg<-c(milkreg,unlist(listres))
}
names(milkreg) <- NULL
names(glen) <- c("chrom","geneLength")

# Adding column stating if SNPs lie within a milk-producing region
SDSres2 <- cbind(SDSres2,vector(mode="numeric",length=dim(SDSres2)[1]))
names(SDSres2)[16] <- "Milk"
SDSres2[,16] <- 0
SDSres2[row.names(SDSres2)%in%milkreg,16] <- 1
SDSres2[!(row.names(SDSres2)%in%milkreg),16] <- 0
SDSres2$Milk <- factor(SDSres2$Milk,levels=unique(SDSres2$Milk))

# General Linear model fit, comparing milk and non-milk regions, Gamma link function. Takes a few minutes to run for each model fit
SDSres2$CHROMOSOME <- factor(SDSres2$CHROMOSOME,levels=unique(SDSres2$CHROMOSOME))
SDSres2$Bin <- factor(SDSres2$Bin,levels=unique(SDSres2$Bin))
SDSres2$Bin <- relevel(SDSres2$Bin,ref="1")
SDSres2$Milk <- relevel(SDSres2$Milk,ref="0")
fitg <- glm(sSDS ~ CHROMOSOME + Bin + Milk,data=SDSres2,family=Gamma(link="inverse"));summary(fitg)
fitg0 <- glm(sSDS ~ CHROMOSOME + Bin,data=SDSres2,family=Gamma(link="inverse"));summary(fitg0)
anova(fitg0,fitg,test="Chisq")
pvalA <- anova(fitg0,fitg)$"F"[2]

# Print summary to file
sink(file=paste0("OutTables/Milk_GLM_Summary_",fname,"N0_FromPermutationTest.txt"))
summary(fitg)
anova(fitg0,fitg,test="Chisq")
sink()

# Now permutation analyses
perms <- 500
mtot <- sum(muse)
bins <- round(mean(mlen)/2)	# Note factor of two since we look both up- and downstream of the target SNPs
#bins <- 12500
pvals <- vector(mode="numeric",length=perms)
for(j in 1:perms){
	nsnps <- rmultinom(1, mtot, prob = (lengthp - 2*bins) )		# Number of pseudo SNPs on each region 
	# Creating pseudo-stature table
	for(i in unique(SDSres$CHROMOSOME)){
		SDST2 <- data.frame(cbind(rep(i,nsnps[row.names(nsnps)==i,]),sort(floor(runif(nsnps[row.names(nsnps)==i,],min=(minp[i]+bins),max=(minp[i]+lengthp[i]-2*bins))))))
		names(SDST2) <- c("CHR","BP")
		SDST2$CHR <- factor(SDST2$CHR,levels=levels(SigSDSRes2$CHROMOSOME))
		SDST2$BP <- as.numeric(as.character(SDST2$BP))
		if(i=="Chr1"){
			milkyP <- SDST2
		}else if(i!="Chr1"){
			milkyP <- rbind(milkyP, SDST2)
		}
	}

	# Classifying regions based on random SNPs
	milkreg <- vector(mode="numeric",length=0)
	for(i in sort(unique(milkyP$CHR))){
		hpss <- as.matrix(cbind(subset(milkyP,CHR==i)[2]-bins,subset(milkyP,CHR==i)[2]+bins))
		SDST <- subset(SDSres2,CHROMOSOME==i)
		milkreg<-c(milkreg,unlist(apply(hpss,1,function(x) row.names(SDST[intersect(which(SDST[SDST$CHROMOSOME==i,5] >= x[1]), which(SDST[SDST$CHROMOSOME==i,5] <= x[2])),]))))
	}
	names(milkreg) <- NULL
	SDSres3 <- SDSres2
	SDSres3[row.names(SDSres3)%in%milkreg,16] <- 1
	SDSres3[!(row.names(SDSres3)%in%milkreg),16] <- 0
	SDSres3$Milk <- factor(SDSres3$Milk,levels=unique(SDSres3$Milk))
	SDSres3$Bin <- factor(SDSres3$Bin,levels=unique(SDSres3$Bin))
	SDSres3$Milk <- relevel(SDSres3$Milk,ref="0")
	fit <- lm(sSDS ~ CHROMOSOME + Bin + Milk,data=SDSres3)
	fit0 <- lm(sSDS ~ CHROMOSOME + Bin,data=SDSres3)
	pvals[j] <- anova(fit0,fit)$"F"[2]	# P-value of comparison
}

write.table(c(pvalA,pvals),file=paste0("PermTest_n",perms,"_",fname,"N0_Milk.dat"), quote = F, row.names = F, col.names = F)

# EOF