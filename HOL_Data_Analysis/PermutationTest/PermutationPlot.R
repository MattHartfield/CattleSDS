# 12th May 2021
# Analysing new permutation test results

library(tidyverse)
library(cowplot)

setwd("/Users/hartfield/Documents/MilkSDS/HOL_Data_Analysis/PermutationTest")

# Reading in tables
datH <- read_table2("CompiledRes/HOL_Permutation_HighN0.dat")
datL <- read_table2("CompiledRes/HOL_Permutation_LowN0.dat")

# Matrix of actual F-values
actpv <- rbind(datH[1,],datL[1,])
#row.names(actpv) <- cbind("HighN0","LowN0")
#actpv <- data.frame(MilkFatF=c(0.5847,0.3616),MilkProtF=c(0.5716,0.1405),Stat6F=c(0.28,0.4408),Stat5F=c(0.04794,0.08695))

# Combining permutations into one table
datH <- as_tibble(cbind("HighN0",datH[2:1001,]))
datL <- as_tibble(cbind("LowN0",datL[2:1001,]))
colnames(datH)[1] <- "N0"
colnames(datL)[1] <- "N0"
dat <- rbind(datH,datL)

# Now to plot
# Milk fat content
tpv1 <-  as_tibble(cbind(unique(dat$N0),actpv[,1]))
colnames(tpv1) <- c("N0","Fval")
PmfpH <- dat %>% filter(N0=="HighN0") %>% summarize(PmfpH = sum(MilkFatF > as.double(tpv1[1,2]))/1000) %>% as.double
PmfpL <- dat %>% filter(N0=="LowN0") %>% summarize(PmfpL = sum(MilkFatF > as.double(tpv1[2,2]))/1000) %>% as.double
tpvt1 <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P-value = ", c(PmfpH, PmfpL)),x=15,y=60)
p1 <- ggplot(dat,aes(x=MilkFatF)) + geom_histogram(binwidth=0.025) + facet_wrap(~N0) + xlab("Permuted F-values, Milk Fat Content") + ylab("Count") + geom_vline(data=tpv1,aes(xintercept=Fval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt1,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5))

# Milk protein content
tpv2 <-  as_tibble(cbind(unique(dat$N0),actpv[,2]))
colnames(tpv2) <- c("N0","Fval")
tpv2$Fval <- as.numeric(as.character(tpv2$Fval))
PmppH <- dat %>% filter(N0=="HighN0") %>% summarize(PmppH = sum(MilkProtF > as.double(tpv2[1,2]))/1000) %>% as.double
PmppL <- dat %>% filter(N0=="LowN0") %>% summarize(PmppL = sum(MilkProtF > as.double(tpv2[2,2]))/1000) %>% as.double
tpvt2 <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P-value = ", c(PmppH, PmppL)),x=12,y=70)
p2 <- ggplot(dat,aes(x=MilkProtF)) + geom_histogram(binwidth=0.025) + facet_wrap(~N0) + xlab("Permuted F-values, Milk Protein Content") + ylab("Count") + geom_vline(data=tpv2,aes(xintercept=Fval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt2,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5))

# Stature 6
tpv3 <-  as_tibble(cbind(unique(dat$N0),actpv[,3]))
colnames(tpv3) <- c("N0","Fval")
Ps6pH <- dat %>% filter(N0=="HighN0") %>% summarize(Ps6pH = sum(Stat6F > as.double(tpv3[1,2]))/1000) %>% as.double
Ps6pL <- dat %>% filter(N0=="LowN0") %>% summarize(Ps6pL = sum(Stat6F > as.double(tpv3[2,2]))/1000) %>% as.double
tpvt3 <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P-value = ", c(Ps6pH, Ps6pL)),x=8.5,y=70)
p3 <- ggplot(dat,aes(x=Stat6F)) + geom_histogram(binwidth=0.025) + facet_wrap(~N0) + xlab("Permuted F-values, Stature (6 of 7 breeds)") + ylab("Count") + geom_vline(data=tpv3,aes(xintercept=Fval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt3,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5))

# Stature 5
tpv4 <-  as_tibble(cbind(unique(dat$N0),actpv[,4]))
colnames(tpv4) <- c("N0","Fval")
tpv4$Fval <- as.numeric(as.character(tpv4$Fval))
Ps5pH <- dat %>% filter(N0=="HighN0") %>% summarize(Ps5pH = sum(Stat5F > as.double(tpv4[1,2]))/1000) %>% as.double
Ps5pL <- dat %>% filter(N0=="LowN0") %>% summarize(Ps5pL = sum(Stat5F > as.double(tpv4[2,2]))/1000) %>% as.double
tpvt4 <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P-value = ", c(Ps5pH, Ps5pL)),x=8.5,y=75)
p4 <- ggplot(dat,aes(x=Stat5F)) + geom_histogram(binwidth=0.025) + facet_wrap(~N0) + xlab("Permuted F-values, Stature (5 of 7 breeds)") + ylab("Count") + geom_vline(data=tpv4,aes(xintercept=Fval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt4,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5))

# All together!
gr <- (1+sqrt(5))/2
baseh = 12
pdf("PermutationPlots.pdf",height=baseh,width=baseh*gr)
plot_grid(p1,p2,p3,p4,labels=c("A","B","C","D"),label_size=24)
dev.off()

# EOF