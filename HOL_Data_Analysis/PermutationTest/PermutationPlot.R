# 12th May 2021
# Analysing new permutation test results

# 28th Sept 2021
# Updated to consider rho correlations; add milk cases without outlier

library(tidyverse)
library(cowplot)

setwd("/Users/hartfield/Documents/MilkSDS/HOL_Data_Analysis/PermutationTest")

# Reading in tables
datH <- read_table2("CompiledRes/HOL_Permutation_HighN0.dat")
datL <- read_table2("CompiledRes/HOL_Permutation_LowN0.dat")

# Matrix of actual rho values
actpv <- rbind(datH[1,],datL[1,])
#row.names(actpv) <- cbind("HighN0","LowN0")
#actpv <- data.frame(MilkFatR=c(0.5847,0.3616),MilkProtR=c(0.5716,0.1405),Stat6R=c(0.28,0.4408),Stat5R=c(0.04794,0.08695))

# Combining permutations into one table
datH <- as_tibble(cbind("HighN0",datH[2:1001,]))
datL <- as_tibble(cbind("LowN0",datL[2:1001,]))
colnames(datH)[1] <- "N0"
colnames(datL)[1] <- "N0"
dat <- rbind(datH,datL)

# Now to plot
# Milk fat content
tpv1 <-  as_tibble(cbind(unique(dat$N0),actpv[,1]))
colnames(tpv1) <- c("N0","Rval")
tpv1 <- dat %>% select(N0,MilkFatR) %>% group_by(N0) %>% summarize(meanMFR=mean(MilkFatR)) %>% full_join(tpv1,.)
tpv1 <- dat %>% select(N0,MilkFatR) %>% group_by(N0) %>% summarize(sdMFR=sd(MilkFatR)) %>% full_join(tpv1,.)
tpv1 <- tpv1 %>% group_by(N0) %>% summarise(zR = (Rval-meanMFR)/sdMFR) %>% inner_join(tpv1,.)
tpv1 <- tpv1 %>% group_by(N0) %>% summarise(plab = ifelse(zR>=0,2*pnorm(zR,lower.tail=F),2*pnorm(zR,lower.tail=T))) %>% inner_join(tpv1,.)
tpvt1 <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P = ", sprintf("%.3f",tpv1$plab)),x=-0.425,y=40)
p1 <- ggplot(dat,aes(x=MilkFatR)) + geom_histogram(binwidth=0.025,fill="gray60") + facet_wrap(~N0) + xlab(bquote("Permuted"~rho~"values, Milk Fat Content")) + ylab("Count") + geom_vline(data=tpv1,aes(xintercept=Rval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt1,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5),plot.margin=margin(r=1,unit="cm"))

# Milk protein content
tpv2 <-  as_tibble(cbind(unique(dat$N0),actpv[,2]))
colnames(tpv2) <- c("N0","Rval")
tpv2 <- dat %>% select(N0,MilkProtR) %>% group_by(N0) %>% summarize(meanMPR=mean(MilkProtR)) %>% full_join(tpv2,.)
tpv2 <- dat %>% select(N0,MilkProtR) %>% group_by(N0) %>% summarize(sdMPR=sd(MilkProtR)) %>% full_join(tpv2,.)
tpv2 <- tpv2 %>% group_by(N0) %>% summarise(zR = (Rval-meanMPR)/sdMPR) %>% inner_join(tpv2,.)
tpv2 <- tpv2 %>% group_by(N0) %>% summarise(plab = ifelse(zR>=0,2*pnorm(zR,lower.tail=F),2*pnorm(zR,lower.tail=T))) %>% inner_join(tpv2,.)
tpvt2 <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P = ", sprintf("%.3f",tpv2$plab)),x=-0.2,y=75)
p2 <- ggplot(dat,aes(x=MilkProtR)) + geom_histogram(binwidth=0.025,fill="gray60") + facet_wrap(~N0) + xlab(bquote("Permuted"~rho~"values, Milk Protein Content")) + ylab("Count") + geom_vline(data=tpv2,aes(xintercept=Rval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt2,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5),plot.margin=margin(r=1,unit="cm"))

# Milk fat content without outlier
tpv1a <-  as_tibble(cbind(unique(dat$N0),actpv[,3]))
colnames(tpv1a) <- c("N0","Rval")
tpv1a <- dat %>% select(N0,MilkFatNoOutR) %>% group_by(N0) %>% summarize(meanMFNR=mean(MilkFatNoOutR)) %>% full_join(tpv1a,.)
tpv1a <- dat %>% select(N0,MilkFatNoOutR) %>% group_by(N0) %>% summarize(sdMFNR=sd(MilkFatNoOutR)) %>% full_join(tpv1a,.)
tpv1a <- tpv1a %>% group_by(N0) %>% summarise(zR = (Rval-meanMFNR)/sdMFNR) %>% inner_join(tpv1a,.)
tpv1a <- tpv1a %>% group_by(N0) %>% summarise(plab = ifelse(zR>=0,2*pnorm(zR,lower.tail=F),2*pnorm(zR,lower.tail=T))) %>% inner_join(tpv1a,.)
tpvt1a <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P = ", sprintf("%.3f",tpv1a$plab)),x=-0.4,y=40)
p1a <- ggplot(dat,aes(x=MilkFatNoOutR)) + geom_histogram(binwidth=0.025,fill="gray60") + facet_wrap(~N0) + xlab(bquote("Permuted"~rho~"values, Milk Fat Content (No Outlier)")) + ylab("Count") + geom_vline(data=tpv1a,aes(xintercept=Rval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt1a,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5),plot.margin=margin(r=1,unit="cm"))

# Milk protein content without outlier
tpv2a <-  as_tibble(cbind(unique(dat$N0),actpv[,4]))
colnames(tpv2a) <- c("N0","Rval")
tpv2a <- dat %>% select(N0,MilkProtNoOutR) %>% group_by(N0) %>% summarize(meanMPNR=mean(MilkProtNoOutR)) %>% full_join(tpv2a,.)
tpv2a <- dat %>% select(N0,MilkProtNoOutR) %>% group_by(N0) %>% summarize(sdMPNR=sd(MilkProtNoOutR)) %>% full_join(tpv2a,.)
tpv2a <- tpv2a %>% group_by(N0) %>% summarise(zR = (Rval-meanMPNR)/sdMPNR) %>% inner_join(tpv2a,.)
tpv2a <- tpv2a %>% group_by(N0) %>% summarise(plab = ifelse(zR>=0,2*pnorm(zR,lower.tail=F),2*pnorm(zR,lower.tail=T))) %>% inner_join(tpv2a,.)
tpvt2a <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P = ", sprintf("%.3f",tpv2a$plab)),x=-0.2,y=75)
p2a <- ggplot(dat,aes(x=MilkProtNoOutR)) + geom_histogram(binwidth=0.025,fill="gray60") + facet_wrap(~N0) + xlab(bquote("Permuted"~rho~"values, Milk Protein Content (No Outlier)")) + ylab("Count") + geom_vline(data=tpv2a,aes(xintercept=Rval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt2a,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5),plot.margin=margin(r=1,unit="cm"))

# Stature 6
tpv3 <-  as_tibble(cbind(unique(dat$N0),actpv[,5]))
colnames(tpv3) <- c("N0","Rval")
tpv3 <- dat %>% select(N0,Stat6R) %>% group_by(N0) %>% summarize(meanS6R=mean(Stat6R)) %>% full_join(tpv3,.)
tpv3 <- dat %>% select(N0,Stat6R) %>% group_by(N0) %>% summarize(sdS6R=sd(Stat6R)) %>% full_join(tpv3,.)
tpv3 <- tpv3 %>% group_by(N0) %>% summarise(zR = (Rval-meanS6R)/sdS6R) %>% inner_join(tpv3,.)
tpv3 <- tpv3 %>% group_by(N0) %>% summarise(plab = ifelse(zR>=0,2*pnorm(zR,lower.tail=F),2*pnorm(zR,lower.tail=T))) %>% inner_join(tpv3,.)
tpvt3 <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P = ", sprintf("%.3f",tpv3$plab)),x=-0.35,y=60)
p3 <- ggplot(dat,aes(x=Stat6R)) + geom_histogram(binwidth=0.025,fill="gray60") + facet_wrap(~N0) + xlab(bquote("Permuted"~rho~"values, Stature (6 of 7 breeds)")) + ylab("Count") + geom_vline(data=tpv3,aes(xintercept=Rval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt3,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5),plot.margin=margin(r=1,unit="cm"))

# Stature 5
tpv4 <-  as_tibble(cbind(unique(dat$N0),actpv[,6]))
colnames(tpv4) <- c("N0","Rval")
tpv4 <- dat %>% select(N0,Stat5R) %>% group_by(N0) %>% summarize(meanS5R=mean(Stat5R)) %>% full_join(tpv4,.)
tpv4 <- dat %>% select(N0,Stat5R) %>% group_by(N0) %>% summarize(sdS5R=sd(Stat5R)) %>% full_join(tpv4,.)
tpv4 <- tpv4 %>% group_by(N0) %>% summarise(zR = (Rval-meanS5R)/sdS5R) %>% inner_join(tpv4,.)
tpv4 <- tpv4 %>% group_by(N0) %>% summarise(plab = ifelse(zR>=0,2*pnorm(zR,lower.tail=F),2*pnorm(zR,lower.tail=T))) %>% inner_join(tpv4,.)
tpvt4 <- data.frame(N0=c("HighN0","LowN0"),plab=paste0("P = ", sprintf("%.3f",tpv4$plab)),x=-0.3,y=60)
p4 <- ggplot(dat,aes(x=Stat5R)) + geom_histogram(binwidth=0.025,fill="gray60") + facet_wrap(~N0) + xlab(bquote("Permuted"~rho~"values, Stature (5 of 7 breeds)")) + ylab("Count") + geom_vline(data=tpv4,aes(xintercept=Rval),linetype="dotted", size=1.25) + theme_bw(base_size=24) + geom_text(data=tpvt4,aes(x=x,y=y,label=plab),size=6) + theme(axis.title.y=element_text(angle=0,vjust=0.5),plot.margin=margin(r=1,unit="cm"))

# All together!
gr <- (1+sqrt(5))/2
baseh = 12
pdf("PermutationPlots.pdf",height=baseh,width=baseh*gr)
plot_grid(p1,p2,p1a,p2a,p3,p4,labels=c("A","B","C","D","E","F"),label_size=24,nrow=3,ncol=2,byrow=T)
dev.off()

# EOF