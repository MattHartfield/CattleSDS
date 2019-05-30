"Script that:
1. Describes the data:
-plots the cM with physical distance 
-plots the cumulative cM info that is captured in
each BinID
-plots the physical distance info that is captured in
each BinID
2. Obtains the derivative of the cubic regression of 
each chromosome
3. For each midbin position, obtain the local recombi rate.
"

setwd("/home/jimenez/Documents/Aarhus/cluster_foulum/Analysis_25_SNPs/LRR/LocalRR/")
getwd()

#### TEST: one chromosome (Chromosome 1)
# Import the data:
chrD=read.table(file ="../Indexing/filtered_recombination_data_withBin")

names(chrD)=c("R_ID", "Chrm","PhysD","Pos3.0+1","Index",
              "BTAU","BTAU+1","cM","ID1","Database",
              "ID2","Type","Informative_Meiosis",
              "Kosambi_distance","BinID")

head(chrD)

chrD <- chrD[chrD[,2]==23,]

# Physical distance and cM
plot(chrD$PhysD,chrD$cM, pch=20, cex=0.25, col="grey30")
x11()
# Coverage of the bins (on cM information):
plot(chrD$cM,chrD$BinID, pch=20, cex=0.25, col="grey30",
     xlab="cM", ylab="Bin ID", main="How much cM info per bin")

# Coverage of the bins (on Physical Distance):
plot(chrD$PhysD,chrD$BinID, pch=20, cex=0.25, col="grey30",
     xlab="Phys.Distance", ylab="Bin ID", main="Coverage of the bins")

plot()
# In previous scripts,we tested that Cubic Regression 
# is the model that fits the data best:
CubReg=lm(cM~PhysD + I(PhysD^2) + I(PhysD^3), data=chrD)
summary(CubReg)

# Derivative of my equation:
CubReg$coefficients
a = CubReg$coefficients[0]
b = CubReg$coefficients[1]
c = CubReg$coefficients[2]
d = CubReg$coefficients[3]

# We need to take the derivative of the cubic regression formula:
CubRegDeriv = (D(expression(a+b*x + c*x^2 + d*x^3), 'x'))

# Check how the derivate of the cubic regression looks like:
numbers = seq(0,10^8,10^3)
CubRegDeriv_Bin = eval(D(expression(a+b*x + c*x^2 + d*x^3), 'x'),
                       list(x=numbers))

plot(CubRegDeriv_Bin)
########################### 

#################################"


# DATA
chrD=read.table(file ="../Indexing/filtered_recombination_data_withBin")

names(chrD)=c("R_ID", "Chrm","PhysD","Pos3.0+1","Index",
              "BTAU","BTAU+1","cM","ID1","Database",
              "ID2","Type","Informative_Meiosis",
              "Kosambi_distance","BinID")

head(chrD)
##############
# Loop over chromosomes/bins:

pdf("local_recombination_rate.pdf")
par(mfrow=c(2,2))
# Step 1: per chromosome...
nb_chrm = c(1:29)

for (chrm in nb_chrm){
  
  # Selection of data from chrom:
  chrD_selectionChrm = chrD[chrD[,2]==chrm,]
  nb_data_in_local_rr <- table(chrD_selectionChrm[,15])
  # Physical distance and cM
  #plot(chrD$PhysD,chrD$cM, pch=20, cex=0.25, col="grey30")
  # Coverage of the bins (on cM information):
  plot(chrD_selectionChrm$cM,chrD_selectionChrm$BinID, pch=20, cex=0.25, col="grey30",
       xlab="cM", ylab="Bin ID", main="How much cM info per bin")
#   plot(chrD_selectionChrm$cM,chrD_selectionChrm$BinID, pch=20, cex=0.25, col="grey30",
#        xlab="cM", ylab="Bin ID", main="2How much cM info per bin",yaxt='n')
# axis(side = 2, at=1:(length(nb_bins_in_local_rr)), labels = F)
# text(seq(1,nb_bins_in_local_rr,by=1), par("usr")[3]-0.05, 
#      srt = 60, adj=1, xpd = TRUE,
#      labels = nb_bins_in_local_rr, cex=0.50)
  # Coverage of the bins (on Physical Distance):
  #plot(chrD$PhysD,chrD$BinID, pch=20, cex=0.25, col="grey30",
       #xlab="Phys.Distance", ylab="Bin ID", main="Coverage of the bins")
  
  # In previous scripts,we tested that Cubic Regression 
  # is the model that fits the data best:
  CubReg=lm(cM~PhysD + I(PhysD^2) + I(PhysD^3), 
            data=chrD_selectionChrm)
  
  summary(CubReg)
  
  # Derivative of my equation:
  CubReg$coefficients
  a = CubReg$coefficients[1]
  b = CubReg$coefficients[2]
  c = CubReg$coefficients[3]
  d = CubReg$coefficients[4]
  
  print(a)
  print(b)
  print(c)
  print(d)
# Step 2: Per bin...
nb_bins = c(1:max(unique(chrD_selectionChrm$BinID))) # we are not interested in
                                       # binID==0 
list_recombin_rate = 0
list_midbins = 0
list_recombin_rate_check=0
  
for (bin in nb_bins){
  print(bin)
  chrD_selectionBin= chrD_selectionChrm[chrD_selectionChrm$BinID==bin,]
  
  x1 = min(chrD_selectionBin$PhysD)
  print(x1)
  x2 = max(chrD_selectionBin$PhysD)
  print(x2)
  
  L = x2-x1
  
  Midbin_perBin = x1+(L/2)
  print(Midbin_perBin)
  # Midbin_perBin is the same as (x1+x2)/2 # Wrong!
  
  list_midbins[bin] = Midbin_perBin
  cubic_regression_formula = expression(a + b*x + c*x^2 + d*x^3)
  # Substituting the midbin point in the derivate function:
  CubRegDeriv_Bin = eval(D(cubic_regression_formula, 'x'),
                         list(x=Midbin_perBin))
  
  cub_check = b + 2*c*Midbin_perBin+ 3*d*Midbin_perBin^2
  print(cub_check)
  list_recombin_rate_check[bin] = cub_check
  list_recombin_rate[bin] = CubRegDeriv_Bin/(10^(-6))  # check!!!
  
  # Copy to file
  cat(paste(chrm, bin, x1, x2, Midbin_perBin, CubRegDeriv_Bin/(10^(-6)), sep=","),
      file="local_recomb_rate_per_bin.txt",append=TRUE, sep='\n')
}

print(list_midbins)
print(list_recombin_rate)
list_recombin_rate = list_recombin_rate
plot(c(1:length(list_midbins)), list_recombin_rate, 
     xlab='Midbin', ylab='Local recombination rate',
     main=paste('Local recombi.rate: chrm', chrm),xaxt='n')
axis(side = 1, at=1:(length(list_midbins)), labels = F)
text(seq(1,length(list_midbins),by=1), par("usr")[3]-0.05, 
       srt = 60, adj=1, xpd = TRUE,
       labels = list_midbins, cex=0.50)
abline(h=0)
}
dev.off()