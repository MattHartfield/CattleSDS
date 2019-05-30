## Examine and parse the estimated local recombination rates from Belen's script.
## TB Oct 7, 2018

library(tidyverse) 
## Data 
# chr: chromosome number
# bin : bin number 
# begin end and mid are the physical postions of the bin a chromosome
# cMperMb: the estimated local recombination rate estimated for the bin

df=read_csv("local_recomb_rate_per_bin_NEW.txt",col_names = c("chr","bin","begin","end","mid", "cMperMb"),progress = T)
names(df)
glimpse(df)
head(df)
dim(df)
summary(df)

df %>%  filter(!is.na(cMperMb)) -> df # excluding bins with no LRR estimated

# an overview of the genomewide distribution of the local recombination rates
ggplot(data = df) + geom_histogram(aes(x=cMperMb),binwidth = 0.1)


df %>% group_by(chr) %>% 
  summarise(mean=mean(cMperMb),
            median=median(cMperMb),
            nBins=n()) ->dfSummaryByChr

dfSummaryByChr
dim(dfSummaryByChr)

# Per chromsome
ggplot(data = df) + geom_histogram(aes(x=cMperMb, fill=chr)) + facet_wrap(~chr, ncol=2)
