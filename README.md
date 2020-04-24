README FOR MILK SDS FILES

Scripts used in the study "Using singleton densities to detect recent selection in Bos taurus". Comments to m.hartfield@ed.ac.uk.

***Gamma Shapes Calculation***

'DAF_Calc' folder contains scripts for computing gamma shapes on server, and the actual files used, for both high and low N0 scenarios. Based on those provided in Yair Field's original scripts (https://github.com/yairf/SDS). Scripts assume that (i) simulpop is available (used via conda in this case), (ii) ms is already present in each directory, it can be obtained from http://home.uchicago.edu/rhudson1/source/mksamples.html.

***Tip Simualtions***

Files in 'Holstein_Demog' contain code for calculating mean tip lengths for Holstein Demography. MSPrime is used to simulate gene trees (https://msprime.readthedocs.io/en/stable/).

***Data analyses***

This folder contains the files used to analyse the raw SDS scores.

The SDS files have to be first downloaded from <DRYAD URL HERE>.

'SDS_DataAnalysis.R' performs the data analyses. It can be run for the command line using 'Rscript SDS_DataAnalysis.R Demog'. 'Demog' is a switch indicating which demography to analyse (1 = high N0, 2 = low N0).

'SigSDSAnnotate.sh' will locate genes close to high SDS regions. Requires use of the Bos taurus GTF gene annotation file (ftp://ftp.ensembl.org/pub/release-96/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.96.gtf.gz) and Bedtools (https://bedtools.readthedocs.io/en/latest/).

The 'Milk_Protein_Genes' subfolder contains data on milk protein gene regions. 'QTLStatureDat' contains info on stature QTLs with positions in the UMD assembly. 'StatureQTLs' contains information on those regions using the ARS-UCD assembly.

'Permutation_Tests' contains script for performing the permutation analyses, set up to be run on a cluster.
