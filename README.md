# README FOR CATTLE SDS FILES

Scripts used in the study "Using singleton densities to detect recent selection in Bos taurus". Comments to m.hartfield@ed.ac.uk.

### Gamma Shapes Calculation

`DAF_Calc` folder contains scripts for computing gamma shapes on a computer cluster, along with the actual files used for both high and low N0 scenarios. Based on those provided in the [original SDS scripts](https://github.com/yairf/SDS "SDS Scripts"). Scripts assume that (i) [simuPOP](http://simupop.sourceforge.net/ "simuPOP") is available, (ii) ms is present in each directory; it can be obtained from [here](http://home.uchicago.edu/rhudson1/source/mksamples.html "ms").

### Tip Simulations

Files in `Holstein_Demog` contain code for calculating mean tip lengths for Holstein Demography. [MSPrime](https://msprime.readthedocs.io/en/stable/ "MSPrime") is used to simulate gene trees.

### Data analyses

This folder contains the files used to analyse the raw SDS scores. The SDS data files (both SDS scores and allele polarisation information) have been archived on Dryad and are currently under curation (provisional DOI [https://doi.org/10.5061/dryad.547d7wm8q](https://doi.org/10.5061/dryad.547d7wm8q)).

`SDS_DataAnalysis.R` performs the data analyses. It can be run for the command line using `Rscript SDS_DataAnalysis.R Demog`. `Demog` is a switch indicating which demography to analyse (1 = high N0, 2 = low N0). You first have to change the `setwd` command on line 180 to point to the directory on your computer where the data files are stored. You also need to make `OutFigures` and `OutTables` folders in that directory to store outputs. You can also run the `SDS_DataAnalyses_Script.sh` script to execute both high- and low-N0 analyses simultaneously.

`SigSDSAnnotate.sh` will locate genes close to high SDS regions. Requires use of the [Bos taurus GTF gene annotation file](https://www.ensembl.org/Bos_taurus/Info/Index "B. taurus GTF") and [Bedtools](https://bedtools.readthedocs.io/en/latest/ "Bedtools").

The `MilkQTL` subfolder contains data on QTLs on milk protein content and fat content. `QTLStatureDat` contains info on stature QTLs with positions in the UMD assembly. `StatureQTLs` contains information on those regions using the ARS-UCD assembly.

`Permutation_Tests` contains script for performing the permutation analyses, set up to be run on a cluster.

### Power Simulations

`backTrajectory_printout.py` is the Python script used to simulate allele trajectories, assuming Holstein demography. It is executed as `python backTrajectory_printout.py Sel Rep`, for `Sel` the homozygote selection coefficient and `Rep` the simulation replicate. It uses simuPOP to simulate trajectories. Results can then be used by [mbs](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-166 "mbs").
