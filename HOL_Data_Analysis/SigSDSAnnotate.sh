#!/bin/bash
# 29th March 2019
# Script to look for genes close to significant SDS scores
for fname in high low
do
	rm OutTables/SigSDSGenesProteinCoding_${fname}N0.csv
	bedtools intersect -a sigsds_${fname}N0.bed -b Bos_taurus.ARS-UCD1.2.96.gtf -wb > SigSDSAnnotate_${fname}N0.dat
	awk '$6=="gene"' SigSDSAnnotate_${fname}N0.dat > SigSDSGenes_${fname}N0.dat
	awk 'BEGIN{print "Chromosome,GeneName,StartPosition,EndPosition"};{print $1,$17,$7,$8};' SigSDSGenes_${fname}N0.dat | sed -e 's/ /,/g' | sed -e "s/;//g" | sed -e "s/GeneName/Gene Name/" | sed -e "s/StartPosition/Start Position/" | sed -e "s/EndPosition/End Position/" | sed -e "s/ensembl/(Unnamed)/" > SigSDSGenes_${fname}N0P1.dat
	awk -F ";" '{print $(NF-1)};' SigSDSGenes_${fname}N0.dat | awk '{print $2}' | sed -e "s/_/ /" | sed -e "s/\\\"//g" | awk 'BEGIN{print "Gene Biotype"};{print $0};' | awk '{ for(i=1;i<=NF;i++) { j=toupper(substr($i,1,1)); printf "%s%s ",j,substr($i,2)} print ""}' > SigSDSGenes_${fname}N0P2.dat
	paste -d ,, SigSDSGenes_${fname}N0P1.dat SigSDSGenes_${fname}N0P2.dat > OutTables/SigSDSGenesProteinCoding_${fname}hN0.csv
	rm SigSDSAnnotate_${fname}N0.dat SigSDSGenes_${fname}N0*.dat
done
