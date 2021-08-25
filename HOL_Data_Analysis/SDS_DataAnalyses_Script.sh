#!/bin/sh
# 21st Apr 2020
# Running analyses for low, high N0
rm -r Rout
mkdir Rout
Rscript SDS_DataAnalysis.R 1 &> Rout/HighN0.out &
Rscript SDS_DataAnalysis.R 2 &> Rout/LowN0.out &
wait
echo 'Both N0 analyses complete'
