#!/bin/bash

# 5th Nov 2018
# Script to (1) Simulate HOL tip-lengths using MSPRIME
# (2) Plot results in R
# Using msprime in conda

# Update 27th Feb
# Simulating and plotting both high and low N0 simulations

source activate msprime-env
rm MSP_Results.dat MSP_Results_Low.dat
python MSP_Holstein.py 1000 > MSP_Results.dat
python MSP_Holstein_LowMod.py 1000 > MSP_Results_Low.dat
Rscript MSP_DatPlot_Both.R
