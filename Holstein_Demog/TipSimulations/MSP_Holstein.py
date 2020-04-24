# 26th September 2018
# Python script for executing msprime simulations of Holstein demography,
# Calculate mean tip lengths to determine age range of mutations
# Execute with python3 MSP_Holstein.py 1000 > MSP_Results.dat
# (Results piped to STDOUT)

# Update 21st Jan 2019
# Using actual data files provided by Simon Boitard

import msprime
import math
import sys
import numpy as np
import scipy as sp

def holstein_demog(n,reps):
	# Using the data from Boitard et al 2016 PLoS Genet, we can recreate the 
	# demographic expansion of the HOL population, tracing back in time.
	# First outline each of the population sizes at each of the step changes
	N_0 = 793
	N_1 = 1077
	N_2 = 1320
	N_3 = 1815
	N_4 = 3892
	N_5 = 6909
	N_6 = 11367
	N_7 = 15603
	N_8 = 20160
	N_9 = 23783
	N_10 = 27155
	N_11 = 36513
	N_12 = 50726
	N_13 = 52882
	N_14 = 63092
	N_15 = 72842
	N_16 = 71209
	N_17 = 79079
	N_18 = 86207
	N_19 = 65485
	N_20 = 31652
	# Now for the times at which step changes occur, in number of generations	
	T_1 = 9
	T_2 = 24
	T_3 = 47
	T_4 = 83
	T_5 = 140
	T_6 = 229
	T_7 = 367
	T_8 = 584
	T_9 = 924
	T_10 = 1455
	T_11 = 2288
	T_12 = 3590
	T_13 = 5629
	T_14 = 8821
	T_15 = 13818
	T_16 = 21639
	T_17 = 33882
	T_18 = 53045
	T_19 = 83043
	T_20 = 130000
	# Setting up population size changes
	demographic_events = [
		msprime.PopulationParametersChange(time=T_1, initial_size=N_1),
		msprime.PopulationParametersChange(time=T_2, initial_size=N_2),
		msprime.PopulationParametersChange(time=T_3, initial_size=N_3),
		msprime.PopulationParametersChange(time=T_4, initial_size=N_4),
		msprime.PopulationParametersChange(time=T_5, initial_size=N_5),
		msprime.PopulationParametersChange(time=T_6, initial_size=N_6),
		msprime.PopulationParametersChange(time=T_7, initial_size=N_7),
		msprime.PopulationParametersChange(time=T_8, initial_size=N_8),
		msprime.PopulationParametersChange(time=T_9, initial_size=N_9),
		msprime.PopulationParametersChange(time=T_10, initial_size=N_10),
		msprime.PopulationParametersChange(time=T_11, initial_size=N_11),
		msprime.PopulationParametersChange(time=T_12, initial_size=N_12),
		msprime.PopulationParametersChange(time=T_13, initial_size=N_13),
		msprime.PopulationParametersChange(time=T_14, initial_size=N_14),
		msprime.PopulationParametersChange(time=T_15, initial_size=N_15),
		msprime.PopulationParametersChange(time=T_16, initial_size=N_16),
		msprime.PopulationParametersChange(time=T_17, initial_size=N_17),
		msprime.PopulationParametersChange(time=T_18, initial_size=N_18),
		msprime.PopulationParametersChange(time=T_19, initial_size=N_19),
		msprime.PopulationParametersChange(time=T_20, initial_size=N_20)
	]
	# Now to simulate samples from this demography
	replicates1 = msprime.simulate(sample_size=n, Ne=N_0, demographic_events=demographic_events, num_replicates=reps)
	idx = np.arange(reps)
	tips1 = np.zeros((reps,n))
	for j, tree_sequence in enumerate(replicates1):
		tree = tree_sequence.first()
		for i in tree_sequence.samples():
			tips1[j,i] = tree.branch_length(i)
	
	# Calculating mean over tips for each sim, then mean of means
	tipmean = np.mean(tips1,axis=1)
	simmean = np.mean(tipmean)
	# Now to calculate CIs via bootstrapping
	# First create 'nbs' (1000 default) random draws between 0 and 'reps'
	nbs = 1000
	ridx = np.random.randint(reps, size=(nbs, reps))
	tipsBS = np.zeros((nbs,reps))
	# Then use indices to resample and recalculate mean from estimates
	for k in range(0,nbs):
		for i in range(0,reps):
			tipsBS[k,i] = tipmean[ridx[k,i]]
	bsres = np.mean(tipsBS,axis=1)	# Distribution of bootstrap values
	simlb = np.quantile(bsres, 0.025)
	simub = np.quantile(bsres, 0.975)
	print("{} {} {} {}".format(val,simmean,simlb,simub))

treps=int(sys.argv[1])
values1 = np.arange(10,100,20)
values2 = np.arange(100,1050,50)
values = np.concatenate((values1,values2),axis=None)
for val in values:
	holstein_demog(val,treps)
