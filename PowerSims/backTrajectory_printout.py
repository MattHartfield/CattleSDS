#!/usr/bin/env python

# Script to simulate partial sweep trajectory using simuPOP and HOL demography
# Based on script from http://bopeng.github.io/simuPOP/userGuide_ch7_sec2.html#backward-time-trajectory-simulations-function-simulatebackwardtrajectory

import sys
import simuPOP as sim
from simuPOP.utils import Trajectory, simulateBackwardTrajectory

num_of_genereations = 130000

# Input command line arguments
sel = float(sys.argv[1])	# Homozygote selection coefficient
rep = int(sys.argv[2])		# Replicate to run

def genBackwards(gen_forward):

    return int(abs(gen_forward-(num_of_genereations+1)))
    
# Loop of Holstein Demography (from 'backward' file)
def baseNt_Holstein(gen):

    my_gen_backwards = genBackwards(gen)
    if my_gen_backwards >= 130000:
        return 31652
    if my_gen_backwards >= 83043:
        return 65485
    if my_gen_backwards >= 53045:
        return 86207
    if my_gen_backwards >= 33882:
        return 79079
    if my_gen_backwards >= 21639:
        return 71209
    if my_gen_backwards >= 13818:
        return 72842
    if my_gen_backwards >= 8821:
        return 63092
    if my_gen_backwards >= 5629:
        return 52882
    if my_gen_backwards >= 3590:
        return 50726
    if my_gen_backwards >= 2288:
        return 36513
    if my_gen_backwards >= 1455:
        return 27155
    if my_gen_backwards >= 924:
        return 23783
    if my_gen_backwards >= 584:
        return 20160
    if my_gen_backwards >= 367:
        return 15603
    if my_gen_backwards >= 229:
        return 11367
    if my_gen_backwards >= 140:
        return 6909
    if my_gen_backwards >= 83:
        return 3892
    if my_gen_backwards >= 47:
        return 1815
    if my_gen_backwards >= 24:
        return 1320
    if my_gen_backwards >= 9:
        return 1077
    if my_gen_backwards >= 1:
        return 793

def fitness(sp):
    return [1, 1+0.5*sp, 1+sp]

# simulate a trajectory backward in time
traj = simulateBackwardTrajectory(N=baseNt_Holstein, fitness=fitness(sel), nLoci=1, endGen=num_of_genereations, endFreq=[0.7])
nzero=baseNt_Holstein(num_of_genereations)
	
# Printing out trajectory file
# NOTE: backward file indexes generation at 1, but I rescale to start from time zero so can be used in coalescent simulation
outfile=open('traj/traj_HOL_' + str(int(rep)) + '_0.dat','w')
for i in range(num_of_genereations,0,-1):
	my_gen = str( (genBackwards(i)-1)/(4*nzero) )			# Scaling time by 4*nzero
	my_gen2 = str( (genBackwards(i-1)-1)/(4*nzero) )
	my_popsize = str(baseNt_Holstein(i)/(nzero))			# Scaling pop by 1*nzero
	my_freq = str(traj.freq(i,0)[0])
	outfile.write( my_gen + "\t" + my_gen2 + "\t" + my_popsize + "\t" + my_freq  +"\n")
	
# Last line going to infinity ('999' in mbs notation)
outfile.write( my_gen2 + "\t 999 \t" + my_popsize + "\t" + my_freq  +"\n")

# EOF