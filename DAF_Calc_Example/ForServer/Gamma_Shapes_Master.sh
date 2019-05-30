#!/bin/bash
#SBATCH -c 1
#SBATCH --mem-per-cpu 128
#SBATCH --partition short
#SBATCH --time=0:20:00

# 26th Feb 2019
# Script for parallelising the creation of 'Gamma-shapes' file
# for SDS analysis

DAFS=($(seq 0.05 0.01 0.95))
NDAF=${#DAFS[@]}

rm -r GSF* *.out
for ((i=0; i<NDAF; i++));
do
	mkdir GSF${DAFS[${i}]}
	cp Makefile ms backward_wHolstein.py GSF${DAFS[${i}]}/
	sed -i -e "s/YFI/${DAFS[${i}]}/" GSF${DAFS[${i}]}/Makefile
	echo '#!/bin/bash' > GS${i}.sh
	echo '#SBATCH -c 1' >> GS${i}.sh
	echo '#SBATCH --mem-per-cpu 256' >> GS${i}.sh
	echo '#SBATCH --partition normal' >> GS${i}.sh
	echo '#SBATCH --time=20:00:00' >> GS${i}.sh
	echo . /home/mhart/miniconda3/etc/profile.d/conda.sh >> GS${i}.sh
	echo 'conda activate simupop-env' >> GS${i}.sh
	echo 'make sim &' >> GS${i}.sh
	echo 'conda deactivate' >> GS${i}.sh
	echo 'wait' >> GS${i}.sh	
	echo 'exit 0' >> GS${i}.sh
	chmod u+x GS${i}.sh
	cp GS${i}.sh GSF${DAFS[${i}]}/
 	rm GS${i}.sh
 	cd GSF${DAFS[${i}]}/ 	
 	sbatch GS${i}.sh
 	cd ../
done

wait
exit 0
