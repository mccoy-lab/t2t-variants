#!/bin/bash

#SBATCH --job-name=XchrPhase3
#SBATCH -N 1
#SBATCH --time=10:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=950GB
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=1-23

cd /scratch4/mschatz1/rmccoy22/chm13_local_ancestry

i=${SLURM_ARRAY_TASK_ID}
if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]; then
        i=X
	./run_rfmix_Xchr.sh ${i}
else
	./run_rfmix_phase3.sh ${i}
fi
