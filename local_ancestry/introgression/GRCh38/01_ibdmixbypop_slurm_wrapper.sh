#!/bin/bash

#SBATCH --job-name=ibdmix
#SBATCH -N 1
#SBATCH --time=10:00:00
#SBATCH --mem=24GB
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=1-23

cd /scratch4/mschatz1/rmccoy22/syan11/ibdmix

i=${SLURM_ARRAY_TASK_ID}
if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
	i=X;
	./run_ibdmix_bypop_chrX.sh ${i} ${POPID};
else
	./run_ibdmix_bypop.sh ${i} ${POPID};
fi
