#!/bin/bash

#SBATCH --job-name=liftoverVindija
#SBATCH -N 1
#SBATCH --time=5:00:00
#SBATCH --mem=24GB
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-23


cd /scratch4/mschatz1/rmccoy22/rmccoy22/liftover_vcfs

i=${SLURM_ARRAY_TASK_ID}
if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]] ; then
	i=X
fi

./liftover_1KGP_vcfs.sh ${i}
