#!/bin/bash

#SBATCH --job-name=cov_mask
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=48:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-24%24

cd /scratch4/mschatz1/rmccoy22/rmccoy22/mask

if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
  chrom="chrX"
elif [[ ${SLURM_ARRAY_TASK_ID} -eq 24 ]]
then
  chrom="chrY"
else
  chrom="chr${SLURM_ARRAY_TASK_ID}"
fi

python coverage_mask.py ${chrom} > coverage_mask_${chrom}.txt
