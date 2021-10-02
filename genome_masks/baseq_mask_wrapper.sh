#!/bin/bash

#SBATCH --job-name=baseq_mask
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=12:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-22%22

cd /scratch4/mschatz1/rmccoy22/rmccoy22/mask

if [[ ${chrom} -eq 23 ]]
then
  chrom="chrX"
elif [[ ${chrom} -eq 24 ]]
then
  chrom="chrY"
else
  chrom="chr${SLURM_ARRAY_TASK_ID}"
fi

python baseq_mask.py ${chrom} > baseq_mask_${chrom}.txt
