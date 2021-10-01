#!/bin/bash

#SBATCH --job-name=mapq_mask
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=12:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-22%22

cd /scratch4/mschatz1/rmccoy22/rmccoy22/mask

chrom="chr${SLURM_ARRAY_TASK_ID}"

python mapq_mask.py ${chrom} > mapq_mask_${chrom}.txt
