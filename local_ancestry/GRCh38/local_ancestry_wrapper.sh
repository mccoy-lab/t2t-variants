#!/bin/bash
#SBATCH --job-name=lai
#SBATCH -n 12
#SBATCH --time=10:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rajiv.mccoy@jhu.edu
#SBATCH --array=1-460%50

module restore
module load r

cd /scratch4/mschatz1/rmccoy22/rmccoy22/grch38_local_ancestry

AGP=`sed "${SLURM_ARRAY_TASK_ID}q;d" agp_file_list.txt`

Rscript --verbose test_local_ancestry.R ${AGP}
