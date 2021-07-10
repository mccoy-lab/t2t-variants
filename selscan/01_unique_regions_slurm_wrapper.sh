#!/bin/bash

### USAGE: sbatch 01_unique_regions_slurm_wrapper.sh

#SBATCH --job-name=uniqueRegions
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --mem=117600MB
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=1-21,23

cd /scratch4/mschatz1/rmccoy22/syan11/snp_selection

# read in chromosome number
i=${SLURM_ARRAY_TASK_ID}
if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
        i=X;
fi

echo "Running ./get_unique_regions.sh on chr"$i
./get_unique_regions.sh $i
