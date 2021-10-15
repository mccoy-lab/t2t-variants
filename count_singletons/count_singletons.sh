#!/bin/bash

#SBATCH --job-name=count_singletons
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=36:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-22%22

#module restore
#module load gcc
#module load r

module load vcftools

cd /scratch4/mschatz1/rmccoy22/rmccoy22/count_singletons/

vcftools \
  --keep unrelated_samples.keep \
  --singletons \
  --gzvcf /scratch4/mschatz1/rmccoy22/1kg-CHM13-recalibrated-PASS/1kgp.chr${SLURM_ARRAY_TASK_ID}.recalibrated.snp_indel.pass.vcf.gz \
  --out chr${SLURM_ARRAY_TASK_ID}

PATH=${PATH}:"~/scratch4-mschatz1/rmccoy22/code/bcftools-1.11/"

bcftools query \
  --samples-file unrelated_samples.keep \
  -i'AC==AN' \
  -f'%CHROM %POS %AC %AN\n' \
  /scratch4/mschatz1/rmccoy22/1kg-CHM13-recalibrated-PASS/1kgp.chr${SLURM_ARRAY_TASK_ID}.recalibrated.snp_indel.pass.vcf.gz \
  > chr${SLURM_ARRAY_TASK_ID}_CHM13_singletons.txt 
  
