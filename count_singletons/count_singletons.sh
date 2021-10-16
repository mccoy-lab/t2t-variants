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
PATH=${PATH}:"~/scratch4-mschatz1/rmccoy22/code/bcftools-1.11/"

cd /scratch4/mschatz1/rmccoy22/rmccoy22/count_singletons/

bcftools view \
  --threads 48 \
  -O z \
  -o /scratch4/mschatz1/rmccoy22/1kg-CHM13-recalibrated-PASS/unrelated/1kgp.chr${SLURM_ARRAY_TASK_ID}.recalibrated.snp_indel.pass.vcf.gz \
  --samples-file  /scratch4/mschatz1/rmccoy22/rmccoy22/count_singletons/unrelated_samples.keep \
  /scratch4/mschatz1/rmccoy22/1kg-CHM13-recalibrated-PASS/1kgp.chr${SLURM_ARRAY_TASK_ID}.recalibrated.snp_indel.pass.vcf.gz 

vcftools \
  --singletons \
  --gzvcf /scratch4/mschatz1/rmccoy22/1kg-CHM13-recalibrated-PASS/unrelated/1kgp.chr${SLURM_ARRAY_TASK_ID}.recalibrated.snp_indel.pass.vcf.gz \
  --out chr${SLURM_ARRAY_TASK_ID}

bcftools query \
  -i'AC==AN' \
  -f'%CHROM %POS %AC %AN\n' \
  /scratch4/mschatz1/rmccoy22/1kg-CHM13-recalibrated-PASS/unrelated/1kgp.chr${SLURM_ARRAY_TASK_ID}.recalibrated.snp_indel.pass.vcf.gz \
  > chr${SLURM_ARRAY_TASK_ID}_CHM13_singletons.txt 
