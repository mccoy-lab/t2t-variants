#!/bin/bash

#SBATCH --job-name=rare_haps_chm13
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=36:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-22

#module restore
#module load gcc
#module load r

cd /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/

# replace with GRCh38 if desired
SAMPLE='CHM13'

mkdir /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}

cd /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}

~/scratch4-mschatz1/rmccoy22/code/plink \
--vcf /scratch4/mschatz1/rmccoy22/1kg-GRCh38-NYGC-highcoverage/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${SLURM_ARRAY_TASK_ID}.filtered.shapeit2-duohmm-phased.vcf.gz \
--make-bed \
--keep-allele-order \
--out /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${SLURM_ARRAY_TASK_ID};

~/scratch4-mschatz1/rmccoy22/code/plink \
--bfile /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${SLURM_ARRAY_TASK_ID} \
--r2 \
--ld-window 100 \
--ld-window-r2 1.0 \
--maf 0.1 \
--out /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${SLURM_ARRAY_TASK_ID};

rm /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${SLURM_ARRAY_TASK_ID}.bed;

cd /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/

Rscript --verbose rare_haps_ref.R ${SAMPLE} ${SLURM_ARRAY_TASK_ID} >& ${SAMPLE}_${SLURM_ARRAY_TASK_ID}.log
