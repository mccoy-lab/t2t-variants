#!/bin/bash

#SBATCH --job-name=rare_haps
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=36:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-100%25

# note that modules must be loaded in environment from which you sbatch

cd /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/

SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" 1KGP_samples.txt`

mkdir /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}

cd /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}

cat ../1KGP_samples.txt | grep -v -w ${SAMPLE} | awk '{print $1"\t"$1}' > ${SAMPLE}.keep

for i in {1..22}
do
~/scratch4-mschatz1/rmccoy22/code/plink \
--vcf /scratch4/mschatz1/rmccoy22/code/htslib-1.11/tabix /scratch4/mschatz1/rmccoy22/1kg-GRCh38-NYGC-highcoverage/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz  \
--make-bed \
--keep-allele-order \
--keep ${SAMPLE}.keep \
--out /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${i};
~/scratch4-mschatz1/rmccoy22/code/plink \
--bfile /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${i} \
--r2 \
--ld-window 100 \
--ld-window-r2 1.0 \
--maf 0.1 \
--out /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${i};
rm /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${i}.bed;
done

cd /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/

Rscript --verbose rare_haps.R ${SAMPLE} >& ${SAMPLE}.log
