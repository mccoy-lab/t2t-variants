#!/bin/bash

#SBATCH --job-name=rare_haps_batch
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=36:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-990%44

#module restore
#module load gcc
#module load r

cd /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/

SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" sbatch_list_2.txt | awk '{print $1}'`
CHROM=`sed "${SLURM_ARRAY_TASK_ID}q;d" sbatch_list_2.txt | awk '{print $2}'`

mkdir -p /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}
mkdir -p /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/output

cd /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}

if [ ! -f ${SAMPLE}.keep ]; then
    cat ../1KGP_samples.txt | grep -v -w ${SAMPLE} | awk '{print $1"\t"$1}' > ${SAMPLE}.keep
fi

~/scratch4-mschatz1/rmccoy22/code/plink \
--vcf /scratch4/mschatz1/rmccoy22/1kg-GRCh38-NYGC-highcoverage/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHROM}.filtered.shapeit2-duohmm-phased.vcf.gz \
--make-bed \
--keep-allele-order \
--snps-only \
--biallelic-only \
--keep ${SAMPLE}.keep \
--out /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${CHROM}

~/scratch4-mschatz1/rmccoy22/code/plink \
--bfile /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${CHROM} \
--r2 \
--ld-window 100 \
--ld-window-r2 1.0 \
--maf 0.1 \
--out /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${CHROM}

rm /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/${SAMPLE}/chr${CHROM}.bed

cd /scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/

Rscript --verbose rare_haps.R ${SAMPLE} ${CHROM} >& ${SAMPLE}_${CHROM}.log
