#!/bin/bash

#SBATCH --job-name=prepVCFs
#SBATCH -N 1
#SBATCH --time=10:00:00
#SBATCH --mem=20G
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=22

cd /scratch4/mschatz1/rmccoy22/syan11/gwas_ld/

############################################################################

### DATA

HTSLIB_PATH=/scratch4/mschatz1/rmccoy22/code/htslib-1.11
BCFTOOLS=/scratch4/mschatz1/rmccoy22/code/bcftools-1.11/bcftools

#VCFPATH=/scratch4/mschatz1/rmccoy22/syan11/gwas_ld/nonsyn_winnowmap_100mer_vcfs
#VCFPATH=/scratch4/mschatz1/rmccoy22/syan11/gwas_ld/dbsnp_lifted
#VCFPATH=/scratch4/mschatz1/rmccoy22/1kg-CHM13-recalibrated-PASS
#VCFPATH=/scratch4/mschatz1/rmccoy22/syan11/gwas_ld/nonsyn_winnowmap_100mer_vcfs
VCFPATH=/scratch4/mschatz1/rmccoy22/syan11/gwas_ld/nonsyn_vcfs

############################################################################


# get chromosome number
i=${SLURM_ARRAY_TASK_ID}
if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]; then
	i=X
fi
# check that file is there
ls $VCFPATH/*chr${i}.*vcf.gz

# annotate 1KGP VCFs with variant IDs for novel regions
#mkdir -p 1kg-CHM13-recalibrated-PASS-annot
#$BCFTOOLS annotate \
#	-O z \
#	-o nonsyn_winnowmap_100mer_vcfs/chr${i}.nonsyn_winnowmap_100mer.annot.vcf.gz \
#	--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
#	$VCFPATH/*chr${i}.*vcf.gz
#${HTSLIB_PATH}/tabix nonsyn_winnowmap_100mer_vcfs/chr${i}.nonsyn_winnowmap_100mer.annot.vcf.gz

# concatenate variant annotated nonsyntenic vcfs from running Ohana
INPUT_VCFS=/scratch4/mschatz1/rmccoy22/syan11/snp_selection/snp_genotypes/chr${i}/for_ohana
# make list of chromosome vcfs to concatenate
ls ${INPUT_VCFS}/chr${i}*vcf.gz > ${VCFPATH}/chr${i}_files.txt
$BCFTOOLS concat \
	-O z \
	-o ${VCFPATH}/chr${i}.1Mb_nonsyn.annot.vcf.gz \
	-f ${VCFPATH}/chr${i}_files.txt \
	--threads 48
${HTSLIB_PATH}/tabix ${VCFPATH}/chr${i}.1Mb_nonsyn.annot.vcf.gz
