#!/bin/bash
#SBATCH --job-name=prepVCFs
#SBATCH -N 1
#SBATCH --time=0:30:00
#SBATCH --mem=117600MB
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu

### USAGE: sbatch prep_split_vcfs.sh <chr_num> <split_num>
### <chr_num>: number of chromosome to run on (ex: `21`)
### <split_num>: number of split VCF for that chromosome
###
### Adds variant IDs to split VCFs and converts into plink ped files for Ohana.

############################################################################

### DATA

# chromosome to run on
i=$1
# number of split VCF
num=$2

# path to bcftools install
BCFTOOLS_PATH=/scratch4/mschatz1/rmccoy22/code/bcftools-1.11
# path to plink install
PLINK=/scratch4/mschatz1/rmccoy22/code/plink

# path to directory where vcfs are
VCFPATH=snp_genotypes/chr$i
# `within` file of sample IDs and populations for plink
WITHIN_FILE=1KGP_within.txt

############################################################################


cd /scratch4/mschatz1/rmccoy22/syan11/snp_selection
echo "running ./prep_split_vcfs.sh on chr"$i "split" $num

# get path to split vcf
input_vcf=${VCFPATH}/chr${i}_unique_split${num}.vcf.gz
# set name of split vcf with variant ID annotations
id_vcf=${VCFPATH}/chr${i}_unique_IDs_split${num}.vcf.gz

# add variant IDs to VCF
${BCFTOOLS_PATH}/bcftools annotate \
	-O z \
	-o $id_vcf \
	--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
	$input_vcf
echo "added variant IDs for chr"$i "split" $num

# calculate population-specific allele frequencies
${PLINK} --vcf $id_vcf \
        --within $WITHIN_FILE \
        --freq gz \
        --keep-allele-order \
        --out $id_vcf
echo "calculated allele frequencies for chr"$i "split" $num

# convert VCF into .ped file
${PLINK} --vcf $id_vcf \
	--out $id_vcf \
	--recode 12 \
	--tab
echo "vcf converted into ped for chr"$i "split" $num