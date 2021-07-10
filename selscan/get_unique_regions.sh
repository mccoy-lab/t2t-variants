#!/bin/bash

### USAGE: ./get_unique_regions.sh <chr>
### <chr>: number of chromosome to run on (ex: `21`)
###
### Gets regions of VCF that are unique to CHM13.

############################################################################

### DATA

# number of chromosome to run on
i=$1

# path to bcftools install
BCFTOOLS_PATH=/scratch4/mschatz1/rmccoy22/code/bcftools-1.11
# path to htslib install
HTSLIB_PATH=/scratch4/mschatz1/rmccoy22/code/htslib-1.11
# path to plink install
PLINK=/scratch4/mschatz1/rmccoy22/code/plink

# file of regions of no synteny between CHM13 and GRCh38
UNIQUE1=chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed

# path to directory with snp vcfs
INPUT_VCFs=/scratch4/mschatz1/rmccoy22/1kg-CHM13-nonrecalibrated
# file with 2,504 core 1000 genomes samples, in same order as Ohana Q matrix
SAMPLES=samples.txt

# path to output modified snp vcfs to
mkdir -p snp_genotypes/chr$i
VCFPATH=snp_genotypes/chr$i

############################################################################


### subset vcf to just CHM13 "unique" regions

# subset to 2504 core 1000 Genomes samples
${BCFTOOLS_PATH}/bcftools view \
	-o ${VCFPATH}/chr${i}_2504samples.vcf.gz \
	-O z \
	-S $SAMPLES \
	--force-samples \
	${INPUT_VCFs}/chr${i}.*.vcf.gz
echo "vcf subsetted"
# index subsetted vcf
${HTSLIB_PATH}/tabix -p vcf ${VCFPATH}/chr${i}_2504samples.vcf.gz
echo "vcf indexed"

# intersect SNP VCF with no synteny CHM13 regions
${HTSLIB_PATH}/tabix -h \
	-R $UNIQUE1 \
	${VCFPATH}/chr${i}_2504samples.vcf.gz \
	> ${VCFPATH}/chr${i}_unique.vcf
echo "intersected with CHM13 unique regions"
# bgzip and index intersected vcf
${HTSLIB_PATH}/bgzip ${VCFPATH}/chr${i}_unique.vcf
${HTSLIB_PATH}/tabix -p vcf ${VCFPATH}/chr${i}_unique.vcf.gz
echo "finished zipping and indexing intersected vcf"


### split vcf into portions of size 3GB for parallelizing

# get VCF header
${BCFTOOLS_PATH}/bcftools view -h \
	${VCFPATH}/chr${i}_unique.vcf.gz \
	> ${VCFPATH}/chr${i}_unique_header.txt
# get variants
zcat ${VCFPATH}/chr${i}_unique.vcf.gz | \
	grep -v "#" \
	> ${VCFPATH}/chr${i}_unique_variants.txt
# split up variants into files of size 1GB (about ~100,000 lines)
split -l 100000 \
        -d \
        ${VCFPATH}/chr${i}_unique_variants.txt \
        ${VCFPATH}/chr${i}_unique_split
# convert split variant files into VCFs
for splitfile in ${VCFPATH}/chr${i}_unique_split*
do
	# combine header and variants
	cat ${VCFPATH}/chr${i}_unique_header.txt $splitfile \
		> ${splitfile}.vcf
	${HTSLIB_PATH}/bgzip ${splitfile}.vcf
	rm $splitfile
done
rm ${VCFPATH}/chr${i}_unique_header.txt ${VCFPATH}/chr${i}_unique_variants.txt
