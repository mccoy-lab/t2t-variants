#!/bin/bash

### USAGE: ./run_ibdmix_bypop_chrX.sh <chr num> <population ID>
### <chr num>: number of chromosome to run on (ex: `21`)
### <population ID>: superpopulation of samples to run on (ex: `EUR`)
###
### Runs ibdmix on GRCh38 `gt` files.

############################################################################

### DATA

# path to bcftools install
BCFTOOLS=/scratch4/mschatz1/rmccoy22/code/bcftools-1.11/bin/bcftools
# path to tabix install
TABIX=/scratch4/mschatz1/rmccoy22/code/htslib-1.11/tabix
# path to IBDmix's `src` directory
IBDMIX_PATH=/scratch4/mschatz1/rmccoy22/rmccoy22/IBDmix/build/src

# path to 1KGP metadata file with sample, superpop, and sex info
METADATA=1KGP_metadata.txt
# path to 1KGP NYGC VCFs (on GRCh38)
NYGC_VCFS=/scratch4/mschatz1/rmccoy22/1kg-GRCh38-NYGC-highcoverage
# path to GRCh38 VCFs (on GRCh38)
GRCH38_VCFS=/scratch4/mschatz1/rmccoy22/chm13_local_ancestry/make_grch38_vcf
# path to Vindija Neanderthal VCFs, lifted over to GRCh38 and gunzipped
VINDIJA_VCFS=/scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/vindija

# path to output directory for `gt` files
GT_PATH=/scratch4/mschatz1/rmccoy22/syan11/ibdmix/grch38_gts
# path to output directory for ibdmix
OUTPATH=/scratch4/mschatz1/rmccoy22/syan11/ibdmix/ibdmix_output

############################################################################


i=X
# read in superpopulation ID
POP=$2
echo "chr"$i $POP


### prepping 1KGP NYGC VCFs

echo "subsetting samples and SNPs from 1KGP VCF..."
# if list of samples from superpop doesn't already exist, then make it
mkdir -p vcfs
if ! [[ -f ${POP}_samples_female.txt  ]]; then
        grep $POP $METADATA | grep "female" | cut -f 1 > vcfs/${POP}_samples_female.txt
fi
# subset samples from 1KGP NYGC VCF
mkdir -p vcfs/1KGP
$BCFTOOLS view \
	-S vcfs/${POP}_samples_female.txt \
	--force-samples \
	-v snps \
	-O z \
	-o vcfs/1KGP/chr${i}_NYGC_${POP}.vcf.gz \
	${NYGC_VCFS}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.eagle2-phased.reheadered.vcf.gz
$TABIX -p vcf \
	vcfs/1KGP/chr${i}_NYGC_${POP}.vcf.gz

### prepping GRCh38 VCFs

echo "subsetting SNPs from GRCh38 VCF..."
mkdir -p vcfs/GRCh38
# subset GRCh38 VCF to just SNPs
$BCFTOOLS view \
	-v snps \
	-O z \
	-o vcfs/GRCh38/chr${i}_${POP}_GRCh38_SNPs.vcf.gz \
	${GRCH38_VCFS}/chr${i}_GRCh38.vcf.gz
$TABIX -p vcf \
	vcfs/GRCh38/chr${i}_${POP}_GRCh38_SNPs.vcf.gz

### merging query VCFs

echo "merging GRCh38 and 1KGP VCFs..."
mkdir -p vcfs/GRCh38_1KGP_merged
# merge GRCh38 and 1KGP VCFs
$BCFTOOLS merge \
	-0 \
	-O z \
	-o vcfs/GRCh38_1KGP_merged/chr${i}_${POP}_GRCh38_NYGC_merged.vcf.gz \
	vcfs/GRCh38/chr${i}_${POP}_GRCh38_SNPs.vcf.gz \
	vcfs/1KGP/chr${i}_NYGC_${POP}.vcf.gz
# remove `chr` prefix
zcat vcfs/GRCh38_1KGP_merged/chr${i}_${POP}_GRCh38_NYGC_merged.vcf.gz | \
	sed 's/chr//g' \
	> vcfs/GRCh38_1KGP_merged/chr${i}_${POP}_GRCh38_NYGC_merged_noChr.vcf


### running ibdmix

echo "generating gt files for ibdmix..."
# generate `gt` files for ibdmix
${IBDMIX_PATH}/generate_gt \
	-a ${VINDIJA_VCFS}/chr${i}_GRCh38.vcf \
	-m vcfs/GRCh38_1KGP_merged/chr${i}_${POP}_GRCh38_NYGC_merged_noChr.vcf \
	-o ${GT_PATH}/chr${i}_${POP}_allGTs_GRCh38gts.gt

echo "running ibdmix..."
# run ibdmix
${IBDMIX_PATH}/ibdmix \
    -g ${GT_PATH}/chr${i}_${POP}_allGTs_GRCh38gts.gt \
    -o ${OUTPATH}/chr${i}_${POP}_GRCh38_ibdmix.out
