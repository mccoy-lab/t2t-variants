#!/bin/bash

### USAGE: ./ibdmix.sh <chr ID>
### <chr ID>: number of chromosome to run on (ex: `21`)
###
### Runs ibdmix on CHM13, using the CEU superpopulation as a reference.

############################################################################

### DATA

# path to bcftools install
BCFTOOLS=/scratch4/mschatz1/rmccoy22/code/bcftools-1.11/bin/bcftools
# path to tabix install
TABIX=/scratch4/mschatz1/rmccoy22/code/htslib-1.11/tabix
# path to IBDmix's `src` directory
IBDMIX_PATH=/scratch4/mschatz1/rmccoy22/rmccoy22/IBDmix/build/src

# path to 1KGP NYGC VCFs (on GRCh38)
NYGC_VCFS=/scratch4/mschatz1/rmccoy22/1kg-GRCh38-NYGC-highcoverage
# path to CHM13 VCFs (on GRCh38)
CHM13_VCFS=/scratch4/mschatz1/rmccoy22/chm13_local_ancestry/liftover_chm13_to_hg38
# path to Vindija Neanderthal VCFs, lifted over to GRCh38
VINDIJA_VCFS=/scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/vindija

############################################################################


### prep CHM13 vcfs

# subset 1KGP NYGC vcfs to just CEU samples
$BCFTOOLS view \
-S /scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/1kg-GRCh38-NYGC-highcoverage/eur.txt \
-v snps \
-O z \
-o /scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/1kg-GRCh38-NYGC-highcoverage/chr${1}_1kgp_nygc_grch38.vcf.gz \
${NYGC_VCFS}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${1}.filtered.shapeit2-duohmm-phased.vcf.gz

$TABIX \
-p vcf \
/scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/1kg-GRCh38-NYGC-highcoverage/chr${1}_1kgp_nygc_grch38.vcf.gz

# subset CHM13 vcf to just snps
$BCFTOOLS view \
-v snps \
-O z \
-o /scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/chm13/chr${1}_chm13_grch38.vcf.gz \
${CHM13_VCFS}/chr${1}_CHM13_joint_sorted.vcf.gz

$TABIX \
-p vcf \
/scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/chm13/chr${1}_chm13_grch38.vcf.gz

# merge query vcfs (CHM13 and 1KGP CEU samples)
$BCFTOOLS merge \
-0 \
-O z \
-o /scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/chm13_1kgp_merged/chr${1}_chm13_1kgp_grch38_merged.vcf.gz \
/scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/chm13/chr${1}_chm13_grch38.vcf.gz \
/scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/1kg-GRCh38-NYGC-highcoverage/chr${1}_1kgp_nygc_grch38.vcf.gz
# remove `chr` prefix
zcat /scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/chm13_1kgp_merged/chr${1}_chm13_1kgp_grch38_merged.vcf.gz \
| sed 's/chr//g' > \
/scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/chm13_1kgp_merged/chr${1}_chm13_1kgp_grch38_merged.vcf

# unzip vindija vcf
gunzip ${VINDIJA_VCFS}/chr${1}_GRCh38.vcf.gz


### run ibdmix

${IBDMIX_PATH}/generate_gt \
-a ${VINDIJA_VCFS}/chr${1}_GRCh38.vcf \
-m /scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/chm13_1kgp_merged/chr${1}_chm13_1kgp_grch38_merged.vcf \
-o /scratch4/mschatz1/rmccoy22/archaic_hominin/gt/chr${1}_merged_GRCh38.gt

${IBDMIX_PATH}/ibdmix \
-g /scratch4/mschatz1/rmccoy22/archaic_hominin/gt/chr${1}_merged_GRCh38.gt \
-o /scratch4/mschatz1/rmccoy22/archaic_hominin/ibdmix/chr${1}_vindija_GRCh38.out