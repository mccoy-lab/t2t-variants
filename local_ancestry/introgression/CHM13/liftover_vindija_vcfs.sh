#!/bin/bash

### USAGE: ./liftover_vindija_vcfs.sh <chr ID>
### <chr ID>: name of chromosome to lift over (ex: `chr21`)
###
### Lift over Vindija Neanderthal vcfs from hg19 to hg38.

############################################################################

### DATA

# path to GATK install
GATKPATH=/scratch4/mschatz1/rmccoy22/code/gatk-4.1.9.0/

# path to VCFs to liftover
VCFPATH=/scratch4/mschatz1/rmccoy22/archaic_hominin/vcfs/vindija/
# path to liftover chain file
CHAIN=/scratch4/mschatz1/rmccoy22/rmccoy22/liftover_vcfs/GRCh37ToGRCh38.over.chain
# path to target reference genome
REFERENCE=/scratch4/mschatz1/rmccoy22/rmccoy22/ref/GRCh38.fa
# name of chromosome to lift over
ID=$1

############################################################################


# generate GRCh38 reference dictionary (needed for GATK)
samtools dict $REFERENCE > ${REFERENCE}.dict

# lift over, asking GATK to fix genotypes if REF/ALT allele is swapped between assemblies
${GATKPATH}/gatk LiftoverVcf \
  -I ${VCFPATH}/chr${ID}_mq25_mapab100.vcf.gz \
  -O ${VCFPATH}/chr${ID}_GRCh38.vcf.gz \
  --CHAIN ${CHAIN} \
  --REJECT ${VCFPATH}/${ID}_rejects.vcf.gz \
  -R $REFERENCE \
  --RECOVER_SWAPPED_REF_ALT true \
  --MAX_RECORDS_IN_RAM 10000
  &>> ${ID}_liftover.log
