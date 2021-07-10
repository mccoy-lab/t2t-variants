#!/bin/bash

### USAGE: ./run_rfmix.sh <chr number>
### <chr number>: number of chromosome to run on (ex: `1`)
###
### Runs RFMix to infer local ancestry for a chromosome of CHM13.

############################################################################

### DATA

# path to rfmix install
RFMIX=/scratch4/mschatz1/rmccoy22/code/rfmix-master/rfmix
# path to directory for output files
outpath=rfmix_output/phase3

# path to 1KGP sample map file
sample_map=/scratch4/mschatz1/rmccoy22/1KGP_population_membership/1KGP_2504_sample_superpop_map.txt
# path to 1KGP genetic map
genetic_map=/scratch4/mschatz1/rmccoy22/genetic_maps/RFmix.whole_genome.GRCh38.map

############################################################################


# take chromosome number and CRF spacing as arguments
i=$1

# path to query VCF (CHM13)
query_vcf=/scratch4/mschatz1/rmccoy22/chm13_local_ancestry/liftover_chm13_to_hg38/chr${i}_CHM13_joint_sorted.vcf.gz
# path to VCF of reference haplotypes (1KGP phase3)
ref_vcf=/scratch4/mschatz1/rmccoy22/1kg-GRCh38-phase3/ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz

ls $query_vcf
ls $ref_vcf
ls $genetic_map
ls $sample_map

$RFMIX -f $query_vcf \
	-r $ref_vcf \
	-m $sample_map \
	-g $genetic_map \
	-o $outpath/chr${i} \
	--chromosome=chr${i}
