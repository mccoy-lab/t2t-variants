#!/bin/bash

### USAGE: ./02_prep_vcfs_wrapper.sh
###
### Wrapper function to submit slurm jobs for prepping split VCFs for Ohana.

############################################################################

### DATA

# path to where directories of split vcfs are
VCFS=snp_genotypes

############################################################################

# submit batch jobs for all split vcfs
for i in {1..22} X; do
	# get split numbers of vcfs
	all_splits=($(ls ${VCFS}/chr${i}/chr${i}_unique_split* | cut -f 1 -d "." | cut -f 3 -d "t"))
	# submit batch job for each split vcf for this chromosome
	for split_num in ${all_splits[@]}; do
		echo $split_num
		sbatch prep_split_vcfs.sh $i $split_num
	done
done
