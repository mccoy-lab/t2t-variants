#!/bin/bash

### USAGE: ./03_ohana_wrapper.sh
###
### Wrapper function to submit slurm jobs for running Ohana on split VCFs.

############################################################################

### DATA

# path to where directories of split vcfs are
VCFS=snp_genotypes

############################################################################

# submit batch jobs for all split vcfs
for i in {1..22} X; do
        # get split numbers of vcfs
        all_splits=($(ls ${VCFS}/chr${i}/chr${i}_unique_IDs_split* | cut -f 1 -d "." | cut -f 3 -d "t"))
        # submit batch job for each split vcf for this chromosome
        for split_num in ${all_splits[@]}; do
                echo "chr"$i "split"$split_num
                sbatch run_ohana.sh $i $split_num
        done
done