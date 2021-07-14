#!/bin/sh

novel_and_nonsyntenic_regions=chm13_v1.0_novel_and_nonsyntenic_and_intersection.all.bed
gwas_positions=gwas_catalog_v1.0-lifted.all.CHM13.sorted.bed
clinvar_positions=clinvar_20210308.CHM13.lifted.positions.noY.sorted.bed
clinvar_pathogenic_positions=clinvar_20210308.CHM13.lifted.pathogenic.positions.noY.sorted.bed

closest_gwas_out=closest_gwas_catalog_to_each_novel_region.tsv
closest_clinvar_out=closest_clinvar_to_each_novel_region.tsv
closest_clinvar_pathogenic_out=closest_pathogenic_clinvar_to_each_novel_region.tsv

bedtools closest -t first -a ${novel_and_nonsyntenic_regions} -b ${gwas_positions} > ${closest_gwas_out}
bedtools closest -t first -a ${novel_and_nonsyntenic_regions} -b ${clinvar_positions} > ${closest_clinvar_out}
bedtools closest -t first -a ${novel_and_nonsyntenic_regions} -b ${clinvar_pathogenic_positions} > ${closest_clinvar_pathogenic_out}