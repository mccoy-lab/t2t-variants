# t2t-variants

Code associated with the manuscript "A complete reference genome improves analysis of human genetic variation". This repository contains downstream analysis code used in various population genetic analyses by the McCoy lab, organized into directories described below:

`gwas_ld` - identifying variants in the GWAS Catalog that are in linkage disequilibrium (LD) with variants in CHM13 non-syntenic regions.

`haplotype_structure` - investigation of LD-discordant alleles (i.e., alleles in negative linkage disequilibrium) which may arise in GRCh38 by consequence of it's construction with BAC clone libraries derived from multiple donors.

`local_ancestry` - superpopulation-level inference of local ancestry and Neanderthal introgression for GRCh38 and CHM13, in comparison to reference haplotypes from the 1000 Genomes Project.

`selscan` - searching for variants in CHM13 non-syntenic and novel regions that show extreme differences in allele frequency across the 1000 Genomes populations (using [Ohana](https://github.com/jade-cheng/ohana)).
