# t2t-variants

Code associated with the manuscript "A complete reference genome improves analysis of human genetic variation". This repository contains downstream analysis code used in various population genetic analyses, organized into directories described below:

`count_singletons` - counting the number of singleton alleles in 1000 Genomes Project samples, as well as in the T2T-CHM13 genome itself (i.e., fixed among 1000 Genomes Project samples).

`distance_to_novel-nonsyntenic` - determining the distance of T2T-CHM13 novel and non-syntenic regions to the closest GWAS Catalog or ClinVar variants and investigating the distribution of the lengths of these regions.

`genome_masks` - identifying regions of the T2T-CHM13 reference genome where SNV calling is relatively reliable.

`gwas_ld` - identifying variants in the GWAS Catalog that are in linkage disequilibrium (LD) with variants in T2T-CHM13 non-syntenic regions.

`haplotype_structure` - investigation of LD-discordant alleles (i.e., alleles in negative linkage disequilibrium) which may arise in GRCh38 by consequence of it's construction with BAC clone libraries derived from multiple donors.

`liftover_vcfs` - liftover variants from GRCh38 to T2T-CHM13, recover biallelic and multiallelic ref/alt switches, and intersect variants that failed to lift over with T2T-CHM13 indels relative to GRCh38.

`local_ancestry` - superpopulation-level inference of local ancestry and Neanderthal introgression for GRCh38 and CHM13, in comparison to reference haplotypes from the 1000 Genomes Project.

`selscan` - searching for variants in T2T-CHM13 non-syntenic and novel regions that show extreme differences in allele frequency across the 1000 Genomes populations (using [Ohana](https://github.com/jade-cheng/ohana)).
