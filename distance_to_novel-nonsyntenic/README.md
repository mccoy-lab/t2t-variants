# Find closest GWAS Catalog and ClinVar variants to T2T-CHM13 Novel and Non-syntenic Regions

1. Generate non-overlapping regions from novel regions and non-syntenic regions in the T2T-CHM13 assembly
	- `intersect_novel_and_nonsyntenic.sh`
	- This script separates novel and non-syntenic regions into a set of regions that are *only* annotated as novel, regions that are *only* annoated as non-syntenic, and regions that are annotated as both novel and non-syntenic

2. Concatenate non-overlapping regions into a single file
	- `python concatenate_regions.py`
	- Labels non-overlapping regions from the first step and concatenates these regions into a single bed file

3. Find the closet GWAS Catalog and ClinVar variant to each novel and non-syntenic region
	- `find_closest_vars.sh`
	- Uses `bedtools closest` to find the closest GWAS Catalog, ClinVar, and pathogenic ClinVar variant

4. Plot distances from novel and non-syntenic regions to closest GWAS Catalog and ClinVar variants
	- `python generate_distance_histogram.py`
	- Plots distances from the output of step 3, annotated by region type (novel only, non-syntenic only, or both)


# Plot lengths of T2T-CHM13 Novel and Non-syntenic Regions

1. Plot lengths of raw novel and non-syntenic regions
	- `python generate_length_histogram.py`
	- Plots the lengths of the raw novel and non-syntenic regions, used as input to step 1 in the previous section
2. Plot lengths of non-overlapping novel and non-syntenic regions
	- `python generate_no_overlap_length_histogram.py`
	- Plots the lengths of the non-overlapping novel and non-sytenic regions that were used as an input to step 3 in the previous section
	 