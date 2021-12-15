# Find closest GWAS Catalog and ClinVar variants to T2T-CHM13 Previously Unresolved and Non-syntenic Regions

**NOTE**: Any references to "novel" regions in this README file and any scripts are referring to regions that were previously unresolved in the GRCh38 assembly.

The regions of the T2T-CHM13 v1.0 assembly that were previously-unresolved (novel) in the GRCh38 assembly are described in this file: `chm13_v1.0_uncoveredByGRCh38WinnowmapAlignments.bed`
* This track is available in the UCSC assembly hub browser, and is described [here](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1233948703_ZFOcFP7TeC49XE5ajuTQc7BgMgal&g=hub_2395475_chm13_uncovered_byGRCh38)
* A BigBed of this track can be downloaded here: [http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/grch38NonSyntenic/chm13_v1.0_uncoveredByGRCh38WinnowmapAlignments.bigBed](http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/grch38NonSyntenic/chm13_v1.0_uncoveredByGRCh38WinnowmapAlignments.bigBed)

The 1Mbp regions of the T2T-CHM13 v1.0 assembly that are non-syntenic with the GRCh38 assembly are described in this file: `chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed`
* This file is available in this repo, and is in the process of being added to the UCSC assembly hub browser.

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
	 
