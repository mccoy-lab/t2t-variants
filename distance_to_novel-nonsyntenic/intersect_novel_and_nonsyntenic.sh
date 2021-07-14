#!/bin/sh

novel_regions=chm13_v1.0_uncoveredByGRCh38WinnowmapAlignments.bed
nonsyntenic_regions=chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp.bed
intersected_file=chm13_v1.0_winnowmap-novel_and_no_synteny_1Mbp-INTERSECTION.bed

novel_regions_only=chm13_v1.0_uncoveredByGRCh38WinnowmapAlignments-ONLY.bed
nonsyntenic_regions_only=chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp-ONLY.bed

bedtools intersect -a ${novel_regions} -b ${nonsyntenic_regions} > ${intersected_file}

bedtools subtract -a ${novel_regions} -b ${intersected_file} > ${novel_regions_only}
bedtools subtract -a ${nonsyntenic_regions} -b ${intersected_file} > ${nonsyntenic_regions_only}

python merge_regions.py
