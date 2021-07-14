#!/usr/bin/env python3

#=====================#
# Define region files #
#=====================#

# These are non-overlapping regions generated from "bedtools intersect" between 
# the novel and non-syntenic regions, and subsequent "bedtools subtract" operations

novel_only_bed = 'chm13_v1.0_uncoveredByGRCh38WinnowmapAlignments-ONLY.bed'
nonsyntenic_only_bed = 'chm13.draft_v1.0_plus38Y.no_snyteny_1Mbp-ONLY.bed'
intersection_bed = 'chm13_v1.0_winnowmap-novel_and_no_synteny_1Mbp-INTERSECTION.bed'


#===============#
# Merge regions #
#===============#

out_file = open('chm13_v1.0_novel_and_nonsyntenic_and_intersection.all.bed','w')

for line in open(novel_only_bed):
	if line.startswith('#'):
		continue
	line = line.rstrip('\n') + '\tNOVEL\n'
	out_file.write(line)

for line in open(nonsyntenic_only_bed):
	if line.startswith('#'):
		continue
	line = line.rstrip('\n') + '\tNONSYNTENIC\n'
	out_file.write(line)

for line in open(intersection_bed):
	if line.startswith('#'):
		continue
	line = line.rstrip('\n') + '\tBOTH\n'
	out_file.write(line)

out_file.close()
