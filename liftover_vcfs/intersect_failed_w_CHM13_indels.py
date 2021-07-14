import gzip
import sys

failed_positions_file = sys.argv[1] # This should be a vcf.gz file
chm13_on_grch38_file = '/scratch4/mschatz1/rmccoy22/CHM13/chm13.202000921_with38Y-align2-GRCh38.dip.vcf.gz' # Change to sys.argv[2]?


def intersects(range1, range2):
	if max(range1[0], range2[0]) < min(range1[1], range2[1]):
		return True
	return False


chm13_indels_by_chr = {}

for line in gzip.open(chm13_on_grch38_file):
	line = line.decode('utf-8')
	if line.startswith('#'):
		continue
	fields = line.rstrip('\n').split('\t')
	chrom = fields[0]
	chm13_indels_by_chr.setdefault(chrom, {})
	start_pos = int(fields[1])
	end_pos = start_pos + len(fields[3])
	if len(fields[3]) != len(fields[4]):
		chm13_indels_by_chr[chrom][(start_pos, end_pos)] = (fields[3], fields[4])

indel_count = 0
total_count = 0
chm13_alt_count = 0

for line in open(failed_positions_file):
	if line.startswith('#'):
		continue
	fields = line.strip().split('\t')
	chrom = fields[0]
	chm13_indels_by_chr.setdefault(chrom, set())
	start_pos = int(fields[1])
	end_pos = start_pos + len(fields[3])

	rsid = fields[2]

	indel_overlap_check = 0
	chm13_alt_check = 0

	ref_allele = fields[3]
	alt_alleles = fields[4].split(',')

	for indel_pos in chm13_indels_by_chr[chrom]:
		if intersects((start_pos, end_pos), indel_pos):
			indel_overlap_check = 1
			for alt_allele in alt_alleles:
				if alt_allele == chm13_indels_by_chr[chrom][indel_pos][1] and ref_allele == chm13_indels_by_chr[chrom][indel_pos][0] and (start_pos, end_pos) == indel_pos:
					chm13_alt_check = 1
					break

	indel_overlap_count += indel_overlap_check
	total_count += 1
	chm13_alt_count += chm13_alt_check

print(f"""\n{total_count} total variants analyzed
	  {indel_overlap_count} of these variants intersect indels in the CHM13-on-GRCh38 vcf file
	  {chm13_alt_count} of the variants that intersect indels in the CHM13-on-GRCh38 vcf file
	  are cases in which CHM13 has (one of) the alternate alleles\n""")