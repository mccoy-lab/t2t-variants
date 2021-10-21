import sys
import gzip


#==================#
# Define Filepaths #
#==================#

recovered_id_file = sys.argv[1]
rejected_vcf_file = sys.argv[2]
recovered_lifted_vcf_file = sys.argv[3]
output_file = open(sys.argv[4], 'w')

#========================#
# Load in recovered info #
#========================#

recovered_info = {}

for line in open(recovered_id_file):
	fields = line.rstrip('\n').split('\t')
	rsid = fields[2]
	if rsid in recovered_info:
		print(rsid, 'repeated in recovered')
	else:
		recovered_info[rsid] = fields

#==================#
# Define functions #
#==================#

def rev_comp(seq):
	mapping = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
	out_seq = ''
	for char in seq:
		out_seq = mapping[char] + out_seq
	return out_seq


#============================#
# Resolve recovered variants #
#============================#

rejected_rsids = set()

for line in gzip.open(recovered_lifted_vcf_file):
	line = line.decode()
	if line.startswith('#'):
		output_file.write(line)

for line in gzip.open(rejected_vcf_file, mode='rt', compresslevel=5):
	# line = line.decode()
	if line.startswith('#'):
		continue
	fields = line.rstrip('\n').split('\t')
	rsid = fields[2]
	if rsid in recovered_info:
		recovered_chrom = recovered_info[rsid][0]
		recovered_pos = recovered_info[rsid][1]
		recovered_ref = recovered_info[rsid][3]
		recovered_alt = recovered_info[rsid][4]

		rejected_ref = fields[3]
		rejected_alt = fields[4].split(',')
		
		rejected_alleles = [rejected_ref] + rejected_alt

		orig_alleles = [field for field in recovered_info[rsid][5].split(';') if field.startswith('AttemptedAlleles')][0][17:].split('*->')
		orig_alleles = set([orig_alleles[0]] + orig_alleles[1].split(','))

		if rejected_ref == 'N' and set([rev_comp(x) for x in rejected_alleles]) == set(orig_alleles):

			rejected_alt = [rev_comp(x) for x in rejected_alt]
			rejected_alt.remove(recovered_ref)
			rejected_alt.append(rev_comp(rejected_ref))

			output_file.write('\t'.join([recovered_chrom, recovered_pos, rsid] + [recovered_ref, ','.join(sorted(rejected_alt))] + fields[5:]) + '\n')

		elif rejected_ref != 'N' and rejected_ref == rev_comp(recovered_alt) : # Checks whether the variant lifts over to the opposite strand

			rejected_alt = [rev_comp(x) for x in rejected_alt]
			rejected_alt.remove(recovered_ref)
			rejected_alt.append(rev_comp(rejected_ref))

			output_file.write('\t'.join([recovered_chrom, recovered_pos, rsid] + [recovered_ref, ','.join(sorted(rejected_alt))] + fields[5:]) + '\n')

		else:
			rejected_alt.remove(recovered_ref)
			rejected_alt.append(rejected_ref)

			output_file.write('\t'.join([recovered_chrom, recovered_pos, rsid] + [recovered_ref, ','.join(sorted(rejected_alt))] + fields[5:]) + '\n')

		if rsid in rejected_rsids:
			print(rsid, 'repeated in recovered')
		rejected_rsids.add(rsid)

output_file.close()
