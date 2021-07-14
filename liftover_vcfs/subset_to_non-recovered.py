import sys
import gzip

recovered_vcf = sys.argv[1]
total_rejected_vcf = sys.argv[2]
out_vcf = open(sys.argv[3], 'w')

recovered_rsids = set()

for line in open(recovered_vcf):
	if line.startswith('#'):
		continue
	fields = line.rstrip('\n').split('\t')
	rsid = fields [2]
	if rsid in recovered_rsids:
		print(rsid, 'repeated in recovered')
	recovered_rsids.add(rsid)


for line in gzip.open(total_rejected_vcf):
	line = line.decode()
	if line.startswith('#'):
		out_vcf.write(line)
		continue
	fields = line.rstrip('\n').split('\t')
	rsid = fields[2]
	failure_mode = fields[6]

	if rsid not in recovered_rsids and failure_mode == 'MismatchedRefAllele':
		out_vcf.write(line)

out_vcf.close()