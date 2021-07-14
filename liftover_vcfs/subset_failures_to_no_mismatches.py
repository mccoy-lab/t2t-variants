import sys
import gzip

for line in gzip.open(sys.argv[1]):
	line = line.decode()
	if line.startswith('#'):
		print(line.rstrip('\n'))
		continue
	fields = line.rstrip('\n').split('\t')
	if fields[6] != 'MismatchedRefAllele':
		print(line.strip('\n'))