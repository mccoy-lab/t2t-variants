import sys

for line in sys.stdin:
	if line.startswith('#'):
		print(line.strip())
		continue
	fields = line.strip().split('\t')
	if ',' in fields[4]:
		for alt in fields[4].split(','):
			print('\t'.join(fields[:4] + [alt] + fields[5:]))
	else:
		print(line.strip())