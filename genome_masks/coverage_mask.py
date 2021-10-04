import sys
import pysam
import pysamstats
import itertools
import statistics

chrom = sys.argv[1]
ref = "/scratch4/mschatz1/CHM13/drafts/t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta"
min_norm_cov = 0.75
max_norm_cov = 1.25

# load dictionary of median autosomal coverage per sample
med_cov = {}
cov_file = open('/scratch4/mschatz1/rmccoy22/rmccoy22/mask/median_cov.txt', 'r')
Lines = cov_file.readlines()
for line in Lines:
  line = line.strip().split()
  med_cov[line[0]] = int(line[1])

samples = list(med_cov.keys())

# adjust expectations if sex chromosome
sample_sex = {}
sex_file = open('/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/samples.info', 'r')
Lines = sex_file.readlines()[1:]
for line in Lines:
  line = line.strip().split()
  sample_sex[line[0]] = line[2]

if (chrom == "chrX"):
  for sample in samples:
    if (sample_sex[sample] == "F"):
      continue;
    elif (sample_sex[sample] == "M"):
      med_cov[sample] = med_cov[sample] / 2;

if (chrom == "chrY"):
  for sample in list(sample_sex.keys()):
    if (sample_sex[sample] == "F"):
      # print(sample, "is female");
      samples.remove(sample);
    elif (sample_sex[sample] == "M"):
      med_cov[sample] = med_cov[sample] / 2;

pys_list = []
for idx, sample in enumerate(samples):
  globals()["f" + str(idx)] = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/" + sample + ".cram", "rb", reference_filename = ref)
  globals()["p" + str(idx)] = pysamstats.stat_coverage(eval("f" + str(idx)), chrom = chrom, truncate = True, pad = True)
  pys_list.append(eval("p" + str(idx)))

for rec in zip(*pys_list):
  norm_cov = []
  for idx, pys in enumerate(rec):
    sample = samples[idx]
    norm_cov_sample = pys['reads_all'] / med_cov[sample]
    norm_cov.append(norm_cov_sample)
  mean_norm_cov = statistics.mean(norm_cov)
  if (mean_norm_cov > 0.75 and mean_norm_cov < 1.25):
    print(pys['chrom'], pys['pos'], pys['pos'] + 1, mean_norm_cov, sep = "\t")
