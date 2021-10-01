import pysam
import pysamstats
import itertools
import statistics

# load dictionary of median autosomal coverage per sample
med_cov = {}
file1 = open('/scratch4/mschatz1/rmccoy22/rmccoy22/mask/median_cov.txt', 'r')
Lines = file1.readlines()
for line in Lines:
  line = line.strip().split()
  med_cov[line[0]] = int(line[1])

samples = list(med_cov.keys())

ref = "/scratch4/mschatz1/CHM13/drafts/t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta"

f00 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG00448.cram", "rb", reference_filename = ref)
f01 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG00449.cram", "rb", reference_filename = ref)
f02 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG00450.cram", "rb", reference_filename = ref)
f03 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG01619.cram", "rb", reference_filename = ref)
f04 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG01620.cram", "rb", reference_filename = ref)
f05 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG01621.cram", "rb", reference_filename = ref)
f06 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG01706.cram", "rb", reference_filename = ref)
f07 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG01707.cram", "rb", reference_filename = ref)
f08 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG01708.cram", "rb", reference_filename = ref)
f09 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG01770.cram", "rb", reference_filename = ref)
f10 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG01771.cram", "rb", reference_filename = ref)
f11 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG01772.cram", "rb", reference_filename = ref)
f12 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG02068.cram", "rb", reference_filename = ref)
f13 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG02069.cram", "rb", reference_filename = ref)
f14 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG02070.cram", "rb", reference_filename = ref)
f15 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG02521.cram", "rb", reference_filename = ref)
f16 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG02522.cram", "rb", reference_filename = ref)
f17 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG02523.cram", "rb", reference_filename = ref)
f18 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG02768.cram", "rb", reference_filename = ref)
f19 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG02770.cram", "rb", reference_filename = ref)
f20 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG03120.cram", "rb", reference_filename = ref)
f21 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG03121.cram", "rb", reference_filename = ref)
f22 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG03122.cram", "rb", reference_filename = ref)
f23 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG03246.cram", "rb", reference_filename = ref)
f24 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG03247.cram", "rb", reference_filename = ref)
f25 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG03248.cram", "rb", reference_filename = ref)
f26 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG03306.cram", "rb", reference_filename = ref)
f27 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG03307.cram", "rb", reference_filename = ref)
f28 = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/HG03308.cram", "rb", reference_filename = ref)


for (p00, p01, p02, p03, p04, p05, p06, p07, p08, p09, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28) in zip(
                        pysamstats.stat_coverage(f00, truncate = True, pad = True),
                        pysamstats.stat_coverage(f01, truncate = True, pad = True),
                        pysamstats.stat_coverage(f02, truncate = True, pad = True),
                        pysamstats.stat_coverage(f03, truncate = True, pad = True),
                        pysamstats.stat_coverage(f04, truncate = True, pad = True),
                        pysamstats.stat_coverage(f05, truncate = True, pad = True),
                        pysamstats.stat_coverage(f06, truncate = True, pad = True),
                        pysamstats.stat_coverage(f07, truncate = True, pad = True),
                        pysamstats.stat_coverage(f08, truncate = True, pad = True),
                        pysamstats.stat_coverage(f09, truncate = True, pad = True),
                        pysamstats.stat_coverage(f10, truncate = True, pad = True),
                        pysamstats.stat_coverage(f11, truncate = True, pad = True),
                        pysamstats.stat_coverage(f12, truncate = True, pad = True),
                        pysamstats.stat_coverage(f13, truncate = True, pad = True),
                        pysamstats.stat_coverage(f14, truncate = True, pad = True),
                        pysamstats.stat_coverage(f15, truncate = True, pad = True),
                        pysamstats.stat_coverage(f16, truncate = True, pad = True),
                        pysamstats.stat_coverage(f17, truncate = True, pad = True),
                        pysamstats.stat_coverage(f18, truncate = True, pad = True),
                        pysamstats.stat_coverage(f19, truncate = True, pad = True),
                        pysamstats.stat_coverage(f20, truncate = True, pad = True),
                        pysamstats.stat_coverage(f21, truncate = True, pad = True),
                        pysamstats.stat_coverage(f22, truncate = True, pad = True),
                        pysamstats.stat_coverage(f23, truncate = True, pad = True),
                        pysamstats.stat_coverage(f24, truncate = True, pad = True),
                        pysamstats.stat_coverage(f25, truncate = True, pad = True),
                        pysamstats.stat_coverage(f26, truncate = True, pad = True),
                        pysamstats.stat_coverage(f27, truncate = True, pad = True),
                        pysamstats.stat_coverage(f28, truncate = True, pad = True)):
  mean_norm_cov = statistics.mean([
          p00['reads_all'] / med_cov[samples[0]], 
          p01['reads_all'] / med_cov[samples[1]], 
          p02['reads_all'] / med_cov[samples[2]],
          p03['reads_all'] / med_cov[samples[3]], 
          p04['reads_all'] / med_cov[samples[4]], 
          p05['reads_all'] / med_cov[samples[5]],
          p06['reads_all'] / med_cov[samples[6]], 
          p07['reads_all'] / med_cov[samples[7]], 
          p08['reads_all'] / med_cov[samples[8]],
          p09['reads_all'] / med_cov[samples[9]], 
          p10['reads_all'] / med_cov[samples[10]], 
          p11['reads_all'] / med_cov[samples[11]],
          p12['reads_all'] / med_cov[samples[12]], 
          p13['reads_all'] / med_cov[samples[13]], 
          p14['reads_all'] / med_cov[samples[14]],
          p15['reads_all'] / med_cov[samples[15]], 
          p16['reads_all'] / med_cov[samples[16]], 
          p17['reads_all'] / med_cov[samples[17]],
          p18['reads_all'] / med_cov[samples[18]], 
          p19['reads_all'] / med_cov[samples[19]], 
          p20['reads_all'] / med_cov[samples[20]],          
          p21['reads_all'] / med_cov[samples[21]], 
          p22['reads_all'] / med_cov[samples[22]], 
          p23['reads_all'] / med_cov[samples[23]],
          p24['reads_all'] / med_cov[samples[24]], 
          p25['reads_all'] / med_cov[samples[25]], 
          p26['reads_all'] / med_cov[samples[26]],
          p27['reads_all'] / med_cov[samples[27]], 
          p28['reads_all'] / med_cov[samples[28]]]);
  if (mean_norm_cov > 0.75 and mean_norm_cov < 1.25):
    print(p00['chrom'], p00['pos'] + 1, mean_norm_cov)
