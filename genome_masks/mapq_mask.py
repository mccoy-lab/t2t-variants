import sys
import pysam
import statistics

chrom = sys.argv[1]
ref = "/scratch4/mschatz1/CHM13/drafts/t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta"

aln = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/" + chrom + "_merged.cram", "rb", reference_filename = ref)

for rec in aln.pileup(chrom, truncate = True):
  mapq = rec.get_mapping_qualities()
  if (len(mapq) > 0):
    mean_mapq = statistics.mean(mapq)
    if (mean_mapq >= 50):
      print(rec.reference_name, rec.pos, rec.pos + 1, mean_mapq, sep = "\t")
