import sys
import pysam
import statistics

chrom = sys.argv[1]
ref = "/scratch4/mschatz1/CHM13/drafts/t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta"

# to do: swap in the merged file below
aln = pysam.AlignmentFile("/scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/" + chrom + "_merged.cram", "rb", reference_filename = ref)

for rec in aln.pileup(chrom, truncate = True):
  baseq = rec.get_query_qualities()
  if (len(baseq) > 0):
    q20_baseq = [i for i in baseq if i > 20]
    q20_prop = len(q20_baseq) / len(baseq)
    if (q20_prop >= 0.9):
      print(rec.reference_name, rec.pos, rec.pos + 1, q20_prop, sep = "\t")
