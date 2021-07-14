# Lifting over variants from GRCh38 to T2T-CHM13

Liftover is performed with [GATK LiftoverVCF](https://gatk.broadinstitute.org/hc/en-us/articles/360037060932-LiftoverVcf-Picard-) tool

1. Perform initial liftover
	- `liftover_vcf.sh <in_vcf> <out_lifted> <out_rejected>`

2. Recover ref/alt switches that failed initial liftover
	- GATK Liftover does not natively support recovering multiallelic alleles, so multiallelic alleles are separated into biallelic variants, lifted over, and compressed
	- `resolve_refalt_switches_in_multialt.sh <initial_rejected_vcf> <out_recovered>`
	- This script uses the `mult_alt_resolver.py` and `parse_mult_alt_alleles.py` scripts

3. Subset reference mismatch variants whose position lifts over, but are not recoverable with GATK Liftover
	- `subset_failures_to_non-recovered.sh <recovered_vcf> <initial_rejected_vcf> <out_refmismatch>`
	- This script uses the `subset_failures_to_non-recovered.py` script

4. Subset variants whose positions failed to lift over
	- `subset_failures_to_no-mismatches.sh <initial_rejected_vcf> <out_position_failed>`
	- This script uses the `subset_failures_to_no-mismatches.py` script

5. Count intersections with CHM13 indels, relative to GRCh38, from variants whose position failed to lift over
	- `python intersect_failed_w_CHM13_indels.py <position_failed_vcf>`
	- This script will output the number of input variants that intersected an indel on CHM13 relative to GRCh38. It will also count the number of these variants whose alternate allele (or one of the possible alternate alleles) *is* the CHM13 indel with which it intersected
	