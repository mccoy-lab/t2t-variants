#!/bin/bash

### USAGE: ./03_calculate_ld_wrapper.sh
###
### Submits slurm jobs to calculate LD between all lifted over GWAS SNPs
### and novel CHM13 snps for a given ancestry.

############################################################################


ml r/4.0.2
ml r-data-table/1.12.8

for ident in ACB.ASW CHS.CHB.JPT BEB.PJL.GIH.ITU.STU FIN.TSI.IBS.CEU.GBR PUR.CLM.PEL.MXL GWD.MSL.ESN.YRI.LWK
do
	sbatch calculate_ld_slurmjob.sh $ident
done
