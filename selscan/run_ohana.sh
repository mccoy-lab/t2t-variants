#!/bin/bash
#SBATCH --job-name=runOhanaBatch
#SBATCH -N 1
#SBATCH --time=7:00:00
#SBATCH --partition=shared
#SBATCH --mem=117600MB
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu

### USAGE: sbatch run_ohana.sh <chr_num> <split_num>
### <chr_num>: number of chromosome to run on (ex: `21`)
### <split_num>: number of split VCF for that chromosome
###
### Adds variant IDs to split VCFs and converts to plink ped files for Ohana.

############################################################################

### DATA

# chromosome to run on
i=$1
# number of split VCF
num=$2

# path to Ohana bin folder
OHANAPATH=~/code/ohana/bin
# path to downsampled Ohana F matrix
F_MATRIX=chr21_pruned_50_F.matrix
# path to downsampled Ohana Q matrix
Q_MATRIX=chr21_pruned_50_Q.matrix
# directory with downsampled Ohana C matrices
C_MATRICES=c_matrices

# path to directory where vcfs are
VCFPATH=snp_genotypes/chr$i

############################################################################


ml gcc
ml openblas
ml lapack

cd /scratch/groups/rmccoy22/syan11/t2t_variants/snp_selection
echo "running ./run_ohana.sh on chr"$i "split" $num

# get path to split vcf with variant IDs
id_vcf=${VCFPATH}/chr${i}_unique_IDs_split${num}.vcf.gz

# convert .ped file into .dgm file (a.k.a. G matrix) for qpas
${OHANAPATH}/convert ped2dgm \
	${id_vcf}.ped \
	${id_vcf}.dgm
echo "vcf converted into dgm for chr"$i "and split" $num

# produce admixture-corrected AFs for full dataset by reading in pre-calculated Q matrix
# -fq: optimizer should not optimize the specified Q matrix
mkdir -p f_matrices
${OHANAPATH}/qpas \
	${id_vcf}.dgm \
	-k 8 \
	-qi $Q_MATRIX \
	-fo f_matrices/chr${i}_unique_split${num}_F.matrix \
	-e 0.0001 \
	-fq \
	-mi 50
echo "Finished inferring F matrix for chr"$i "and split" $num

# test for selection with selscan
mkdir -p selscan_out
for p_num in {1..8}
do
	${OHANAPATH}/selscan ${id_vcf}.dgm \
		f_matrices/chr${i}_unique_split${num}_F.matrix \
                ${C_MATRICES}/chr21_pruned_50_C.matrix \
                -cs ${C_MATRICES}/chr21_pruned_50_C_p${p_num}.matrix \
                > selscan_out/chr${i}_split${num}_selscan_k8_p${p_num}.out
        echo "Finished running selscan on ancestry component "$p_num
done
