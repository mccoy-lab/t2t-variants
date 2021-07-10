#!/bin/bash
#SBATCH --job-name=calculateLD
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu

cd /scratch4/mschatz1/rmccoy22/syan11/gwas_ld

ml r/4.0.2

ident=$1
Rscript calculate_ld.R $ident
