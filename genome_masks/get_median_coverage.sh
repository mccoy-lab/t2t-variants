#!/bin/bash

#SBATCH --job-name=median_cov
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=2:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=rmccoy22@jhu.edu
#SBATCH --array=1-29%29

cd /scratch4/mschatz1/rmccoy22/rmccoy22/mask

SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" sample_list.txt | awk '{print $1}'`

/scratch4/mschatz1/rmccoy22/code/mosdepth \
    -n \
    -t 48 \
    --use-median \
    -f /scratch4/mschatz1/CHM13/drafts/t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta \
    -b /scratch4/mschatz1/rmccoy22/rmccoy22/mask/autosomes.bed \
    ${SAMPLE} \
    /scratch4/mschatz1/mschatz/T2T/2021.09.29.mask/${SAMPLE}.cram

cat ${SAMPLE}.mosdepth.region.dist.txt | awk '{if ($1 == "total" && $3 > 0.5) print}' | head -1 | awk -v "sample=${SAMPLE}" '{print sample"\t"$2}' > ${SAMPLE}.median

rm ${SAMPLE}*.txt ${SAMPLE}*.csi ${SAMPLE}*.gz
