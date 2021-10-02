The approach below generates three genome masks, according to the three criteria described in the following article: https://www.illumina.com/science/genomics-research/articles/identifying-genomic-regions-with-high-quality-single-nucleotide-.html

These masks define regions of the genome where SNV discovery and genotyping is relatively reliable. The three masks are then intersected to create one combined mask. Note that step 1 must be run separately for the autosomes and sex chromosomes.

```
module load samtools

cd /scratch4/mschatz1/rmccoy22/rmccoy22/mask
```

## 1. Normalized depth mask

### Compute median autosomal coverage
```
sbatch get_median_coverage.sh
```
Technically, the above value may differ slightly from the median, as `mosdepth` does not report the median itself, but the cumulative distribution of coverages. We extract the first reported coverage for which the reported cumulative proportion of total bases exceeds 0.5. 

### Compare observed coverage to median autosomal coverage

For sex chromosomes, adjust the median expectation based on the expected dosage of the chromosome for that sample. For example, for females (XX), the median autosomal coverage is used for the X, while for males (XY), the median is divided by 2.

When complete:
```
cat *.median >  median_cov.txt

sbatch coverage_mask_wrapper.sh
```

When complete:
```
# concatenate and bedtools merge individual chromosome masks
for i in {1..22}
do
    cat coverage_mask_chr${i}.txt | bedtools merge;
done > coverage_mask.bed

mkdir coverage_mask_per_chr
mv coverage_mask_chr*.txt coverage_mask_per_chr

# merge the CRAM files in preparation for steps 2 and 3
for i in {1..22}
do
samtools merge \
  -R chr${i} \
  -h HG00448.cram \
  -O CRAM \
  --reference /scratch4/mschatz1/CHM13/drafts/t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta \
  -@ 4 \
  chr${i}_merged.cram \
  HG*.cram &
done

wait

for i in {1..22}
do
  samtools index -@ 4 chr${i}_merged.cram &
done
```

## 2. Mapping quality mask

```
sbatch mapq_mask_wrapper.sh
```
When complete:
```
# concatenate and bedtools merge individual chromosome masks
```


## 3. Base quality mask

```
sbatch baseq_mask_wrapper.sh
```
When complete:
```
# concatenate and bedtools merge individual chromosome masks
```

## 4. Combined mask
```
bedtools intersect -a coverage_mask.bed -b mapq_mask.bed baseq_mask.bed > combined_mask.bed
```
