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

When complete:
```
cat *.median >  median_cov.txt

sbatch coverage_mask_wrapper.sh
```

When complete:
```
# concatenate and bedtools merge individual chromosome masks

# merge the CRAM files in preparation for steps 2 and 3
samtools merge \
  -h HG00448.cram \
  -O CRAM \
  --reference /scratch4/mschatz1/CHM13/drafts/t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta \
  -@ 48 \
  merged.cram \
  *.cram

samtools index -@ 48 merged.cram
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