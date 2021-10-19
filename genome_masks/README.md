The approach below generates three genome masks, according to the three criteria described in the following article: https://www.illumina.com/science/genomics-research/articles/identifying-genomic-regions-with-high-quality-single-nucleotide-.html

These masks define regions of the genome where SNV discovery and genotyping is relatively reliable. The three masks are then intersected to create one combined mask.

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

For sex chromosomes, adjust the median expectation based on the dosage of the chromosome for that sample. For example, for females (XX), we simply use the median autosomal coverage for the X. In contrast, for males (XY), the median autosomal coverage is divided by 2 for both the X and the Y. Coverage for female samples is not considered when producing the coverage mask for the Y chromosome.

When complete:
```
cat *.median >  median_cov.txt

sbatch coverage_mask_wrapper.sh
```

When complete:
```
# concatenate and bedtools merge individual chromosome masks
for i in {1..22} X Y;
do
    cat coverage_mask_chr${i}.txt | bedtools merge >> coverage_mask.bed;
done

mkdir coverage_mask_per_chr
mv coverage_mask_chr*.txt coverage_mask_per_chr
```

## 2. Mapping quality mask
```
# merge the CRAM files in preparation for steps 2 and 3
cd /scratch4/mschatz1/mschatz/T2T/2021.09.29.mask
for i in {1..22} X Y
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

# index the CRAM files
for i in {1..22} X Y
do
  samtools index -@ 4 chr${i}_merged.cram &
done
```

```
sbatch mapq_mask_wrapper.sh
```
When complete:
```
# concatenate and bedtools merge individual chromosome masks
for i in {1..22} X Y;
do
    cat mapq_mask_chr${i}.txt | bedtools merge >> mapq_mask.bed;
done

mkdir mapq_mask_per_chr
mv mapq_mask_chr*.txt mapq_mask_per_chr
```


## 3. Base quality mask

```
sbatch baseq_mask_wrapper.sh
```
When complete:
```
# concatenate and bedtools merge individual chromosome masks
for i in {1..22} X Y;
do
    cat baseq_mask_chr${i}.txt | bedtools merge >> baseq_mask.bed;
done

mkdir baseq_mask_per_chr
mv baseq_mask_chr*.txt baseq_mask_per_chr
```

## 4. Combined mask
```
bedtools multiinter \
  -header \
  -i baseq_mask.bed coverage_mask.bed mapq_mask.bed \
  -names baseq coverage mapq |\
  awk '{if ($4 == 3) print $1"\t"$2"\t"$3}' > combined_mask.bed
```
