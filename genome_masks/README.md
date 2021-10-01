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
```

## 2. Mapping quality mask


## 3. Base quality mask

## 4. Combined mask
```
# bedtools intersect
```
