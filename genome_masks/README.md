```
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
# concatenate individual chromosome masks
# bedtools merge
```

## 2. Mapping quality mask


## 3. Base quality mask

## 4. Combined mask
```
# bedtools intersect
```
