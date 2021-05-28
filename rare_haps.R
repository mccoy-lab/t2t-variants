library(data.table)
library(magrittr)
library(pbmcapply)

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]

# load LD data from PLINK (generated using code snippet above)
ld <- rbindlist(lapply(1:22, function(x) fread(paste0("/scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/", sample_id, "/chr", x, ".ld"))))

# compute haplotype frequency for a pair of alleles observed in specified genome
check_ld_inconsistency <- function(chrom, v1_pos, v2_pos, genome, hap = as.character(NA)) {
  
  # query SNP 1 in 1000 Genomes
  v1 <- system(paste0("/scratch4/mschatz1/rmccoy22/code/htslib-1.11/tabix ~/scratch4-mschatz1/rmccoy22/1kg-GRCh38-phase3/ALL.chr", chrom, ".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz chr", chrom, ":", v1_pos, "-", v1_pos), intern = TRUE) %>%
    strsplit(., "\t", fixed = TRUE) %>%
    unlist()
  ref_v1 <- v1[4]
  alt_v1 <- v1[5]
  v1 <- v1[10:length(v1)]
  
  # query SNP 2 in 1000 Genomes
  v2 <- system(paste0("/scratch4/mschatz1/rmccoy22/code/htslib-1.11/tabix ~/scratch4-mschatz1/rmccoy22/1kg-GRCh38-phase3/ALL.chr", chrom, ".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz chr", chrom, ":", v2_pos, "-", v2_pos), intern = TRUE) %>%
    strsplit(., "\t", fixed = TRUE) %>%
    unlist()
  ref_v2 <- v2[4]
  alt_v2 <- v2[5]
  v2 <- v2[10:length(v2)]
  
  # split the haplotypes
  dt <- data.table(v1, v2)
  dt[, c("v1_h1", "v1_h2") := tstrsplit(v1, "|", fixed = TRUE)]
  dt[, c("v2_h1", "v2_h2") := tstrsplit(v2, "|", fixed = TRUE)]
  dt_hap <- data.table(v1 = c(dt$v1_h1, dt$v1_h2), v2 = c(dt$v2_h1, dt$v2_h2))
  
  sample_id <- fread(cmd = paste0("/scratch4/mschatz1/rmccoy22/code/htslib-1.11/tabix -H ~/scratch4-mschatz1/rmccoy22/1kg-GRCh38-phase3/ALL.chr21.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz chr21 | tail -1 | tr '\t' '\n'")) %>%
    .[-(1:8),] %>%
    setnames(., "sample")
  
  dt_hap[, sample_id := rep(sample_id$sample, 2)]
  dt_hap$haplotype <- c(rep("h1", length(sample_id$sample)), rep("h2", length(sample_id$sample)))
  
  # query the test sample genome (GRCh38 is always 0 by definition)
  if (genome == "GRCh38") {
    
    return(nrow(dt_hap[(v1 == 0 & v2 == 0)]) / nrow(dt_hap))
  
  } else if (genome == "CHM13") {
    
    v1_query <- system(paste0("/scratch4/mschatz1/rmccoy22/code/htslib-1.11/tabix /scratch4/mschatz1/rmccoy22/CHM13/chm13.202000921_with38Y-align2-GRCh38.dip.vcf.gz chr", chrom, ":", v1_pos, "-", v1_pos), intern = TRUE) %>%
      strsplit(., "\t", fixed = TRUE) %>%
      unlist() %>%
      .[10:length(.)]
    v2_query <- system(paste0("/scratch4/mschatz1/rmccoy22/code/htslib-1.11/tabix /scratch4/mschatz1/rmccoy22/CHM13/chm13.202000921_with38Y-align2-GRCh38.dip.vcf.gz chr", chrom, ":", v2_pos, "-", v2_pos), intern = TRUE) %>%
      strsplit(., "\t", fixed = TRUE) %>%
      unlist() %>%
      .[10:length(.)]
      
    if (is.null(v1_query)) {
      v1_query <- 0
    } else if (grepl("^1|1", v1_query)) {
      v1_query <- 1
    }
    
    if (is.null(v2_query)) {
      v2_query <- 0
    } else if (grepl("^1|1", v2_query)) {
      v2_query <- 1
    }
    
    # return the haplotype frequency
    return(nrow(dt_hap[(v1 == v1_query & v2 == v2_query)]) / nrow(dt_hap))
    
  } else if (genome %in% dt_hap$sample_id) {
  
      v1_query <- dt_hap[sample_id == genome & haplotype == hap]$v1
      v2_query <- dt_hap[sample_id == genome & haplotype == hap]$v2
      
      return(nrow(dt_hap[sample_id != genome][(v1 == v1_query & v2 == v2_query)]) / nrow(dt_hap[sample_id != genome]))
  
  } else stop("Select a genome from: GRCh38, CHM13, or 1KGP samples.")
}

hap_freq_h1 <- pbmcmapply(check_ld_inconsistency, ld$CHR_A, ld$BP_A, ld$BP_B, genome = sample_id, hap = "h1", mc.cores = 96L)
hap_freq_h2 <- pbmcmapply(check_ld_inconsistency, ld$CHR_A, ld$BP_A, ld$BP_B, genome = sample_id, hap = "h2", mc.cores = 96L)

ld[, hap_freq_h1 := hap_freq_h1]
ld[, hap_freq_h2 := hap_freq_h2]

fwrite(ld, file = paste0("/scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/output/", sample_id, "_rare_haps.txt"),
       quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
