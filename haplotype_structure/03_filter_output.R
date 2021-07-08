library(data.table)
library(tidyverse)
library(bedr)
library(pbmcapply)

# make sure to export bedtools path
# export PATH=${PATH}:"/scratch4/mschatz1/rmccoy22/code/"

setwd("/scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/output/")

count_rare_haps <- function(sample_id, haplotype) {
  
  dt <- rbindlist(lapply(1:22, function(x) fread(paste0(sample_id, "_", x, "_", haplotype, "_rare_haps.txt")))) %>%
    .[hap_freq == 0] %>%
    .[, chr := paste0("chr", CHR_A)]
  
  dt_A <- dt[, c("chr", "BP_A", "BP_A")] %>%
    setnames(., c("chr", "start", "end")) %>%
    .[, index := .I]
  dt_B <- dt[, c("chr", "BP_B", "BP_B")] %>%
    setnames(., c("chr", "start", "end")) %>%
    .[, index := .I]
  
  # overlap with indel intervals - to be removed
  
  indels <- fread("/scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/filter_bed/GRCh38_CHM13_indels.dip.bed") %>%
    setnames(., c("chr", "start", "end"))
  
  setkey(indels, chr, start, end)
  
  olaps_A <- foverlaps(dt_A, indels, type = "within", which = TRUE, mult = "first")
  olaps_B <- foverlaps(dt_B, indels, type = "within", which = TRUE, mult = "first")
  
  dt[, olaps_indels_A := olaps_A]
  dt[, olaps_indels_B := olaps_B]
  
  # overlap with dip.bed callable intervals - to be kept
  
  dip <- fread("/scratch4/mschatz1/rmccoy22/rmccoy22/rare_haplotypes/filter_bed/chm13.202000921_with38Y-align2-GRCh38.dip.bed") %>%
    setnames(., c("chr", "start", "end"))
  
  setkey(dip, chr, start, end)
  
  olaps_A <- foverlaps(dt_A, dip, type = "within", which = TRUE, mult = "first")
  olaps_B <- foverlaps(dt_B, dip, type = "within", which = TRUE, mult = "first")
  
  dt[, olaps_dip_A := olaps_A]
  dt[, olaps_dip_B := olaps_B]
  
  dt <- dt[is.na(olaps_indels_A) & is.na(olaps_indels_B) & !is.na(olaps_dip_A) & !is.na(olaps_dip_B)]
  
  dt_bed <- paste0(dt$chr, ":", dt$BP_A, "-", dt$BP_B)
  dt_bed_merged <- bedr.merge.region(dt_bed, verbose = FALSE)
  
  return(
    data.table(sample_id = sample_id, 
               haplotype = haplotype, 
               n_discordant_snps = length(dt_bed), 
               n_discordant_haplotypes = length(dt_bed_merged))
  )
}

file_list_h1 <- table(unlist(map(strsplit(list.files(pattern = "_h1_"), "_"), 1)))
sample_list_h1 <- names(file_list_h1[file_list_h1 == 22])
h1_res <- rbindlist(pbmclapply(sample_list_h1, function(x) count_rare_haps(x, "h1"), mc.cores = getOption("mc.cores", 48L)))

file_list_h2 <- table(unlist(map(strsplit(list.files(pattern = "_h2_"), "_"), 1)))
sample_list_h2 <- names(file_list_h2[file_list_h2 == 22])
h2_res <- rbindlist(pbmclapply(sample_list_h2, function(x) count_rare_haps(x, "h2"), mc.cores = getOption("mc.cores", 48L)))

results <- rbind(h1_res, h2_res) %>%
  setorder(sample_id, haplotype)

fwrite(results, file = "~/rh_summary.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
