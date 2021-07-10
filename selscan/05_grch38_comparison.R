library(data.table)
library(dplyr)
library(tidyr)
library(pbmcapply)

setwd("/scratch/groups/rmccoy22/syan11/t2t_variants/snp_selection/")

############################################################################

### DATA

# path to file of selscan combined and filtered output
selscan_filt_path <- "selscan_filt_genes.txt.gz"

# path to ucsc liftover tool
liftover_path <- "~/code/liftOver"
# chain file for CHM13 v1 -> GRCh38
chm13_to_hg38_chain <- "/scratch/groups/rmccoy22/syan11/t2t_variants/t2t-chm13-v1.0.hg38.over.chain"

# path to directory of selscan output files from GRCh38 variants
grch38_selscan_path <- "/scratch/groups/rmccoy22/syan11/sv_selection/ohana/k8/1KGP_SNP/"

# output directory for GRCh38 selscan results
grch38_combined_dir <- "grch38_selscan_combined/"
# path to htslib install
htslib_path <- "~/code/htslib-1.11/"

############################################################################


### lift over top hits to compare to GRCh38 results

selscan_perc_filt <- fread(selscan_filt_path)

# convert selscan output into a bed-like file for liftover
top_hits <- selscan_perc_filt[perc > 0.999] %>%
  .[, snp_start := pos - 1] %>%
  .[, snp_end := pos]
fwrite(top_hits, "selscan_perc_filt.bed",
       sep = "\t", na = "NA", col.names = FALSE, quote = FALSE)

# attempt to lift to grch38 to see what fails
system(paste(liftover_path,
             "selscan_perc_filt.bed",
             chm13_to_hg38_chain,
             "selscan_perc_filt_lifted.bed",
             "selscan_perc_filt_unmapped.bed",
             "-bedPlus=3"),
       intern = TRUE) # 2144 out of 5575 failed

# failed variants
failed <- fread(cmd = "grep -v '#' top_hits_failed.txt",
                col.names = c("chr", "snp_start", "snp_end", "ID", "step", "lle_ratio",
                              "ancestry_component", "perc", "gene_id", "gene_type", "gene_name"))
# successfully lifted variants
lifted <- fread("selscan_perc_filt_lifted.bed",
                col.names = names(selscan_perc_filt))


############################################################################

### read in selscan results from grch38

get_grch38_snps <- function(p, chr) {
  selscan_snp <- fread(paste0(grch38_selscan_path, "selscan_rerun/chr",
                              chr, "_selscan_k_p", p, ".out"))
  snp_map <- fread(paste0(grch38_selscan_path, "vcfs/chr",
                          chr, ".map"), header = FALSE) %>%
    setnames(., c("chr", "ID", "drop", "pos")) %>%
    .[, -3, with = FALSE]
  selscan_snp <- cbind(snp_map, selscan_snp)
  selscan_snp <- selscan_snp[`global-lle` > -1000] %>%
    setorder(., -`lle-ratio`) %>%
    setnames(., "lle-ratio", "lle_ratio") %>%
    .[, ancestry_component := p] %>%
    .[, chr := paste0("chr", chr)]
  return(selscan_snp)
}
# wrapper function for combining SNPs for each chromosome across ancestry components
get_grch38_snps_wrapper <- function(chr) {
  combined <- pbmclapply(1:8,
                         function(x) get_grch38_snps(x, chr)) %>%
    rbindlist()
  outfile <- paste0(grch38_combined_dir,
                    "chr", chr, "_selscan_combined.out")
  # write combined chromosome data to file
  fwrite(combined, file = outfile,
         sep = "\t")
  # bgzip combined chromosome data
  system(command = paste0(htslib_path, "bgzip ",
                          outfile),
         intern = TRUE)
}
# apply to all chromosomes
pbmclapply(as.list(c(1:22, "X")),
           function(x) get_grch38_snps_wrapper(x))


### compare CHM13 SNP LRS values to SNPs in the same region on GRCh38

# get set of nearby GRCh38 SNPs for a given CHM13 SNP
compare_snps <- function(row_index, chm13_in, grch38_in, window) {
  # get snp of interest
  roi <- chm13_in[row_index, ]
  # get snps in GRCh38 output that are within some window of the CHM13 snp
  selscan_subset <- grch38_in[(chr == roi$chr) & (pos >= roi$end - window) & (pos <= roi$end + window)] %>%
    # add info about the CHM13 snp
    .[, chm13_ID := roi$ID] %>%
    .[, chm13_LRS := roi$lle_ratio] %>%
    .[, chm13_ancestry := roi$ancestry_component] %>%
    .[, chm13_gene := roi$gene_name]
  return(selscan_subset)
}
# wrapper function for applying to all SNPs on one chromosome
compare_snps_wrapper <- function(chr, window, selscan_in) {
  grch38_snps <- fread(paste0(grch38_combined_dir,
                              "chr", chr, "_selscan_combined.out.gz"))
  chrom <- paste0("chr", chr)
  chm13_snps <- selscan_in[chr == chrom]
  results <- pbmclapply(1:nrow(chm13_snps),
                        function(x) compare_snps(x, chm13_snps, grch38_snps, window)) %>%
    rbindlist()
  return(results)
}
# apply to all chromosomes
adjacent_snps <- pbmclapply(as.list(c(1:22, "X")),
                            function(x) compare_snps_wrapper(x, 1000, lifted)) %>%
  rbindlist()

fwrite(adjacent_snps, file = "adjacent_snps.txt",
       sep = "\t")


############################################################################

### stats

# how many lifted over SNPs have other SNPs within 1000bp on each side?
length(unique(adjacent_snps$chm13_ID)) # 2302, vs. 3038 unique liftovers

# how many adjacent SNPs are above the LRS threshold in the right ancestry component?
lrs_threshold <- 100
# how many CHM13 snps have GRCh38 snps that were significant in the right ancestry component?
length(unique(adjacent_snps[ancestry_component == chm13_ancestry & lle_ratio >= lrs_threshold]$chm13_ID)) # 345
no_adjacent_snps <- adjacent_snps[(ancestry_component != chm13_ancestry) | (lle_ratio < lrs_threshold)] %>%
  setorder(., -chm13_LRS)
unique(no_adjacent_snps$chm13_gene)

# how many CHM13 snps have a GRCh38 snp within some distance and with equal or greater LRS?
length(unique(adjacent_snps[(`lle_ratio` + 10) >= chm13_LRS & chm13_ancestry == ancestry_component]$chm13_ID))
# 782 for >= chm13_LRS; 1108 for >= chm13_LRS - 5; 1255 for  >= chm13_LRS - 10
# 1255 / 3038 * 100 # 41.3% of all lifted over variants
