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
# path to bedtools install
bedtools_path <- "~/miniconda3/bin/bedtools"

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


############################################################################

### plotting

### prep data of successful liftovers

# subset to grch38 snps with a matching ancestry component
adjacent_filt <- adjacent_snps[!(ancestry_component != chm13_ancestry)] %>%
  # add column of chm13 positions for overlapping
  separate(chm13_ID, into = c("chm13_chr", "chm13_pos", NA, NA),
           sep = "_", remove = FALSE) %>%
  setorder(., -chm13_LRS) %>%
  .[, chm13_end := chm13_pos]

# write snps to file for bedtools
adjacent_filt[, chm13_start := chm13_end - 1]
fwrite(unique(adjacent_filt[, c("chm13_chr", "chm13_start", "chm13_end")]) %>%
         setorderv(c("chm13_chr", "chm13_start", "chm13_end")),
       "adjacent_filt_snps.bed", sep = "\t", col.names = FALSE)
# use bedtools to combine adjacent SNPs into CHM13 regions for plotting
system(paste(bedtools_path, "merge",
             "-d 5000",
             "-i adjacent_filt_snps.bed",
             "> adjacent_filt_merged.bed"),
       intern = TRUE)
adjacent_filt_regions <- fread("adjacent_filt_merged.bed",
                               col.names = c("chm13_chr", "chm13_start", "chm13_end")) %>%
  # make IDs for distinguishing regions
  .[, ident := paste(chm13_chr, chm13_start, chm13_end, sep = "_")]

# annotate snps with the region they're in
setkey(adjacent_filt_regions, chm13_chr, chm13_start, chm13_end)
adjacent_filt_annot <- foverlaps(adjacent_filt, adjacent_filt_regions,
                                 by.x = c("chm13_chr", "chm13_pos", "chm13_end"),
                                 by.y = key(adjacent_filt_regions),
                                 nomatch = NULL) %>%
  .[, -c("chm13_start", "chm13_end", "i.chm13_end")]

# separate out adjacent chm13 and grch38 snps into different lines for plotting
make_adj_plot_dt <- function(id, adjacent_in, lifted_in) {
  # read in one chm13 region of interest
  subset <- adjacent_in[ident == id]
  
  # get neighboring grch38 snps for this region
  grch38_snps <- subset %>%
    .[, c("chr", "ID", "pos", "lle_ratio", "ancestry_component", "chm13_gene", "ident")] %>%
    .[chm13_gene == "", chm13_gene := NA] %>%
    .[, assembly := "grch38"]
  # get grch38 positions and LRS for chm13 query snps (from the successful liftover file)
  chm13_snps <- lifted_in[ID %in% unique(subset$chm13_ID)] %>%
    .[, c("chr", "ID", "end", "lle_ratio", "ancestry_component", "gene_name")] %>%
    .[, ident := id] %>%
    .[, assembly := "chm13"]
  
  return(rbindlist(list(grch38_snps, chm13_snps),
                   use.names = FALSE))
}
adjacent_plot <- pbmclapply(as.list(unique(adjacent_filt_annot$ident)), # 160 loci total
                            function(x) make_adj_plot_dt(x, adjacent_filt_annot, lifted)) %>%
  rbindlist()
adjacent_plot$chr <- factor(adjacent_plot$chr,
                            levels = c("chr1", "chr2", "chr3", "chr4", "chr6", "chr7",
                                       "chr9", "chr10", "chr13", "chr16", "chr17",
                                       "chr20", "chr21", "chr22"))


### plot each locus of interest

make_locusplot <- function(id, adjacent_in, selscan_in) {
  # subset to just region of interest
  subset <- adjacent_in[ident == id]
  
  # get order of ancestry components for this locus and correct plotting colors
  color_assign <- c("#f8766d", "#cc9602", "#7cae01", "#00be67",
                    "#01bfc4", "#00a9ff", "#c77cff", "#ff61cc")
  order <- sort(unique(subset$ancestry_component))
  color_vals <- sapply(order,
                       function(x) color_assign[x])
  
  ### plot with CHM13 and GRCh38 SNPs
  p1 <- ggplot() +
    # grch38 snps
    geom_point_rast(data = subset[assembly == "grch38"],
                    aes(x = pos / 1000, y = lle_ratio),
                    color = "black", size = 0.75) +
    # chm13 snps
    geom_point_rast(data = subset[assembly == "chm13"],
                    aes(x = pos / 1000, y = lle_ratio,
                        color = factor(ancestry_component,
                                       levels = order)),
                    size = 0.5) +
    scale_color_manual(values = color_vals) +
    # annotate top chm13 snp for each ancestry component (if it exists)
    geom_text_repel(data = subset[, .SD[which.max(lle_ratio)],
                                  by = ancestry_component],
                    aes(x = pos / 1000, y = lle_ratio,
                        label = chm13_gene),
                    fontface = "italic", color = "black", size = 3, max.overlaps = 30) +
    # separate facets for each chromosome in case liftover was weird
    facet_grid(~ factor(chr), scales = "free") +
    theme_bw() +
    labs(title = id,
         x = "\nPosition in chromosome on GRCh38 (kbp)",
         y = "Log-likelihood ratio statistic\n",
         color = "Ancestry component") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.spacing = unit(1, "lines"))
  
  ### plot with locus in CHM13 coordinates
  # find min and max location of chm13 snps in this locus
  chm13_pos <- lapply(as.list(unique(subset[assembly == "chm13"]$ID)),
                      function(x) as.numeric(strsplit(x, "_")[[1]][2])) %>%
    unlist()
  minpos <- min(chm13_pos) - 2500
  maxpos <- max(chm13_pos) + 2500
  selscan_subset <- selscan_in[(pos >= minpos) & (pos < maxpos)] %>%
    .[perc < 0.999, ancestry_component := NA]
  selscan_order <- sort(unique(selscan_subset$ancestry_component))
  selscan_color_vals <- sapply(selscan_order,
                               function(x) color_assign[x])
  p2 <- ggplot(data = selscan_subset,
               aes(x = pos / 1000, y = lle_ratio,
                   color = factor(ancestry_component,
                                  levels = selscan_order))) +
    geom_point_rast(size = 0.5, alpha = 0.7) +
    scale_color_manual(values = selscan_color_vals, na.value = "black") +
    # annotate top chm13 snp for each ancestry component (if it exists)
    geom_text_repel(data = selscan_subset[, .SD[which.max(lle_ratio)],
                                          by = ancestry_component],
                    aes(x = pos / 1000, y = lle_ratio,
                        label = gene_name),
                    fontface = "italic", color = "black", size = 3, max.overlaps = 30) +
    theme_bw() +
    labs(x = "\nPosition in chromosome on CHM13 (kbp)",
         y = "Log-likelihood ratio statistic\n",
         color = "Ancestry component") +
    theme(panel.spacing = unit(1, "lines"))
  
  ggsave(paste0("liftover_region_plots/", id, "_locusplot.pdf"),
         p1 / p2,
         width = 9, height = 8)
}
pbmclapply(as.list(unique(adjacent_plot$ident)),
           function(x) make_locusplot(x, adjacent_plot, selscan_perc_filt))