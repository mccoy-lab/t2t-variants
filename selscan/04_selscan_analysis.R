library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(pbmcapply)

setwd("/scratch/groups/rmccoy22/syan11/t2t_variants/snp_selection/")

############################################################################

### DATA

# directory with selscan output
selscan_path <- "selscan_out/"
# directory with VCFs of variants
vcf_path <- "snp_genotypes/"

# path to gff3 file of CHM13 gene annotations
genes_path <- "CHM13.combined.v4.gff3.gz"
# path to files that list IDs of non-PASS variants by chromosome
nopass_path <- "nopass_IDs/"
# path to files that list IDs of variants in novel regions by chromosome
novel_ids_path <- "nonsyn_winnowmap_100mer_ids/"

# path to htslib install
htslib_path <- "/scratch4/mschatz1/rmccoy22/code/htslib-1.11/"

############################################################################


### read in and filter selscan output

# read in AF calculations, for sanity checking and filtering
afs <- pbmclapply(as.list(c(1:22, "X")),
                  function(x) fread(paste0("chr", x, ".frq.strat"),
                                    header = FALSE,
                                    col.names = c("CHR", "SNP", "CLST", "A1", "A2", "MAF", "MAC", "NCHROBS"))) %>%
  rbindlist()
# save(afs, file = "afs_all.RData")
# load("afs_all.RData")

### read in selscan output files

# get set of common, well genotyped SNPs
gtrate <- afs %>%
  group_by(SNP) %>%
  summarize(TOTOBS = sum(as.numeric(NCHROBS))) %>%
  as.data.table() %>%
  # only variants genotyped in at least 95% of samples
  .[TOTOBS > 0.95 * 2 * 2504]
afs_gtrate <- afs[SNP %in% gtrate$SNP]
common <- unique(afs_gtrate[MAF > 0.05 & MAF < 0.95]$SNP)

# read in `selscan.out` files and filter for common, well genotyped variants
filter_selscan <- function(chr, p, num, id_dt) {
  selscan <- fread(paste0(selscan_path, "chr", chr,
                          "_split", num, "_selscan_k8_p", p, ".out")) %>%
    .[, ancestry_component := p]
  selscan <- cbind(id_dt, selscan)
  setnames(selscan, "lle-ratio", "lle_ratio")
  
  # subset to variants common in any 1KGP population
  selscan_filt <- selscan[ID %in% common]
  # filter out SNPs with extreme global LLE values
  selscan_filt <- selscan_filt[`global-lle` > -900] %>%
    .[step < 99]
  return(selscan_filt)
}
# get all ancestry component selscan results for one split VCF
get_selscan_bysplit <- function(chr, num) {
  ids <- fread(cmd = paste0("zcat ", vcf_path, "chr", chr, "/",
                            "chr", chr, "_unique_IDs_split", num, ".vcf.gz",
                            " | cut -f 3 | grep -v '#'"),
               header = TRUE)
  selscan_res <- pbmclapply(1:8,
                            function(x) filter_selscan(chr, x, num, ids)) %>%
    rbindlist()
  return(selscan_res)
}
# wrapper function for applying to all chromosomes and split VCFs
get_selscan_wrapper <- function(chr) {
  split_nums <- fread(cmd = paste0("ls ", selscan_path, "chr", chr, "_split*_p1.out",
                                   " | cut -f 3 -d '_' | cut -f 2 -d 't'"),
                      colClasses = c("character"))$V1
  selscan_res <- pbmclapply(as.list(split_nums[]),
                            function(x) get_selscan_bysplit(chr, x)) %>%
    rbindlist()
  setorder(selscan_res, -lle_ratio)
  return(selscan_res)
}
selscan_all <- pbmclapply(as.list(c(1:22, "X")),
                          function(x) get_selscan_wrapper(x)) %>%
  rbindlist()

# separate fields of variant IDs to get chr and pos
selscan_all <- separate(selscan_all, ID,
                        into = c("chr", "pos", NA, NA),
                        sep = "_", remove = FALSE)
setorder(selscan_all, -lle_ratio)
# save(selscan_all, file = "selscan_all.RData")
load("selscan_all.RData")


############################################################################

### annotate selscan output files

### merge with CHM13 gene annotations

# CHM13 genes v4
genes <- fread(cmd = paste("zcat", genes_path, "| grep -w 'gene'"),
               col.names = c("chr", "source", "type", "start", "end", "score", "strand", "phase", "info"))
# parse column with gene information
genes_sep <- separate(genes, info,
                      into = c("ID", "gene_id", "gene_type", "gene_name"),
                      sep = ";", remove = TRUE) %>%
  separate(ID, into = c(NA, "ID"), sep = "=", remove = TRUE) %>%
  separate(gene_id, into = c(NA, "gene_id"), sep = "=", remove = TRUE) %>%
  separate(gene_type, into = c(NA, "gene_type"), sep = "=", remove = TRUE) %>%
  separate(gene_name, into = c(NA, "gene_name"), sep = "=", remove = TRUE) %>%
  .[, c("chr", "start", "end", "gene_id", "gene_type", "gene_name")]

# set position in selscan output to numeric, for overlapping with genes
selscan_perc_filt$pos <- as.numeric(selscan_perc_filt$pos)
selscan_perc_filt[, end := pos]
# overlap SNP positions with genes
setkey(genes_sep, chr, start, end)
selscan_genes <- foverlaps(selscan_perc_filt[, -c("gene_id", "gene_type", "gene_name")],
                           genes_sep,
                           by.x = c("chr", "pos", "end"),
                           by.y = key(genes_sep)) %>%
  .[, -c("i.end", "i.start")]


### remove non-PASS variants by ID

nopass <- lapply(as.list(c(1:22, "X")),
                 function(x) fread(paste0(nopass_path, "chr", x, "_nopass_IDs.txt"),
                                   sep = "\t", col.names = c("ID"))) %>%
  rbindlist()
selscan_filt <- selscan_genes[!(ID %in% nopass$ID)]


### calculate LRS percentiles for all variants to identify outliers

# get LRS percentile for variants in each ancestry component
get_lle_perc <- function(p) {
  subset <- selscan_filt[ancestry_component == p] %>%
    setorder(., -lle_ratio)
  subset <- subset %>%
    mutate(perc = percent_rank(lle_ratio)) %>%
    as.data.table()
  return(subset)
}
selscan_perc_filt <- pbmclapply(as.list(1:8),
                                function(x) get_lle_perc(x)) %>%
  rbindlist()

# write filtered nonsyntenic variants to file
fwrite(selscan_perc_filt, "selscan_filt_genes.txt",
       sep = "\t")
# zip file
system(command = paste0(htslib_path, "bgzip ",
                        "selscan_filt_genes.txt"))


### subset to novel snps

select_novel_snps <- function(chr_num, selscan_in) {
  chr_ids <- fread(paste0(novel_ids_path,
                          "chr", chr_num, ".nonsyn_winnowmap_100mer.annot_IDs.txt"),
                   header = FALSE)$V1
  selscan_chr <- selscan_in[chr == paste0("chr", chr_num)]
  subset <- selscan_chr[ID %in% chr_ids]
  return(subset)
}
selscan_filt_novel <- pbmclapply(as.list(c(1:22, "X")),
                                 function(x) select_novel_snps(x, selscan_perc_filt)) %>%
  rbindlist()

# write filtered novel variants to file
fwrite(selscan_filt_novel, "selscan_filt_novel_genes.txt",
       sep = "\t")
# zip file
system(command = paste0(htslib_path, "bgzip ",
                        "selscan_filt_novel_genes.txt"))


############################################################################

### stats for selection results

# number of variants in 99.9% percentile across all ancestry components
length(unique(selscan_perc_filt[perc >= 0.999]$ID)) # 5154
# number of these variants overlapping with genes
length(unique(selscan_genes[gene_name != "" & perc > 0.999]$ID)) # 989


### numbers - overlap with exons

# CHM13 genes annotated by liftoff
exons <- fread(cmd = paste("zcat", genes_path, "| grep -w 'exon'"),
               col.names = c("chr", "source", "type", "start", "end", "score", "strand", "phase", "info"))
# parse column with exon information
exons_sep <- separate(exons, info,
                      into = c("exon_id", NA, NA, "transcript_id", "gene_type",
                               "gene_name", "transcript_name", "exon_number", NA),
                      sep = ";", remove = TRUE) %>%
  separate(exon_id, into = c(NA, "exon_id"), sep = "=", remove = TRUE) %>%
  separate(transcript_id, into = c(NA, "transcript_id"), sep = "=", remove = TRUE) %>%
  separate(gene_type, into = c(NA, "gene_type"), sep = "=", remove = TRUE) %>%
  separate(gene_name, into = c(NA, "gene_name"), sep = "=", remove = TRUE) %>%
  separate(transcript_name, into = c(NA, "transcript_name"), sep = "=", remove = TRUE) %>%
  separate(exon_number, into = c(NA, "exon_number"), sep = "=", remove = TRUE) %>%
  .[, c("chr", "start", "end", "exon_id", "exon_number",
        "transcript_id", "transcript_name", "gene_type", "gene_name")]

# set position in selscan output to numeric, for overlapping with exons
selscan_perc_filt$pos <- as.numeric(selscan_perc_filt$pos)
selscan_perc_filt[, end := pos]
# overlap SNP positions with exons
setkey(exons_sep, chr, start, end)
selscan_exons <- foverlaps(selscan_perc_filt[, -c("gene_id", "gene_type", "gene_name")],
                           exons_sep,
                           by.x = c("chr", "pos", "end"),
                           by.y = key(exons_sep)) %>%
  .[, -c("i.end", "i.start")]
fwrite(selscan_exons, "selscan_exons.tsv",
       sep = "\t")

# how many of the top hits overlap with exons?
length(unique(selscan_exons[perc > 0.999 & !(is.na(exon_id))]$ID)) # 195


############################################################################

### plotting

### locus-specific manhattan plots

make_locusplot <- function(chrom, position, window, id) {
  ident <- paste0(chrom, "_", position, "_", ident)
  ident <- "chr16_37828623_A_T"
  subset <- selscan_perc_filt[chr == chrom] %>%
    .[pos > position - window & pos < position + window]
  
  # get order of ancestry components for this locus and correct plotting colors
  color_assign <- c("#f8766d", "#cc9602", "#7cae01", "#00be67",
                    "#01bfc4", "#00a9ff", "#c77cff", "#ff61cc")
  order <- sort(unique(subset[perc >= 0.999]$ancestry_component))
  color_vals <- sapply(order,
                       function(x) color_assign[x])
  
  p <- ggplot() +
    geom_point_rast(data = subset[perc < 0.999],
                    aes(x = pos, y = lle_ratio),
                    size = 0.5, color = "black") +
    geom_point_rast(data = subset[perc >= 0.999],
                    aes(x = pos, y = lle_ratio, color = as.factor(ancestry_component)),
                    size = 0.5) +
    scale_color_manual(values = color_vals) +
    theme_bw() +
    labs(title = ident,
         x = "Position", y = "Log-likelihood ratio statistic (LRS)",
         color = "Ancestry component") +
    xlim(position - window, position + window) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("locusplots/", ident, ".pdf"), p,
         width = 9, height = 5)
}
make_locusplot("chrX", 36684515, 1e6, 7)
make_locusplot("chr16", 37823759, 1e6, 5)


### manhattan plot of all data, split by ancestry component

# prep table of all snps for plotting
all_plot <- selscan_perc_filt %>%
  # use chromosome numbers for ordering
  .[, chr_num := gsub("chr", "", chr)] %>%
  .[chr == "chrX", chr_num := 23] %>%
  .[, chr_num := as.numeric(chr_num)] %>%
  # order SNPs along X axis for plotting
  setorder(., chr_num, pos) %>%
  .[, index := 1:.N, by = ancestry_component]

# label snps in novel regions based on alternating chr number
novel_ids <- unique(selscan_filt_novel$ID)
all_plot %>%
  .[, plot_id := paste0(chr_num %% 2, "_", "all")] %>%
  .[ID %in% novel_ids, plot_id := paste0(chr_num %% 2, "_", "novel")] %>%
  setorder(., -lle_ratio) %>%
  .[, gene_anc := paste(gene_name, ancestry_component, sep = "_")] %>%
  setorder(., gene_anc)
# use ID to get rid of duplicate gene annotations for 
selscan_plotting[duplicated(selscan_plotting$gene_anc), gene_name := ""]

# make chromosome ticks on the bottom of the plot
ticks <- group_by(all_plot[ancestry_component == 8], chr_num) %>%
  summarize(., median_index = median(index))

p <- ggplot() +
  geom_point_rast(data = all_plot[plot_id == "0_all" | plot_id == "1_all"],
                  aes(x = index, y = lle_ratio, color = plot_id),
                  size = 0.3) +
  geom_point_rast(data = all_plot[plot_id == "0_novel" | plot_id == "1_novel"],
                  aes(x = index, y = lle_ratio, color = plot_id),
                  size = 0.3) +
  # geom_text_repel(data = all_plot[perc > 0.999 & (plot_id == "0_novel" | plot_id == "1_novel")],
  #                 aes(x = index, y = lle_ratio, label = gene_name),
  #                 fontface = "italic", color = "black", size = 3,
  #                 force = 1.5, max.overlaps = 100000,
  #                 min.segment.length = 0) +
  theme_bw() +
  theme(legend.position = "none", panel.spacing.x = unit(0, "lines"),
        panel.border = element_blank(), panel.grid = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values = c("#FCC5C0", "#fb8072", "#C6DBEF", "#6baed6")) +
  facet_grid(ancestry_component ~ ., scales = "free") +
  ylab("Likelihood ratio statistic (LRS)") +
  scale_x_continuous(breaks = ticks$median_index, labels = c(1:22, "X"))
ggsave("selscan_ancestry.pdf", p,
       width = 8, height = 7)