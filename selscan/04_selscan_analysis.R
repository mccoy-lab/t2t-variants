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