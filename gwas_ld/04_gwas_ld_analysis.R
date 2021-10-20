library(data.table)
library(dplyr)

setwd("/scratch/groups/rmccoy22/syan11/t2t_variants/gwas_ld")

############################################################################

### DATA

# directory with ld output files by population
ld_path <- "ld_output/"

# path to vcf of subset gwas variants on CHM13
vcfpath <- "gwas_hits_bypop/"
# path to subset of gwas variants used for LD analysis
gwas_subset_path <- "gwas_subset.txt"
# path to gwas catalog file
gwas_path <- "gwas_all_associations_v1.0.2.tsv"

# path to htslib install
htslib_path <- "/scratch4/mschatz1/rmccoy22/code/htslib-1.11/"

############################################################################


### read in ld output files

file_list <- as.list(c("ACB.ASW.novel_snp_ld.tsv", "CHS.CHB.JPT.novel_snp_ld.tsv",
                       "FIN.TSI.IBS.CEU.GBR.novel_snp_ld.tsv",
                       "GWD.MSL.ESN.YRI.LWK.novel_snp_ld.tsv",
                       "PUR.CLM.PEL.MXL.novel_snp_ld.tsv"))
ld <- rbindlist(lapply(file_list,
                       function(x) {fread(paste0(ld_path, x)) %>%
                           .[, pop := x]})) %>%
  # remove SNPs that are in LD with themselves
  .[SNP_A != SNP_B] %>%
  .[, pop := gsub(".novel_snp_ld.tsv", "", pop)]

# get rsIDs because variant names are in `chr_pos_ref_alt` format
vcf_list <- as.list(c("ACB.ASW_snps.vcf.gz", "CHS.CHB.JPT_snps.vcf.gz",
                      "FIN.TSI.IBS.CEU.GBR_snps.vcf.gz",
                      "GWD.MSL.ESN.YRI.LWK_snps.vcf.gz",
                      "PUR.CLM.PEL.MXL_snps.vcf.gz"))
vcf <- rbindlist(lapply(vcf_list,
                        function(x) {fread(cmd = paste0("zcat ",vcfpath, x,
                                                        " | grep -v '##'")) %>%
                            .[, pop := x]})) %>%
  .[, pop := gsub("_snps.vcf.gz", "", pop)] %>%
  .[, c(1:5, 9)] %>%
  .[, alt_list := strsplit(ALT, ",")] %>%
  unnest(., alt_list) %>%
  as.data.table() %>%
  .[, SNP_A := paste(`#CHROM`, POS, REF, alt_list, sep = "_")]
ld_merged <- merge(ld, vcf[, c("SNP_A", "pop", "ID")], by = c("SNP_A", "pop"))


### merge with metadata from gwas catalog files

# subset of gwas variants used for ld analysis
gwas <- fread(gwas_subset_path) %>%
  setnames(., "SNPS", "ID")
# all variants in gwas catalog
all_assoc <- fread(gwas_path) %>%
  .[, c("PUBMEDID", "DISEASE/TRAIT")] %>%
  .[!duplicated(PUBMEDID)]
# merge with full gwas catalog to get pvalues
gwas <- merge(gwas, all_assoc, by = "PUBMEDID")
ld_merged <- merge(ld_merged, gwas, by = "ID") %>%
  .[, -"pop"] %>%
  .[, `P-VALUE` := as.numeric(`P-VALUE`)] %>%
  setorder(., `P-VALUE`)

# write to file
fwrite(ld_merged, "ld_output.txt",
       sep = "\t")
# zip file
system(command = paste0(htslib_path, "bgzip ",
                       "ld_output.txt"))

# write supplementary table of unique gwas hits in LD
ld_merged_supp <- ld_merged %>%
  setorder(., ID, -R2) %>%
  # remove rows with duplicate rsIDs, keeping only the first that appears
  # (the highest r^2 value between the gwas snp and a nonsyntenic snp)
  .[!duplicated(ID), ]
setnames(ld_merged_supp,
         c("ID", "SNP_A",  "SNP_B", "PUBMEDID",
           "P-VALUE", "broad_ancestry", "DISEASE/TRAIT"),
         c("snp_id", "chm13_id", "nonsyntenic_snp_in_strongest_ld", "pubmed_study_id",
           "pvalue", "gwas_catalog_ancestry", "phenotype"))
fwrite(ld_merged_supp[, c("snp_id", "chm13_id", "nonsyntenic_snp_in_strongest_ld",
                          "R2", "pubmed_study_id", "gwas_pvalue",
                          "gwas_catalog_ancestry", "kgp_ancestry", "phenotype")],
       "gwas_ld_results.csv")


### stats

# how many unique snps with LD > 0.5 in a novel region?
length(unique(ld$SNP_A)) # 113
# we started with 22474 unique variants - calculated with
# `cut -f 4 gwas_subset.txt | sort | uniq | wc -l`
# calculate percentage of starting variants that 
113 / 22474 * 100
