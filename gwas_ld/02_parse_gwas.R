library(data.table)
library(dplyr)
library(pbmcapply)

### Subset GWAS Catalog to variants in studies that correspond broadly to
### populations in the 1000 Genomes Project.

setwd("/scratch4/mschatz1/rmccoy22/syan11/gwas_ld")

############################################################################

### DATA

# path to ancestry file for gwas catalog
ancestry_path <- "gwas_catalog-ancestry_r2021-05-19.tsv"
# path to gwas catalog file
gwas_path <- "gwas_all_associations_v1.0.2.tsv"

# path to gatk install
gatk_path <- "/scratch4/mschatz1/rmccoy22/code/gatk-4.1.9.0"
# path to lifted over dbSNP vcf
dbsnp_path <- "/scratch4/mschatz1/rmccoy22/syan11/gwas_ld/dbsnp_lifted/GCF_000001405.CHM13.all.lifted.vcf.gz"
# path to 1KGP within file with population assignments
kgp_meta_path <- "1KGP_within.txt"

############################################################################


### parse gwas catalog and ancestry metadata to generate lists of SNPs for each population

ancestry <- fread(ancestry_path) %>%
  .[, c("PUBMEDID", "BROAD ANCESTRAL CATEGORY", "ADDITONAL ANCESTRY DESCRIPTION")] %>%
  setnames(c("pubmedid", "broad_ancestry", "additional_ancestry"))
# remove studies that have multiple ancestry annotations because then we can't calculate LD
ancestry <- ancestry[ancestry$pubmedid %in% names(which(table(ancestry$pubmedid) == 1))]
gwas <- fread(gwas_path)

# get most common ancestry categories to convert into 1KGP populations
ancestry_table <- as.data.frame(table(ancestry$broad_ancestry)) %>%
  setorder(., -Freq) %>%
  setDT()
# look at additional categories to figure out what that ancestry entails
table(ancestry[broad_ancestry == "Native American"]$additional_ancestry)
# assign population IDs to ancestry names
ancestry[broad_ancestry == "European", kgp_ancestry := "FIN TSI IBS CEU GBR"] %>%
  .[broad_ancestry == "East Asian", kgp_ancestry := "CHS CHB JPT"] %>%
  .[broad_ancestry == "African American or Afro-Caribbean", kgp_ancestry := "ACB ASW"] %>%
  .[broad_ancestry == "South Asian", kgp_ancestry := "BEB PJL GIH ITU STU"] %>%
  .[broad_ancestry == "Hispanic or Latin American", kgp_ancestry := "PUR CLM PEL MXL"] %>%
  .[broad_ancestry == "Hispanic or Latin American, Native American", kgp_ancestry := "PUR CLM PEL MXL"] %>%
  .[broad_ancestry == "Native American", kgp_ancestry := "PUR CLM PEL MXL"] %>%
  .[broad_ancestry == "Sub-Saharan African", kgp_ancestry := "GWD MSL ESN YRI LWK"]
# get only studies with easily interpretable ancestries
ancestry_annot <- ancestry[!(is.na(kgp_ancestry))]

# subset gwas catalog to easily interpretable ancestries
gwas_subset <- merge(gwas[, c("PUBMEDID", "CHR_ID", "CHR_POS", "SNPS", "P-VALUE")],
                     ancestry_annot[, c("pubmedid", "broad_ancestry", "kgp_ancestry")],
                     by.x = "PUBMEDID", by.y = "pubmedid")
# get rid of weird chromosome annotations
gwas_subset <- gwas_subset[!grepl(";", gwas_subset$CHR_ID) & !grepl("x", gwas_subset$CHR_ID)]
fwrite(gwas_subset, "gwas_subset.txt", sep = "\t")


############################################################################

### write files for subsetting SNPs and samples

# read in 1000 genomes metadata
kgp_meta <- fread(kgp_meta_path, col.names = c("id", "id2", "pop"))

select_variants <- function(population) {
  subset <- gwas_subset[kgp_ancestry == population]
  pop_id <- gsub(" ", ".", subset$kgp_ancestry[1])
  
  # write snp IDs for the ancestry to file, for selecting variants with GATK
  fwrite(subset[, c("SNPS")], paste0(pop_id, "_snps.list"),
         col.names = FALSE)
  # run gatk
  system(paste(gatk_path, "SelectVariants",
               "-V", dbsnp_path,
               "-O", paste0(pop_id, "_snps.vcf.gz"),
               "--keep-ids", pop_id, "_snps.list"))
  
  # generate list of samples for plink to calculate LD in
  samples <- kgp_meta[pop %in% unlist(strsplit(population, split = " "))]
  fwrite(samples[, c("id", "id2")], paste0(pop_id, "_samples.txt"),
         col.names = FALSE, sep = "\t")
}
pbmclapply(as.list(unique(gwas_subset$kgp_ancestry)),
           function(x) select_variants(x))
