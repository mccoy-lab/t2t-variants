library(data.table)
library(dplyr)
library(pbmcapply)

### USAGE: Rscript calculate_ld.R <population ident>
### <population ident>: identification for set of populations to calculate LD within
###
### For every lifted over GWAS SNP, generate the SNP ID that would be in the 1KGP VCF
### and return SNPs with LD in the 1KGP samples.

setwd("/scratch4/mschatz1/rmccoy22/syan11/gwas_ld")

############################################################################

### DATA

# path to directory of GATK vcfs with snps of interest
gwas_vcf_path <- "gwas_hits_bypop/"

# path to directory of 1KGP VCFs called on CHM13
kgp_vcf_path <- "1kg-CHM13-recalibrated-PASS-annot/"
# path to directory of VCFs with 1KGP novel SNPs
novel_kgp_path <- "/scratch4/mschatz1/rmccoy22/syan11/gwas_ld/nonsyn_vcfs/"

# path to list of samples in 1KGP novel vcfs
novel_samples_path <- "nonsyn_samples.txt"

# path to htslib install
htslib_path <- "/scratch4/mschatz1/rmccoy22/code/htslib-1.11/"
# path to bcftools install
bcftools_path <- "/scratch4/mschatz1/rmccoy22/code/bcftools-1.11/bcftools"
# path to plink install
plink_path <- "/scratch4/mschatz1/rmccoy22/code/plink"

############################################################################


# read in population identifier as command-line argument
args <- commandArgs(trailingOnly = TRUE)
ident <- args[1]
# vcf of lifted over gwas variants
vcf_file <- paste0(ident, "_snps.vcf.gz")
# file of 1KGP samples for plink to calculate LD between
samples_path <- paste0(gwas_vcf_path, ident, "_samples.txt")

# calculate ld between gwas variant and all snps within some window
get_ld <- function(index, window_size, r2_thresh, vcf_in) {
  # get row of interest
  subset <- vcf_in[index, ]
  chr <- subset$`#CHROM`
  pos <- as.numeric(subset$POS)
  
  # generate SNP ID that should be in samantha's 1KGP VCF
  snp_id <- paste(chr, pos, subset$REF, subset$ALT, sep = "_")
  #print(snp_id)
  
  # get interval to select with tabix
  wide_interval <- paste0(chr, ":", pos - window_size/2, "-", pos + window_size/2)
  # get path to vcf of novel 1KGP variants
  chr_novel_kgp_path <- paste0(novel_kgp_path, chr,
			       ".1Mb_nonsyn.annot.vcf.gz")
  # use tabix to select window around SNP of interest
  row_count <- as.numeric(system(command = paste(paste0(htslib_path, "tabix"),
                                                 chr_novel_kgp_path,
                                                 wide_interval,
                                                 "| wc -l"),
                                 intern = TRUE))
  
  # if there are novel SNPs within the window, run the rest of the function
  if (row_count > 0) {
    # use tabix to select window around SNP of interest
    window_vcf_out <- paste0("ld_output/", snp_id, ".", ident, "_window.vcf")
    system(command = paste(paste0(htslib_path, "tabix"), "-h",
                           chr_novel_kgp_path,
                           wide_interval,
                           ">", window_vcf_out),
           intern = TRUE)
    
    # get interval of just the SNP to select with tabix
    snp_interval <- paste0(chr, ":", pos, "-", pos)
    # get path to VCF of 1KGP variants
    chr_kgp_path <- paste0(kgp_vcf_path, "1kgp.", chr, ".recalibrated.snp_indel.pass.annot.vcf.gz")
    # use tabix to select SNP of interest from main 1KGP VCF
    snp_vcf_out <- paste0("ld_output/", snp_id, ".", ident, "_SNP.vcf")
    system(command = paste(paste0(htslib_path, "tabix"), "-h",
                           chr_kgp_path,
                           snp_interval,
                           # subset to the same samples as the novel regions vcf
                           "|", bcftools_path, "view",
                           "-S", novel_samples_path,
                           ">", snp_vcf_out),
           intern = TRUE)
    
    # combine snp-specific and interval vcfs for calculating LD
    concat_vcf_out <- paste0("ld_output/", snp_id, ".", ident, ".vcf.gz")
    system(command = paste(bcftools_path, "concat",
                           snp_vcf_out, window_vcf_out,
                           "|", bcftools_path, "sort",
                           "-O z",
                           "-o", concat_vcf_out),
           intern = TRUE)
    
    # calculate LD between snp of interest within the window, using plink
    ld_out <- paste0("ld_output/", snp_id, ".", ident)
    system(command = paste(plink_path,
                           "--vcf", concat_vcf_out,
                           "--keep", samples_path,
                           "--r2",
                           "--ld-snp", snp_id,
                           "--ld-window", window_size,
			   #"--ld-window-kb", window_size,
                           "--ld-window-r2", r2_thresh,
                           "--out", ld_out),
           intern = TRUE, ignore.stderr = TRUE)
    # if ld output file exists, read it in and return it
    if (file.exists(paste0(ld_out, ".ld"))) {
      ld_output <- fread(paste0(ld_out, ".ld"))
      return(ld_output)
    }
    
    # remove intermediate files
    #system(command = paste("rm",
                           #window_vcf_out, snp_vcf_out,
                           #concat_vcf_out,
                           #paste0(ld_out, ".log"), paste0(ld_out, ".nosex")))
  }
}

# read in gwas vcf
gwas_vcf <- fread(cmd = paste0("zcat ", gwas_vcf_path, vcf_file,
                               " | grep -v '##'"))

#ld_combined <- pbmclapply(1:5,
ld_combined <- pbmclapply(1:nrow(gwas_vcf),
                          function(x) get_ld(x, 1e6, 0.5, gwas_vcf)) %>%
  rbindlist()
fwrite(ld_combined,
       file = paste0("ld_output/", ident, ".novel_snp_ld.tsv"),
       sep = "\t")
