library(data.table)
library(magrittr)
library(tidyr)
library(pbmcapply)

############################################################################

### DATA

args <- commandArgs(trailingOnly = TRUE)

# path to agp file with clone boundaries
agp_path <- paste0("/scratch4/mschatz1/rmccoy22/rmccoy22/grch38_local_ancestry/agp_split/", args[1])
# path to 1kgp sample metadata
meta_path <- "/scratch4/mschatz1/rmccoy22/syan11/grch38_ancestry_clust/igsr_samples.tsv"
# path to tabix install
tabix_path <- "/scratch4/mschatz1/rmccoy22/code/htslib-1.11/tabix"
# path to directory with 1KGP phase 3 vcfs
kgp_path <- "/scratch4/mschatz1/rmccoy22/1kg-GRCh38-phase3/"

############################################################################


agp <- fread(agp_path) %>%
  setnames(., c("object", "object_beg", "object_end", "part_number", 
                "component_type", "component_id", 
                "component_beg", "component_end", "orientation"))

agp[, object_len := object_end - object_beg]
agp <- agp[!(component_type %in% c("O", "N")) & object_len > 1e3]

sample_metadata <- fread(meta_path) %>%
  .[, c(1, 4, 6), with = FALSE] %>%
  setnames(., c("sample_id", "pop", "superpop"))

knn_ancestry <- function(agp_file, row_index) {

  agp_subset <- agp_file[row_index]
  chrom <- agp_subset$object
  start <- agp_subset$object_beg
  end <- agp_subset$object_end
  chrom_nochr <- gsub("chr", "", chrom)

  results <- tryCatch(
    {
      if (chrom == "chrX") {
        command <- paste0(tabix_path, " -h ", kgp_path, "ALL.", chrom,
                          "_GRCh38_liftover.genotypes.20170504.vcf.gz ",
                          chrom_nochr, ":", start, "-", end)
      } else {
        command <- paste0(tabix_path, " -h ", kgp_path, "ALL.", chrom,
                          ".shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz ",
                          chrom, ":", start, "-", end)
      }

      haplotypes <- fread(cmd = command, skip = "#CHROM", fill = TRUE, colClasses = 'character') %>%
        pivot_longer(!c(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT), names_to = "sample", values_to = "gt") %>%
        as.data.table() %>%
        .[, c("h1", "h2") := tstrsplit(gt, "|", fixed = TRUE)] %>%
        .[, c("POS", "sample", "h1", "h2"), with = FALSE] %>%
        pivot_longer(!c(POS, sample), names_to = "haplotype", values_to = "gt") %>%
        as.data.table() %>%
        .[, hap := paste(sample, haplotype, sep = "_")] %>%
        .[, sample := NULL] %>%
        .[, haplotype := NULL] %>%
        pivot_wider(names_from = hap, values_from = gt) %>%
        as.data.table()

      haplotypes <- haplotypes[, colSums(is.na(haplotypes)) < nrow(haplotypes), with = FALSE]

      hap_mat <- as.matrix(t(haplotypes[, -1]))
      class(hap_mat) <- "numeric"
      #set.seed(1)
      #hap_mat <- as.matrix(hap_mat[sample(1:nrow(hap_mat), replace = FALSE),])
      hap_mat <- rbind(GRCh38 = rep(0, ncol(hap_mat)), hap_mat)
      hap_dists <- dist(hap_mat, method = "manhattan")

      min_dist <- min(as.matrix(hap_dists)[, 1][-1])
      nn <- names(as.matrix(hap_dists)[, 1][-1][which(as.matrix(hap_dists)[, 1][-1] == min_dist)])
      nn <- gsub("_h1", "", nn)
      nn <- gsub("_h2", "", nn)

      sample_metadata$superpop <- factor(sample_metadata$superpop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS"))
      superpop_count <- table(unlist(sample_metadata[match(nn, sample_metadata$sample_id),]$superpop))
      nsnp <- nrow(haplotypes)
      results <- data.table(chrom, start, end,
                            nsnp = nsnp,
                            min_dist = min_dist,
                            component_id = agp_subset$component_id,
                            AFR = superpop_count[1],
                            AMR = superpop_count[2],
                            EAS = superpop_count[3],
                            EUR = superpop_count[4],
                            SAS = superpop_count[5])
    },
    error=function(cond) {
      results <- data.table(chrom, start, end, nsnp = NA, component_id = agp_subset$component_id,
                            AFR = NA,
                            AMR = NA,
                            EAS = NA,
                            EUR = NA,
                            SAS = NA)
      return(results)
    },
    warning=function(cond) {
      results <- data.table(chrom, start, end, nsnp = NA, component_id = agp_subset$component_id,
                            AFR = NA,
                            AMR = NA,
                            EAS = NA,
                            EUR = NA,
                            SAS = NA)
      return(results)
    },
    finally={
      message(paste("Processed clone:", agp_subset$component_id))
    }
  )
  return(results)
}

output <- rbindlist(pbmclapply(1:nrow(agp), function(x) knn_ancestry(agp, x), mc.cores = 12), fill = TRUE)
fwrite(output, file = paste0("/scratch4/mschatz1/rmccoy22/rmccoy22/grch38_local_ancestry/", args[1], ".out"), 
       sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
