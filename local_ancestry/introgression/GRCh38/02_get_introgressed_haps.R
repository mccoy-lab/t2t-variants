library(data.table)
library(dplyr)

setwd("/scratch/groups/rmccoy22/syan11/t2t_variants/ibdmix/grch38_introgression")

############################################################################

### DATA

# directory with ibdmix output files, naming format: `chr{}_{superpop}_GRCh38_ibdmix.out`
input_path <- "ibdmix_output/"
# local ancestry calls for grch38
ancestry_file <- "grch38_local_ancestry.txt"

############################################################################


### subsetting introgression calls to only regions of matched ancestry

# read in grch38 introgression calls from each superpopulation
read_ibdmix_output <- function(chr_id, pop) {
  gts <- fread(paste0(input_path, "chr", chr_id, "_", pop, "_GRCh38_ibdmix.out")) %>%
    # haplotypes must be 1Mb
    .[end - start > 50000] %>%
    # lod score must be > 4
    .[slod > 4] %>%
    .[ID == "GRCh38"] %>%
    .[, superpop := pop] %>%
    .[, -c("ID")] %>%
    .[, chrom := paste0("chr", chrom)]
  return(gts)
}
read_ibdmix_wrapper <- function(pop) {
  neand <- lapply(as.list(c(1:22, "X")),
                  function(x) read_ibdmix_output(x, pop)) %>%
    rbindlist()
}
grch38_neand_all <- lapply(as.list(c("EAS", "SAS", "EUR", "AMR", "AFR")),
                           function(x) read_ibdmix_wrapper(x)) %>%
  rbindlist()

# read in local ancestry
ancestry <- fread(ancestry_file) %>%
  # filter for confident ancestry regions
  .[pmax > 0.5]
# only used introgression calls in regions of matched ancestry
match_ancestry <- function(pop) {
  neand_sub <- grch38_neand_all[superpop == pop]
  ancestry_sub <- ancestry[ancestry == pop]
  setkey(ancestry_sub, chrom, start, end)
  overlap <- foverlaps(neand_sub, ancestry_sub,
                       by.x = c("chrom", "start", "end"),
                       by.y = key(ancestry_sub),
                       nomatch = NULL)
  # get start and end positions of overlapped region
  overlap %>%
    .[, olap_start := unlist(lapply(1:nrow(overlap),
                                    function(x) max(overlap[x, start],
                                                    overlap[x, i.start])))] %>%
    .[, olap_end := unlist(lapply(1:nrow(overlap),
                                  function(x) min(overlap[x, end],
                                                  overlap[x, i.end])))] %>%
    .[, length := olap_end - olap_start]
  return(overlap)
}
neand_ancestry_match <- lapply(as.list(c("EAS", "SAS", "EUR", "AMR", "AFR")),
                               function(x) match_ancestry(x)) %>%
  rbindlist()


############################################################################

### stats

# grch38

# count how much introgressed sequence there is per superpop
neand_anc_count <- neand_ancestry_match %>%
  group_by(superpop) %>%
  summarize(neand_total = sum(length),
            slod_avg = mean(slod)) %>%
  as.data.table()
# count how much introgressed sequence there is in total
sum(neand_anc_count$neand_total) / 1e6 # 26.68476

# how much of the total sequence for each superpop is introgressed?
anc_count <- ancestry %>%
  .[, length := end - start] %>%
  group_by(ancestry) %>%
  summarize(total = sum(length)) %>%
  as.data.table()
neand_anc_count <- cbind(neand_anc_count, anc_count[, c("total")]) %>%
  .[, neand_prop := neand_total/total]