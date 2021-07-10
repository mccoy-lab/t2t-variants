library(data.table)
library(dplyr)

setwd("/scratch/groups/rmccoy22/syan11/t2t_variants/ibdmix/")

############################################################################

### DATA

# path to ibdmix output files, split by chromosome
ibdmix_path <- "/scratch4/mschatz1/rmccoy22/archaic_hominin/ibdmix/"

# local ancestry calls for grch38
ancestry_file <- "grch38_local_ancestry.txt"

# path to ucsc liftover tool
liftover_path <- "~/code/liftOver"
# liftover chain file grch38 --> chm13
hg38_to_chm13_chain <- "/scratch/groups/rmccoy22/syan11/t2t_variants/hg38.t2t-chm13-v1.0.over.chain"

############################################################################


# read in ibdmix output for CHM13
read_ibdmix_output <- function(chr_id) {
  gts <- fread(paste0(input_path, "chr", chr_id, "_vindija_GRCh38.out")) %>%
    # haplotypes must be 1Mb
    .[end - start > 50000] %>%
    # lod score must be > 4
    .[slod > 4] %>%
    # subset to CHM13 introgressed haplotypes only
    .[ID == "CHM13"] %>%
    .[, -c("ID")] %>%
    .[, chrom := paste0("chr", chrom)]
  return(gts)
}
chm13_neand <- lapply(as.list(c(1:22, "X")),
                      function(x) read_ibdmix_output(x, pop)) %>%
  rbindlist()


############################################################################

### compare to introgression results from grch38

# grch38 local ancestry for each clone
ancestry <- fread(ancestry_file) %>%
  # filter for confident ancestry regions
  .[pmax > 0.5]

### lift over clone boundaries to CHM13

# get rid of spaces in CloneName field because liftOver makes it a separate column
ancestry$CloneName <- gsub(" ", "", ancestry$CloneName)
# write confident clones to bed-like file for liftover
fwrite(ancestry, "grch38_good_clones.bed",
       col.names = FALSE, sep = "\t")
# liftover high confidence clone boundaries to CHM13 coords
system(paste(liftover_path,
             "grch38_good_clones.bed",
             hg38_to_chm13_chain,
             "grch38_local_ancestry_CHM13lifted.txt",
             "grch38_good_clones_unmapped.bed",
             "-bedPlus=3"),
       intern = TRUE)
# 588 clones unmapped - how much sequence is that?
fread(cmd = ("grep -v '#' grch38_good_clones_unmapped.bed"),
      col.names = names(ancestry)) %>%
  summarize(sum(length)) # 31489952 bp, or 31.5 Mb (out of 2348 Mb mapped, or 1.3%)

### overlap lifted over clone boundaries with CHM13 introgressed haplotypes

# read in lifted over clones
ancestry_lifted <- fread(cmd = "cut -f 1-9 grch38_local_ancestry_CHM13lifted.txt",
                         col.names = names(ancestry))
# find only the exactly overlapping regions in the CHM13 ibdmix output
setkey(ancestry_lifted, chrom, start, end)
overlap <- foverlaps(chm13_neand, ancestry_lifted,
                     by.x = c("chrom", "start", "end"),
                     by.y = key(ancestry_lifted),
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

# how much total introgressed sequence is in these regions in CHM13?
sum(overlap$length) / 1e6 # 43.56194 Mb
# how much total introgressed sequence without subsetting to these regions?
sum(chm13_neand$end - chm13_neand$start) / 1e6 # 50.43486 Mb