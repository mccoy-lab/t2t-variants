library(data.table)
library(dplyr)
library(ggplot2)

### Plot local ancestry inferred by RFMix for all 23 chromosomes of CHM13.

setwd("/Users/syan/Documents/mccoy-lab/t2t_variants/local_ancestry/chm13")

############################################################################

### DATA

# path to rfmix output files
rfmix_path <- "rfmix_output/phase3/"
# path to ibdmix introgression tracts
ibdmix_path <- "chm13_neand_GRCh38coords.bed"

# path to bed file of annotated centromere regions for each chromosome
cent_path <- "../centromeres.tsv"
# path to table of chromosome lengths
chrom_path <- "../chromosome_lengths.tsv"
# path to table of 1KGP masked regions
mask_path <- "../1kgp_masked_grch38.bed"

############################################################################


### Reading in data

get_rfmix_results <- function(chr_num) {
  chrID <- paste0("chr", chr_num)
  
  if (chrID == "chrX") {
    rfmix <- fread(paste0(rfmix_path, chrID, "_female.msp.tsv"),
                   skip = 2,
                   col.names = c("#chm","spos","epos","sgpos","egpos","n_snps","syndip.0","syndip.1"))
  } else {
    rfmix <- fread(paste0(rfmix_path, chrID, ".msp.tsv"),
                   skip = 2,
                   col.names = c("#chm","spos","epos","sgpos","egpos","n_snps","syndip.0","syndip.1"))
  }
  
  # scale positions to kb
  rfmix[, spos := spos/1000]
  rfmix[, epos := epos/1000]
  
  # add column with chromosomes as numbers, for plotting in order
  rfmix$chr_num <- gsub("chr", "", rfmix$`#chm`)
  rfmix[chr_num == "X", chr_num := 23]
  
  # convert ancestry column from numbers to population names
  ancestry_assignments <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  rfmix[, syndip.0 := ancestry_assignments[(syndip.0) + 1]]
  
  return(rfmix)
}
# get rfmix results for all chromosomes
rfmix_res <- lapply(as.list(c(1:22, "X")),
                    function(x) get_rfmix_results(x)) %>%
  rbindlist()

# get introgression regions
ibdmix <- fread(ibdmix_path) %>%
  # scale positions to kb
  .[, V2 := V2 / 1000] %>%
  .[, V3 := V3 / 1000]
# convert chromosome names to numbers for plotting
ibdmix$V1 <- gsub("chr", "", ibdmix$V1)
ibdmix[V1 == "X", V1 := 23]

# get table of chromosome lengths
chrom_lengths <- fread(chrom_path) %>%
  # get rid of Y and M chromosomes, since they're not in dataset
  .[!(chr_name %in% c("chrY", "chrM")), ] %>%
  # scale positions to kb
  .[, length := length/1000]

# get table of centromere positions for chromosomes
centromeres <- fread(cent_path) %>%
  # get rid of Y chromosome, since it's not in dataset
  .[chrom != "chrY", ] %>%
  # get smallest start and largest end position for centromere annotations
  # so that there's only one centromere annotation for chromosome
  group_by(chrom) %>%
  summarize(start = min(chromStart), end = max(chromEnd)) %>%
  setDT()
# convert chromosome names to numbers for plotting
centromeres$chrom <- gsub("chr", "", centromeres$chrom)
centromeres[chrom == "X", chrom := 23]
# scale positions to kb
centromeres[, start := start/1000]
centromeres[, end := end/1000]

# get table of 1KGP masked regions
masked <- fread(mask_path) %>%
  # scale positions to kb
  .[, V2 := V2 / 1000] %>%
  .[, V3 := V3 / 1000] %>%
  # get rid of extra chromosomes
  .[nchar(V1) < 6] %>%
  .[V1 != "chrY" & V1 != "chrM"]
# manually add in biggest gaps in the chrX vcf, which includes the PARS
masked <- rbind(masked,
                data.frame(V1 = "chrX",
                           V2 = 2781457 / 1000,
                           V3 = (2781457 + 152922355) / 1000,
                           V4 = "PAR_gap"))
# convert chromosome names to numbers for plotting
masked$V1 <- gsub("chr", "", masked$V1)
masked[V1 == "X", V1 := 23]


############################################################################

### Plotting

# width of one chromosome bar
width <- 0.5
# data in order of appearance on "y" axis of plot
order <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
           "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
           "chr22", "chrX")

p <- ggplot(rfmix_res) +
  
  ### adding data
  # plot ancestry blocks for each chromosome
  geom_rect(aes(xmin = as.numeric(chr_num) - width/2,
                xmax = as.numeric(chr_num) + width/2,
                ymin = spos, ymax = epos,
                fill = syndip.0)) +
  scale_fill_manual(name = "Inferred\nancestry",
                    values = c("#67c2a5", "#fc8d62", "#e68ac3", "#8d9fcb", "#a6d853")) +
                    # values = c("#44AA99", "#CC6677", "#332288", "#88CCEE", "#DDCC77")) +
  # plot introgressed regions
  geom_rect(data = ibdmix,
            aes(xmin = as.numeric(V1) - width/2,
                xmax = as.numeric(V1) + width/2,
                ymin = V2, ymax = V3),
            fill = "#ffd92f") +
  # plot masked regions and gaps
  geom_rect(data = masked,
            aes(xmin = as.numeric(V1) - width/2,
                xmax = as.numeric(V1) + width/2,
                ymin = V2, ymax = V3),
            fill = "#d3d3d3") +
  # plot centromere positions
  geom_rect(data = centromeres,
            aes(xmin = as.numeric(chrom) - width/2,
                xmax = as.numeric(chrom) + width/2,
                ymin = start, ymax = end),
            fill = "white", color = "black", size = 0.3) +
  # plot outlined chromosome bar with full chromosome lengths
  geom_rect(data = chrom_lengths,
            aes(xmin = as.numeric(index) - width/2,
                xmax = as.numeric(index) + width/2,
                ymin = 0, ymax = length),
            fill = NA, color = "black", size = 0.4) +
  
  ### plot reformatting
  # flip x and y axes
  coord_flip() +
  theme_classic() +
  # change "y" axis ticks to discrete scale so each chromosome gets its own tick
  scale_x_discrete(limits = order,
                   # spacing between "y" axis and start of chromosome bars
                   expand = c(0, width*1.5)) +
  # spacing between "x" axis and bottom of chr1 bar
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Position (Kbp)", x = "Chromosome") +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  theme(legend.position = c(0.85, 0.7),
        legend.key.size = unit(0.8, 'cm'),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))

ggsave("chm13_ancestry_phase3data_masked.pdf", p,
       width = 8, height = 6)
