library(data.table)
library(tidyverse)
library(cowplot)

sample_list <- fread("~/Downloads/sample_list.txt", header = FALSE) %>%
  setnames(., "sample_id")
igsr <- fread("~/Downloads/igsr_samples.tsv")

# autosomes
table(igsr[`Sample name` %in% sample_list$sample_id]$`Superpopulation code`) * 2

# X
table(igsr[`Sample name` %in% sample_list$sample_id & Sex == "female"]$`Superpopulation code`) * 2 + 
  table(igsr[`Sample name` %in% sample_list$sample_id & Sex == "male"]$`Superpopulation code`)

loc <- fread("~/Downloads/grch38_clust.txt", fill = TRUE) %>%
  setnames(., "component_id", "accession")

lib <- fread("https://ftp.ncbi.nih.gov/repository/clone/reports/Homo_sapiens/clone_acstate_9606.out") %>%
  setnames(., "Accession", "accession")

loc <- merge(loc, lib)
loc[, clone_id := paste(chrom, start, end, sep = "_")]
loc <- loc[!duplicated(clone_id)]
loc[, nnn := AFR + AMR + EAS + EUR + SAS]
loc[, AFR := as.numeric(AFR)]
loc[, AMR := as.numeric(AMR)]
loc[, EAS := as.numeric(EAS)]
loc[, EUR := as.numeric(EUR)]
loc[, SAS := as.numeric(SAS)]
loc[chrom != "chrX", AFR := AFR / 1342]
loc[chrom != "chrX", AMR := AMR / 696]
loc[chrom != "chrX", EAS := EAS / 1030]
loc[chrom != "chrX", EUR := EUR / 1044]
loc[chrom != "chrX", SAS := SAS / 984]
loc[chrom == "chrX", AFR := AFR / 1019]
loc[chrom == "chrX", AMR := AMR / 526]
loc[chrom == "chrX", EAS := EAS / 782]
loc[chrom == "chrX", EUR := EUR / 795]
loc[chrom == "chrX", SAS := SAS / 722]
loc[, sum := AFR + AMR + EAS + EUR + SAS]
loc[, AFR := AFR / sum]
loc[, AMR := AMR / sum]
loc[, EAS := EAS / sum]
loc[, EUR := EUR / sum]
loc[, SAS := SAS / sum]

set.seed(1)
plurality <- colnames(loc[, c("AFR", "AMR", "EAS", "EUR", "SAS")])[max.col(loc[, c("AFR", "AMR", "EAS", "EUR", "SAS")], ties.method = "random")]
loc[, ancestry := plurality]
loc[, pmax := pmax(AFR, AMR, EAS, EUR, SAS)]

loc[, len := end - start]
loc$chrom <- factor(loc$chrom, levels = paste0("chr", c(1:22, "X")))

a <- ggplot(data = loc[!is.na(pmax)], aes(x = nsnp, fill = pmax > 0.5)) + 
  geom_histogram(position = "dodge") + 
  xlab("Number of SNPs") + 
  ylab("Number of clones") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_brewer(palette = "Set2")

b <- ggplot(data = loc[!is.na(pmax)], aes(x = len, fill = pmax > 0.5)) + 
  geom_histogram(position = "dodge") + 
  xlab("Clone length (bp)") + 
  ylab("Number of clones") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_fill_brewer(palette = "Set2")

plot_grid(a, b)

ggplot(data = loc[pmax > 0.5], aes(x = LibAbbr, y = stat(count), fill = ancestry)) +
  geom_bar() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_brewer(palette = "Set2") +
  ylab("Number of clones") +
  xlab("Library")

clone_prop <- group_by(loc, LibAbbr) %>%
  summarize(., tot_len = sum(len)) %>%
  as.data.table() %>%
  setorder(., -tot_len)

clone_prop_abbrev <- rbind(clone_prop[1:7],
      data.table(LibAbbr = "Other", tot_len = sum(clone_prop[8:nrow(clone_prop)]$tot_len)))

clone_prop_abbrev$LibAbbr <- factor(clone_prop_abbrev$LibAbbr, levels = clone_prop_abbrev$LibAbbr)

my_palette <- c("#d885a0",
                "#8ebf7c",
                "#c79ad6",
                "#c4b06f",
                "#799ad6",
                "#d58e72",
                "#6ccab4",
                "#6bbbd3")

a <- ggplot(data = clone_prop_abbrev, aes(x = "", y = tot_len, fill = LibAbbr)) + 
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme(legend.position = "none") +
  scale_fill_manual(values = my_palette, name = "Library") +
  theme_minimal() +
  theme(panel.grid = element_blank(), panel.background = element_blank(), axis.text = element_blank()) +
  xlab("") +
  ylab("")

loc_abbrev_other <- loc[!(LibAbbr %in% clone_prop_abbrev$LibAbbr)]
loc_abbrev_other[, LibAbbr := "Other"]
loc_abbrev <- rbind(loc[LibAbbr %in% clone_prop_abbrev$LibAbbr], loc_abbrev_other)
loc_abbrev$LibAbbr <- factor(loc_abbrev$LibAbbr, levels = clone_prop_abbrev$LibAbbr)

b <- ggplot(data = loc_abbrev[!is.na(ancestry)], 
       aes(x = LibAbbr, y = stat(count), fill = ancestry)) +
  geom_bar(position = "fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid = element_blank()) +
  scale_fill_manual(values = c('#66c2a5','#fc8d62','#e78ac3', '#8da0cb','#a6d854'), name = "Inferred\nancestry") +
  ylab("Proportion of clones") +
  xlab("Library")

plot_grid(a, b, rel_widths = c(1, 1.3))

ggplot(data = loc[pmax > 0.5], 
       aes(x = start, xend = end, y = chrom, yend = chrom, color = ancestry)) +
  geom_segment(size = 5) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c('#66c2a5','#fc8d62','#e78ac3', '#8da0cb','#a6d854'), name = "") +
  ylab("") +
  xlab("Position (bp)")

setorder(loc, chrom, start, end)
sum(loc[!is.na(chrom) & !is.na(ancestry) & pmax > 0.5]$len)

fwrite(loc[!is.na(chrom) & !is.na(ancestry), c("chrom", "start", "end", "accession", "CloneName", "LibAbbr", "ancestry", "nsnp", "pmax")],
       file = "~/Downloads/grch38_local_ancestry.txt",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE,
       col.names = TRUE)
