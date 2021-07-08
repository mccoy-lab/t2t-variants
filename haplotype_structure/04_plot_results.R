library(data.table)
library(tidyverse)
library(ggbeeswarm)

rh <- fread("~/Downloads/rh_summary.txt", header = T)

igsr <- fread("~/Downloads/igsr_samples.txt")[, c(1, 4, 6)] %>%
  setnames(., c("sample_id", "pop", "superpop"))

rh <- merge(rh, igsr, by = "sample_id", all.x = TRUE)

rh[sample_id == "GRCh38", superpop := "GRCh38"]
rh[sample_id == "CHM13", superpop := "CHM13"]

rh$superpop <- factor(rh$superpop, 
                      levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "GRCh38", "CHM13"))

ggplot(data = rh[!is.na(superpop)], aes(x = superpop, y = n_discordant_snps, color = superpop)) + 
  geom_beeswarm() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_color_brewer(palette = "Set2", name = "") +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100, 1000, 10000), limits = c(100, 20000)) +
  xlab("") +
  ylab("Number of private recombinant SNPs")

ggplot(data = rh[!is.na(superpop)], aes(x = superpop, y = n_discordant_snps, color = superpop)) + 
  geom_beeswarm() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_color_brewer(palette = "Set2", name = "") +
  xlab("") +
  ylab("Number of private recombinant SNPs")

ggplot(data = rh[!is.na(superpop)], aes(x = superpop, y = n_discordant_haplotypes, color = superpop)) + 
  geom_beeswarm() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_color_brewer(palette = "Set2", name = "") +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100, 1000), limits = c(10, 2000)) +
  xlab("") +
  ylab("Number of private recombinant haplotypes")

ggplot(data = rh[!is.na(superpop)], aes(x = superpop, y = n_discordant_haplotypes, color = superpop)) + 
  geom_beeswarm() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_color_brewer(palette = "Set2", name = "") +
  xlab("") +
  ylab("Number of private recombinant haplotypes")




