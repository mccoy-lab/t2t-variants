library(data.table)
library(tidyverse)
library(ggbeeswarm)

dt <- fread("~/Downloads/singleton_total.txt") %>%
  setnames(., c("singleton_count", "sample_id"))

igsr <- fread("~/Downloads/igsr_samples.txt")[, c(1, 4, 6)] %>%
  setnames(., c("sample_id", "pop", "superpop"))

dt <- merge(dt, igsr, by = "sample_id", all.x = TRUE)

dt[, haploid_singleton_count := singleton_count]
dt[sample_id != "CHM13", haploid_singleton_count := singleton_count / 2]

dt[sample_id == "CHM13", superpop := "CHM13"]

set.seed(1)
dt <- dt[, .SD[sample(.N, min(100,.N))], by = superpop]

dt$superpop <- factor(dt$superpop, levels = c("AFR", "AMR", "EAS", "EUR", "SAS", "CHM13"))

ggplot(data = dt[!is.na(superpop)], aes(x = superpop, y = haploid_singleton_count, color = superpop)) + 
  geom_beeswarm() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_color_brewer(palette = "Set2", name = "") +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100, 1000, 10000), limits = c(100, 20000)) +
  xlab("") +
  ylab("Ploidy-Adjusted Count of Singleton Alleles")
