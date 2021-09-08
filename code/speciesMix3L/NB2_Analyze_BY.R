library(Matrix)
library(argparse)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrastr)
source("~/scripts/GB_src/io_functions.R")
source("~/scripts/GB_src/ParseBCtoWells.R")
basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191117_3Level_barnyard5/"
matpath = paste0(basepath, "pipeline_hill/analysis_barnyard/make_matrices/")
out_dir = paste0(basepath,"analysis/results/NB2/")
dir.create(out_dir)
setwd(out_dir)

# load hash assignments
df = read.table(file = paste0(basepath, "analysis/results/NB1/hashCellAssignments.txt"), head = T, stringsAsFactors = FALSE)

# change infinite hash enrichment values to Hash UMI count (since there is no second best hash)
df = mutate(df, pseudo_enrichment = ifelse(!is.infinite(top_to_second_best_ratio), 
                                             top_to_second_best_ratio, hash_umis))

#####################################################
# plots 

# Plot Faceted by annealing type
filter(df, hg19_frags + mm9_frags > 1000) %>% 
  ggplot() +
  geom_point_rast(aes(x =  hg19_frags, y = mm9_frags, color = celltype), size = 2, alpha = 0.5) +
  facet_wrap(~captureCode+OligoCode) +
  xlab("Human reads/cell") +
  ylab("Mouse reads/cell") +
  ylim(0,20000) + 
  xlim(0,20000) + 
  theme_bw() 
ggsave(filename = paste0(out_dir,"BarnyardScatter_facet_rast.pdf"), width = 5, height = 3.5)

filter(df, hg19_frags + mm9_frags > 1000) %>% 
  ggplot() +
  geom_point_rast(aes(x =  hg19_frags, y = mm9_frags, color = doublet), size = 2, alpha = 0.5) +
  facet_wrap(~captureCode+OligoCode) +
  xlab("Human reads/cell") +
  ylab("Mouse reads/cell") +
  ylim(0,20000) + 
  xlim(0,20000) +
  scale_color_manual(values=c("red", "black")) +
  theme_bw() 
ggsave(filename = paste0(out_dir,"BarnyardScatter_facet_HashDoublets_rast.pdf"), width = 5, height = 3.5)

# Hash Enrichment plots 
# Hash UMIs
filter(df, hg19_frags + mm9_frags > 1000) %>% 
  ggplot() +
  geom_histogram(aes(x = log10(hash_umis))) +
  facet_wrap(~captureCode+OligoCode) +
  xlab("log10(Hash UMIs)") +
  ylab("Frequency") +
  theme_bw() 
  #geom_vline(xintercept=log10(10), color = "red")
  ggsave(filename = paste0(out_dir,"Histogram_hashUMIs_facet.pdf"), width = 3.5, height = 3.5)

# Hash Enrichment
filter(df, hg19_frags + mm9_frags > 1000) %>% 
  ggplot() +
  geom_histogram(aes(x = log10(pseudo_enrichment))) +
  facet_wrap(~captureCode+OligoCode) +
  xlab("log10(Enrichment Ratio)") +
  ylab("Frequency") +
  theme_bw() 
  #geom_vline(xintercept=log10(5), color = "red")
  ggsave(filename = paste0(out_dir,"Histogram_hashEnrichment_facet.pdf"), width = 3.5, height = 3.5)

#####################

# Plot polydT + LNA cells only (for main figure) 
filter(df, hg19_frags + mm9_frags > 1000, captureCode == "polydT", OligoCode == "LNA", hash_umis > 10, top_to_second_best_ratio > 5) %>% 
  ggplot() +
  geom_point_rast(aes(x =  hg19_frags, y = mm9_frags, color = celltype), size = 0.6, alpha = 0.3) +
  xlab("Human reads/cell") +
  ylab("Mouse reads/cell") +
  ylim(0,20000) + 
  xlim(0,20000) + 
  theme_bw() 
ggsave(filename = paste0(out_dir,"BarnyardScatter_polydTLNA_Only_wColor_rast.pdf"), width = 3.5, height = 2.25)

filter(df, hg19_frags + mm9_frags > 1000, captureCode == "polydT", OligoCode == "LNA") %>% 
  ggplot() +
  geom_point_rast(aes(x =  hg19_frags, y = mm9_frags, color = doublet), size = 0.6, alpha = 0.3) +
  xlab("Human reads/cell") +
  ylab("Mouse reads/cell") +
  ylim(0,20000) + 
  xlim(0,20000) +
  scale_color_manual(values=c("red", "black")) +
  theme_bw() 
ggsave(filename = paste0(out_dir,"BarnyardScatter_polydTLNA_Only_HashDoublets_rast.pdf"), width = 3.5, height = 2.25)

# Hash Enrichment plots 
# Hash UMIs
filter(df, hg19_frags + mm9_frags > 1000, captureCode == "polydT", OligoCode == "LNA") %>% 
  ggplot() +
  geom_histogram(aes(x = log10(hash_umis))) +
  xlab("log10(Hash UMIs)") +
  ylab("Frequency") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"Histogram_hashUMIs_polydTLNA.pdf"), width = 2.25, height = 2.25)

# Hash Enrichment
filter(df, hg19_frags + mm9_frags > 1000, captureCode == "polydT", OligoCode == "LNA") %>% 
  ggplot() +
  geom_histogram(aes(x = log10(pseudo_enrichment))) +
  xlab("log10(Enrichment Ratio)") +
  ylab("Frequency") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"Histogram_hashEnrichment_polydTLNA.pdf"), width = 2.25, height = 2.25)













