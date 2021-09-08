library(Matrix)
library(argparse)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrastr)
source("~/scripts/GB_src/io_functions.R")
source("~/scripts/GB_src/ConvertBCs_ToWells_V3.R")
basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190421_barnyard7/"
matpath = paste0(basepath, "pipeline_hill/analysis_barnyard/make_matrices/")
out_dir = paste0(basepath,"analysis/results/NB2/")
dir.create(out_dir)
setwd(out_dir)

################################

DOUBLET_PERCENTAGE_THRESHOLD = 0.9

df = read.table(file = paste0(basepath, "analysis/results/NB1/hashCellAssignments.txt"), head = T, stringsAsFactors = FALSE)

# change infinite hash enrichment values to Hash UMI count (since there is no second best hash)
df = df %>% mutate(df, pseudo_enrichment = ifelse(!is.infinite(top_to_second_best_ratio), 
                                                  top_to_second_best_ratio, hash_umis))


# filter low umi cells: 
umi_cutoff = 500
df = mutate(df, nFrags = hg19_frags + mm9_frags, ) %>% 
  dplyr::filter(nFrags >= umi_cutoff) %>% 
  rename(human = hg19_frags, mouse = mm9_frags)

#################################################################
# Barnyard of all cells conditions pooled 
ggplot(df) +
  geom_point_rast(aes(x =  human, y = mouse, color = cell_type), size = 1, alpha = 0.5) +
  xlab("Human reads/cell") +
  ylab("Mouse reads/cell") + 
  labs(color = "Hash label") +
  ylim(0,20000) + 
  xlim(0,20000) + 
  theme_bw() 
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = paste0(out_dir,"Barnyard_byCellType_pooled_rast.pdf"), width = 3.5, height = 2.25)

# Barnyard color doublets
ggplot(df) +
  geom_point_rast(aes(x =  human, y = mouse, color = doublet), size = 1, alpha = 0.5) +
  xlab("Human reads/cell") +
  ylab("Mouse reads/cell") + 
  labs(color = "Hash-based") +
  scale_color_manual(values=c("red", "black")) + 
  ylim(0,20000) + 
  xlim(0,20000) + 
  theme_bw() 
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = paste0(out_dir,"Barnyard_isdoublet_pooled_rast.pdf"), width = 3.5, height = 2.25)

# Valley Plot
ggplot(df) +
  geom_histogram(aes(x = fracHuman), bins = 50) +
  xlab("fraction human reads") +
  ylab("cells") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"BarnyardValley_pooled.pdf"), width = 2.5, height = 2.5)

# Hash Enrichment Scores
ggplot(df) +
  geom_histogram(aes(x = log10(top_to_second_best_ratio)), bins = 50) +
  xlab("log10(Enrichment Score)") +
  ylab("cells") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"EnrichmentScores_pooled.pdf"), width = 2.5, height = 2.5)

# Hash umis
ggplot(df) +
  geom_histogram(aes(x = log10(hash_umis)), bins = 50) +
  xlab("log10(hash UMIs)") +
  ylab("cells") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"HashUMIs_pooled.pdf"), width = 2.5, height = 2.5)



##############################  
#Facet Barnyards by hashing condition
ggplot(df) +
  facet_wrap(~Tn5_WellTreat_coding) + 
  geom_point_rast(aes(x =  human, y = mouse, color = cell_type), size = 4, alpha = 0.5) +
  xlab("Human reads/cell") +
  ylab("Mouse reads/cell") +
  labs(color = "Hash label") +
  ylim(0,10000) + 
  xlim(0,10000) + 
  theme_bw() 
ggsave(filename = paste0(out_dir,"Barnyard_byCellType_facet_rast.pdf"), width = 5.8, height = 3.2)

# faceted barnyard by doublet
ggplot(df) +
  geom_point_rast(aes(x =  human, y = mouse, color = doublet), size = 4, alpha = 0.5) +
  facet_wrap(~Tn5_WellTreat_coding) + 
  xlab("Human reads/cell") +
  ylab("Mouse reads/cell") +
  labs(color = "Hash-based") +
  scale_color_manual(values=c("red", "black")) + 
  ylim(0,10000) + 
  xlim(0,10000)+
  theme_bw() 
ggsave(filename = paste0(out_dir,"Barnyard_isDoublet_facet_rast.pdf"), width = 5.8, height = 3.2)

# faceted Valley Plot
ggplot(df) +
  facet_wrap(~Tn5_WellTreat_coding) + 
  geom_histogram(aes(x = fracHuman), bins = 50) +
  xlab("fraction human reads") +
  ylab("cells") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"BarnyardValley_facet.pdf"), width = 5, height = 3.2)

# faceted Hash Enrichment Scores
ggplot(df) +
  facet_wrap(~Tn5_WellTreat_coding) + 
  geom_histogram(aes(x = log10(top_to_second_best_ratio)), bins = 50) +
  xlab("log10(Enrichment Score)") +
  ylab("cells") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"EnrichmentScores_facet.pdf"), width = 5, height = 3.2)

##################################
##################################
# compare library properties accross conditions. 
# nFrags/cell
ggplot(df) +
  geom_boxplot(aes(x = Tn5_WellTreat_coding, y = log10(nFrags), 
                   fill = Tn5_WellTreat_coding)) +
  xlab("Condition") +
  ylab("log10(Frags)") +
  labs(fill = "Condition") +
  #scale_fill_viridis(discrete = TRUE, option = "D")+
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw() +
  theme(axis.text.x = element_blank()) 
ggsave(filename = paste0(out_dir,"FragsPerCell.pdf"), width = 3.5, height = 2.25)

# FRIP
#ggplot(df) +
#  geom_boxplot(aes(x = Tn5_WellTreat_coding, y = FRIP, 
#                   fill = Tn5_WellTreat_coding)) +
#  xlab("Condition") +
#  ylab("FRIP") +
#  labs(fill = "Condition") +
#  scale_fill_brewer(palette = "Dark2") + 
#  theme_bw() +
#  theme(axis.text.x = element_blank()) 
#ggsave(filename = paste0(out_dir,"FRIP.pdf"), width = 3.5, height = 2.25)

# Hash UMIs
ggplot(df) +
  geom_boxplot(aes(x = Tn5_WellTreat_coding, y = log10(hash_umis), 
                   fill = Tn5_WellTreat_coding)) +
  xlab("Condition") +
  ylab("log10(hash umis)") +
  labs(fill = "Condition") +
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw() +
  theme(axis.text.x = element_blank()) 
ggsave(filename = paste0(out_dir,"hash_UMIs.pdf"), width = 3.5, height = 2.25)

# Hash Enrichment
ggplot(df) +
  geom_boxplot(aes(x = Tn5_WellTreat_coding, y = log10(top_to_second_best_ratio), 
                   fill = Tn5_WellTreat_coding)) +
  xlab("Condition") +
  ylab("log10(hash Enrich.)") +
  labs(fill = "Condition") +
  scale_fill_brewer(palette = "Dark2") + 
  theme_bw() +
  theme(axis.text.x = element_blank()) 
ggsave(filename = paste0(out_dir,"hash_Enrichment.pdf"), width = 3.5, height = 2.25)

################################
# summarize hash vs chromatin doublet calling.

df1 =
  df %>% 
  mutate(chrom_doublet = 
           ifelse(fracHuman > 0.9, "singlet", 
                  ifelse(fracHuman < 0.1, "singlet", "doublet")), 
         chrom_cellType = ifelse(human > mouse, "A549", "3T3"))

# get number of correctly and incorrectly identified species 
df_summary = filter(df1, chrom_doublet == "singlet") %>% 
  mutate(label_correct = ifelse(cell_type == chrom_cellType, TRUE, FALSE)) %>% 
  group_by(label_correct) %>% 
  summarise(n())

print(df_summary)
  
# get number of species doublets (chromatin based) 
# correctly and incorrectly identified by hash
df_summary_1 = filter(df1, chrom_doublet == "doublet") %>% 
  mutate(label_correct = ifelse(doublet == chrom_doublet, TRUE, FALSE)) %>% 
  group_by(label_correct) %>% 
  summarise(n())

print(df_summary_1)

# get number of intra-species doublets (hash based only)

df_summary_2 = filter(df1, chrom_doublet == "singlet") %>% 
  group_by(doublet) %>% 
  summarise(cells = n(), 
            med_umi = median(human + mouse))
print(df_summary_2)













