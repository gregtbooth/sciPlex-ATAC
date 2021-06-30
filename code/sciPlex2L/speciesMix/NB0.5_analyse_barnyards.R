
suppressPackageStartupMessages({
  library(Matrix)
  library(argparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrastr)
})

source("~/scripts/GB_src/io_functions.R")
source("~/scripts/GB_src/ConvertBCs_ToWells_V3.R")

basepath = "github/"
matpath = paste0(basepath, "pipeline_hill/analysis_barnyard/make_matrices/")
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/barnyard"))
setwd(paste0(out_dir, "results/barnyard"))

# load hash assignments
HashAsmnts = read.table(file = paste0(out_dir, "hash/hashCellAssignments.txt"), head = T, stringsAsFactors = FALSE)

# load genome window matrix
sample_window_matrix = load_mtx_file(paste0(matpath, "barnyard.window_matrix.mtx.gz"))

# count human and mouse reads per cell barcode
hg19_counts = Matrix::colSums(sample_window_matrix[grepl('hg19', rownames(sample_window_matrix)),])
mm9_counts = Matrix::colSums(sample_window_matrix[grepl('mm9', rownames(sample_window_matrix)),])

if (sum(hg19_counts) == 0) {
  stop('No hg19 counts found. Chromosomes must start with hg19 in name (e.g. hg19chr19) for human and mm9 for mouse.')
}
if (sum(mm9_counts) == 0) {
  stop('No mm9 counts found. Chromosomes must start with mm9 in name (e.g. mm9chr19) for mouse and hg19 for human.')
}

barnyard_df = data.frame(cell=names(mm9_counts), hg19=hg19_counts, mm9=mm9_counts)

# merge cells with hashing info 
by_df = inner_join(barnyard_df, HashAsmnts, by = "cell") %>% 
  filter(treatment %in% c("human", "mouse", "PREmixed")) %>% 
  mutate(fracHuman = (hg19/(hg19+mm9)), 
         nFrags = hg19+mm9, 
         mixing = ifelse(treatment == "PREmixed", treatment, "POSTmixed"), 
         hashEnrichment = ifelse(is.finite(top_to_second_best_ratio), 
                                 top_to_second_best_ratio, hash_umis)) %>%
  select(cell, hash_umis, pval, qval, top_to_second_best_ratio, 
         hashEnrichment, cell_type, treatment, mixing, doublet,
         hg19_frags = hg19, mm9_frags = mm9,fracHuman, 
         nFrags)


write.table(by_df, file = "SpeciesMix_hashCellAssignments.txt", 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

##################################
# barnyard plots  

#by_df = read.table(file = "SpeciesMix_hashCellAssignments.txt", head = TRUE)

filter(by_df, nFrags >= 500) %>% 
  ggplot() +
    facet_wrap(~mixing)+
    geom_point_rast(aes(x = hg19_frags/1000, y = mm9_frags/1000, color = treatment), 
                    size = 2, alpha = 0.75) +
    xlab("Human reads/cell (1k)") +
    ylab("Mouse reads/cell (1k)") + 
    ylim(0,60) + 
    xlim(0,60) + 
    theme_bw()
ggsave(filename = "BarnyardScatter_Facet.pdf",
       width = 5, height = 2.25)

# Valley Plots
filter(by_df, nFrags >= 500) %>% 
  ggplot() +
  facet_wrap(~mixing)+
  geom_histogram(aes(x = fracHuman), bins = 100)+
  xlab("fraction human reads") +
  ylab("cells") +
  theme_bw()
ggsave(filename = "valleyPlot_Facet.pdf", 
       width = 4.5, height = 2.25)

# Hash Enrichment Scores
filter(by_df, nFrags >= 500) %>% 
  ggplot() +
  geom_histogram(aes(x = log10(hashEnrichment)))+
  facet_wrap(~mixing)+
  xlab("log10(Enrichment Score)") +
  ylab("cells") +
  theme_bw()
ggsave(filename = "hashEnrichment_Facet.pdf", 
       width = 4.5, height = 2.25)

# get collision rates for pre and post hash samples
c_rate = filter(by_df, nFrags >= 500) %>%
  group_by(mixing) %>% 
  summarise(tot_cells = n(),
            collisions = sum(fracHuman < 0.9 & fracHuman > 0.1)) %>% 
  mutate(fracCollision = collisions/tot_cells)

c_rate


