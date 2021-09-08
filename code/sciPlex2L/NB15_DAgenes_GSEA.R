basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190521_scichem2_AllPlates/"
bin_directory = paste0(basepath, "analysis/bin/")
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB15"))
setwd(paste0(out_dir, "results/NB15"))

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(monocle3)
  library(cicero)
  library(tidymodels)
  library(devtools)
  library(monocle3)
  library(furrr)
  library(viridis)
  library(piano)
  library(UpSetR)
  library(snowfall)
  library(pheatmap)
  plan(multicore)
})

# source some code from sanjay (sciPlex paper)
suppressPackageStartupMessages({
  source(paste0(out_dir, "external_resources/bin/GSA_helper_functions.R"))
  source(paste0(out_dir, "external_resources/bin/loadGSCSafe.R"))
})

options(dplyr.summarise.inform = FALSE) 

set.seed(2017)

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered/"),
                       showLogo = FALSE)

##########################
########## GSEA ##########
##########################

# load RNA-based DEGs
simple_model_test_res = read.table(file = paste0(out_dir, "results/NB5/sciPlex_RNA_A549_small_screen_DEGs.txt"), head = TRUE)

# or 

# load ATAC-based DAGs
#simple_model_test_res = read.table(file = paste0(out_dir, "results/NB5/sciPlex_ATAC_A549_small_screen_DAgeneScores.txt"), head = TRUE) %>% 
#  mutate(id = gene_id)


simple_model_hit_table = simple_model_test_res %>%
  filter(grepl("dose_", term)) %>%
  mutate(hit = ifelse(q_value < 0.05, 1, 0)) %>%
  dplyr::select(id, term, hit) %>%
  spread(term, hit, fill = 0)


colnames(simple_model_hit_table) = stringr::str_replace(colnames(adjusted_model_hit_table), "dose_", "")
pdf(file = "RNA_DEG_hits_upset.pdf",
    onefile = FALSE,
    height = 3) # or other device
upset(
  as.data.frame(simple_model_hit_table),
  mainbar.y.label = "DEGs in intersection",
  sets.x.label = "Total DEGs"
)
dev.off()

########## GSEA ##########

## Load Gene Set Collection
hallmarks <-loadGSCSafe(file= paste0(out_dir, "external_resources/GeneSets/h.all.v6.0.symbols.gmt", sep = ""))

# This function is designed to collect GSA stats on each drug + cell_type combo
drug_effect_gsa <- function(drug_gene_test_res, gsc, min_genes_in_set=5, ncores=1){
  drug_gene_test_res = drug_gene_test_res %>% 
    group_by(gene_short_name) %>%
    filter(abs(test_val) == max(abs(test_val))) 
  directions=data.frame(row.names=as.character(drug_gene_test_res$gene_short_name), stats=drug_gene_test_res$test_val)
  geneLevelStats=data.frame(row.names=as.character(drug_gene_test_res$gene_short_name), stats=drug_gene_test_res$test_val)
  tryCatch({
    gsares <- runGSA(geneLevelStats=geneLevelStats, direction=geneLevelStats, geneSetStat="mean", gsc=gsc, nPerm=10000, ncpu=ncores, gsSizeLim=c(min_genes_in_set, Inf), verbose=T)
    return (gsares)
  }, error = function(e) { print(e); return (NA) } )
}
Sys.setenv(OMP_NUM_THREADS = 1)

simple_model_gsa_res = simple_model_test_res %>% filter(grepl("dose_", term)) %>% group_by(treatment, term) %>% 
  nest(-term) %>% 
  mutate(
    gsares = purrr::map(.f = drug_effect_gsa, .x = data, hallmarks, ncores=4)
  ) 

# Extract p values and GSA test statistics into a table
simple_model_gsa_summary = 
  simple_model_gsa_res %>%
  mutate(gsa_summary = purrr::map(.f = GSAsummaryTable, .x = gsares)) %>% 
  unnest(gsa_summary)


## Make a heatmap showing the effects of all the drugs on one cell line:
simple_model_gsa_summary = simple_model_gsa_summary %>%
  mutate(up_stat = `Stat (mix.dir.up)`, dn_stat = `Stat (mix.dir.dn)`) %>%
  dplyr::select(term,
                Name,
                `Genes (tot)`,
                `p (mix.dir.up)`,
                up_stat,
                `p (mix.dir.dn)`,
                dn_stat)
simple_model_gsa_summary = simple_model_gsa_summary %>% mutate(gsa_stat = ifelse((is.na(up_stat) == FALSE &
                                                                                    `p (mix.dir.up)` < `p (mix.dir.dn)`),
                                                                                 up_stat,
                                                                                 ifelse(is.na(dn_stat) == FALSE,-dn_stat, 0)
))

simple_model_gsa_stat = simple_model_gsa_summary %>% dplyr::select(Name, term, gsa_stat) %>% spread(term, gsa_stat)

gsa_effects = simple_model_gsa_stat

# Get the top gene sets for both
num_top_gene_sets = 8

top_simple_gene_sets =
  simple_model_gsa_summary %>%
  group_by(term) %>%
  top_n(num_top_gene_sets, abs(gsa_stat))

top_simple_gene_sets =
  top_simple_gene_sets %>%
  arrange(term, gsa_stat) %>%
  dplyr::select(Name, term, gsa_stat) %>%
  mutate(Name = stringr::str_split_fixed(Name, "%", 3)[, 1]) %>% as.data.frame
sig_pathway_names = unique(as.character(top_simple_gene_sets$Name))


# Make into a matrix
gsa_effect_matrix = gsa_effects %>% dplyr::select(-Name) %>% as.matrix
gsa_effect_matrix[is.na(gsa_effect_matrix)] = 0
row.names(gsa_effect_matrix) = stringr::str_split_fixed(gsa_effects$Name, "%", 3)[, 1]

# Make into a data frame
gsa_effect_df =
  gsa_effects %>%
  gather(key = term, value = gsa_stat, contains("dose"))

gsa_effect_df[is.na(gsa_effect_df)] = 0
gsa_effect_df$Name = stringr::str_split_fixed(gsa_effects$Name, "HALLMARK_", 2)[, 2]
gsa_effect_df$Name = stringr::str_replace(gsa_effect_df$Name, "_", " ")


gsa_effect_matrix[gsa_effect_matrix < -10] = -10
gsa_effect_matrix[gsa_effect_matrix > 10] = 10

matrix_to_plot = gsa_effect_matrix[sig_pathway_names, ]

colnames(matrix_to_plot)  = (colnames(matrix_to_plot) %>% stringr::str_split_fixed("_", 2))[, 2]

rn =   stringr::str_split_fixed(rownames(matrix_to_plot), pattern = "HALLMARK_", n =  2)[, 2]
rownames(matrix_to_plot) = stringr::str_replace_all(rn, "_", " ")
pheatmap(
  matrix_to_plot,
  cex = 0.5,
  filename = "sciPlex_RNA_DEG_GOheatmap.pdf",
  width = 3,
  #height = 2.25,
  height = 5,
  cellwidth = 8,
  cellheight = 8,
  fontsize = 10,
  fontsize_row = 25,
  fontsize_col = 25,
  legend = F,
  treeheight_row = 15,
  treeheight_col = 15,
  border_color = "black"
)
