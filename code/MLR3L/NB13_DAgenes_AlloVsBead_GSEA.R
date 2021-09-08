basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB13"))
setwd(paste0(out_dir, "results/NB13"))

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(monocle3)
  library(cicero)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
  library(Gviz)
  library(biomaRt)
  library(RColorBrewer)
  library(viridis)
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

prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 


#######################
#######################

# load PeakMatrix into memory
gMat = getMatrixFromArrow(ArrowFile = paste0(out_dir, "mlr_filtered_annotated/ArrowFiles/MLR.arrow"),
                          ArchRProj = prj, 
                          cellNames= prj$cellNames, # restricts to filtered cells
                          useMatrix = "GeneScoreMatrix", 
                          binarize = FALSE)

# Create CDS from gene score Matrix (SummarizedExperiment)
g_rd = data.frame(rowData(gMat)) %>% 
  dplyr::mutate(gene_coords = paste0(seqnames, "_", start, "_", end)) %>% 
  dplyr::rename(gene_short_name = name) %>% 
  dplyr::select(gene_coords, gene_short_name, idx)
row.names(g_rd) <- g_rd$gene_short_name
row.names(gMat) <- g_rd$gene_short_name
cds_geneScores = monocle3::new_cell_data_set(assays(gMat)$GeneScoreMatrix, 
                                             cell_metadata = colData(gMat), 
                                             gene_metadata = g_rd)

#broadly define allo and bead stimulations
cdat = colData(cds_geneScores) %>% 
  data.frame() %>%
  tibble::rownames_to_column(var = "cell") %>% 
  mutate(stim_type = 
           ifelse(Responder == "stimAlone", "stimAlone", 
                  ifelse(Stimulator == "Bead", "bead", 
                         ifelse(Stimulator == "noStim", "noStim",
                                ifelse(stringr::str_sub(Stimulator, -1) == stringr::str_sub(Responder, -1), "auto", "allo")))))

colData(cds_geneScores)$stim_type = cdat$stim_type

# Attach Pseudotime values for cells from peak cds
cds_peak = readRDS(file = paste0(out_dir,"results/NB7/cds_NaiveActive_Tcells_PT"))
cdat_peak = colData(cds_peak) %>% 
  data.frame() %>% 
  dplyr::select(cell, Pseudotime)

cdat = left_join(cdat, cdat_peak, by = "cell")
colData(cds_geneScores)$Pseudotime = cdat$Pseudotime

# count Activated Tcells from each broad stim_type
# including no stim etc. 

T_act_counts = filter(cdat, cellType_broad == "Activated_Tcell") %>% 
  group_by(stim_type) %>% 
  summarise(Tcells_act = n()) 

TcellCounts = group_by(cdat, stim_type) %>% 
  summarise(Tcells_total = n()) %>% 
  left_join(T_act_counts, by = "stim_type")

############
# Isolate only the bead or allo-Activated Tcells 
idx_cells = 
  (colData(cds_geneScores)$cellType_broad == "Activated_Tcell" & 
     (colData(cds_geneScores)$stim_type == "allo" | colData(cds_geneScores)$stim_type == "bead"))

cds_a = cds_geneScores[,idx_cells]


#######################
#######################

# perform DA analysis on activated T-cells 2 ways
# 1) Genes (accessibility) influenced by Bead vs Allo stim (activated Tcells only).
# 2) Genes (accessibility) influenced by position in pseudotime (activated Tcells only). 

# Only need to run fit_models once: load the saved output in next section 

# subset to peaks open in > 0.5% cells
cds_a <- detect_genes(cds_a)
cutoff = ncol(cds_a)*.005
cds_f = cds_a[rowData(cds_a)$num_cells_expressed > cutoff, ]

# 1)
# fit regression models for each peak across with stimulation type as predictor 
GA_fits_stimType <- fit_models(
  cds = cds_f, 
  model_formula_str = "~stim_type + nFrags", 
  expression_family = "quasipoisson")

# get coeficients from fit_models
GA_fit_coefs_ST = coefficient_table(GA_fits_stimType)

GA_fit_coefs_ST_f = dplyr::select(
  GA_fit_coefs_ST, gene_coords, gene_short_name, 
  num_cells_expressed, status, term,estimate, 
  std_err, test_val, p_value,  normalized_effect, 
  model_component, q_value)

write.table(GA_fit_coefs_ST_f, 
            file =  "Tcell_StimType_DAgene_coefs.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")

# 2)
# fit regression models for each gene accessibility across with Pseudotime as predictor 
GA_fits_PT <- fit_models(
  cds = cds_f, 
  model_formula_str = "~Pseudotime + nFrags", 
  expression_family = "quasipoisson")

# get coeficients from fit_models
GA_fit_coefs_PT = coefficient_table(GA_fits_PT)

GA_fit_coefs_PT_f = dplyr::select(
  GA_fit_coefs_PT, gene_coords, gene_short_name,
  num_cells_expressed, status, term, estimate, 
  std_err, test_val, p_value, normalized_effect, 
  model_component, q_value)

write.table(GA_fit_coefs_PT_f, 
            file =  "Tcell_Pseudotime_DAgene_coefs.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")

##########################
########## GSEA ##########
##########################

# load DA coefficients comparing Allo vs. Bead Stim. T-cells. 
# DA based on Stimulation type (ST)
# NOTE:  Here positive coefs mean open in Bead relative to Allo
coefs_ST = read.table("Tcell_StimType_DAgene_coefs.txt", head = TRUE) %>% 
  filter(term == "stim_typebead") 

# repeat analysis with DA genes based on pseudotime
# DA based on Pseudotime (PT) %>% 
coefs_PT = read.table("Tcell_Pseudotime_DAgene_coefs.txt", head = TRUE) %>% 
  filter(term == "Pseudotime")

## Load Gene Set Collection
hallmarks <-loadGSCSafe(file= paste0(out_dir,"external_resources/GeneSets/h.all.v6.0.symbols.gmt"))

# This function is designed to collect GSA stats on DAgenes 
get_gsa_res <- function(DAgene_test_res, gsc, min_genes_in_set=5, ncores=1){
  DAgene_test_res = DAgene_test_res %>% 
    group_by(gene_short_name) %>%
    filter(abs(test_val) == max(abs(test_val))) 
  directions=data.frame(row.names=as.character(DAgene_test_res$gene_short_name), stats=DAgene_test_res$test_val)
  geneLevelStats=data.frame(row.names=as.character(DAgene_test_res$gene_short_name), stats=DAgene_test_res$test_val)
  tryCatch({
    gsares <- runGSA(geneLevelStats=geneLevelStats, direction=geneLevelStats, geneSetStat="mean", gsc=gsc, nPerm=10000, ncpu=ncores, gsSizeLim=c(min_genes_in_set, Inf), verbose=T)
    return (gsares)
  }, error = function(e) { print(e); return (NA) } )
}
Sys.setenv(OMP_NUM_THREADS = 1)


# run GSEA using DA genes (comparing Allo vs Bead directly)
DAgene_gsa_res = coefs_PT %>% 
  nest(-term) %>% 
  mutate(
    gsares = purrr::map(.f = get_gsa_res, .x = data, hallmarks, ncores=4)
  ) 

# Extract p values and GSA test statistics into a table
DAgene_gsa_summary = 
  DAgene_gsa_res %>%
  mutate(gsa_summary = purrr::map(.f = GSAsummaryTable, .x = gsares)) %>% 
  unnest(gsa_summary)


## Make a heatmap showing the effects of all the drugs on one cell line:
simple_model_gsa_summary = DAgene_gsa_summary %>%
  mutate(up_stat = `Stat (mix.dir.up)`, dn_stat = `Stat (mix.dir.dn)`) %>%
  dplyr::select(term,
                Name,
                `Genes (tot)`,
                `p (mix.dir.up)`,
                up_stat,
                `p (mix.dir.dn)`,
                dn_stat)
simple_model_gsa_summary = simple_model_gsa_summary %>% 
  mutate(gsa_stat = ifelse((is.na(up_stat) == FALSE & `p (mix.dir.up)` < `p (mix.dir.dn)`),
                           up_stat, ifelse(is.na(dn_stat) == FALSE, -dn_stat, 0)
))

simple_model_gsa_stat = simple_model_gsa_summary %>% 
  dplyr::select(Name, term, gsa_stat) %>% 
  spread(term, gsa_stat)

gsa_effects = simple_model_gsa_stat

# Get the top gene sets for both
num_top_gene_sets = 10

top_simple_gene_sets =
  simple_model_gsa_summary %>%
  group_by(term) %>%
  top_n(num_top_gene_sets, abs(gsa_stat))

top_simple_gene_sets =
  top_simple_gene_sets %>%
  arrange(term, gsa_stat) %>%
  dplyr::select(Name, term, gsa_stat) %>%
  mutate(Name = stringr::str_split_fixed(Name, "%", 3)[, 1], 
         direction = ifelse(gsa_stat > 0, "up", "down"), 
         Name_short = stringr::str_split_fixed(Name, "HALLMARK_", 2)[, 2]) %>% 
  as.data.frame
sig_pathway_names = unique(as.character(top_simple_gene_sets$Name))

# prepare barplot of Gene set enrichments
top_simple_gene_sets$Name_fact <- 
  factor(top_simple_gene_sets$Name_short, levels = top_simple_gene_sets$Name_short[order(top_simple_gene_sets$gsa_stat)])

ggplot(top_simple_gene_sets, aes(x = Name_fact, y = gsa_stat, fill = direction)) +
  geom_col() +
  scale_fill_brewer(palette = "Set1") +
  xlab("pathway") + 
  ylab("directionality") +
  coord_flip() +
  theme_bw()

ggsave(file = "StimType_GSEA_barplot.pdf", width = 5, height = 3)



##### GSEA heatmap #############
# Make into a matrix
gsa_effect_matrix = gsa_effects %>% dplyr::select(-Name) %>% as.matrix
gsa_effect_matrix[is.na(gsa_effect_matrix)] = 0
row.names(gsa_effect_matrix) = stringr::str_split_fixed(gsa_effects$Name, "%", 3)[, 1]

# Make into a data frame
gsa_effect_df =
  gsa_effects %>%
  gather(key = term, value = gsa_stat, contains("stim"))

gsa_effect_df[is.na(gsa_effect_df)] = 0
gsa_effect_df$Name = stringr::str_split_fixed(gsa_effects$Name, "HALLMARK_", 2)[, 2]
gsa_effect_df$Name = stringr::str_replace(gsa_effect_df$Name, "_", " ")


gsa_effect_matrix[gsa_effect_matrix < -10] = -10
gsa_effect_matrix[gsa_effect_matrix > 10] = 10

matrix_to_plot = as.matrix(gsa_effect_matrix[sig_pathway_names, ])

colnames(matrix_to_plot)  = (colnames(matrix_to_plot) %>% stringr::str_split_fixed("_", 2))[, 2]

rn =   stringr::str_split_fixed(rownames(matrix_to_plot), pattern = "HALLMARK_", n =  2)[, 2]
rownames(matrix_to_plot) = stringr::str_replace_all(rn, "_", " ")
pheatmap(
  matrix_to_plot,
  cex = 0.5,
  filename = "Pseudotime_GOheatmap.pdf",
  width = 3,
  height = 2.25,
  cellwidth = 8,
  cellheight = 8,
  fontsize = 10,
  fontsize_row = 25,
  fontsize_col = 25,
  legend = F,
  cluster_cols = FALSE,
  treeheight_row = 15,
  treeheight_col = 15,
  border_color = "black"
)
###########################
