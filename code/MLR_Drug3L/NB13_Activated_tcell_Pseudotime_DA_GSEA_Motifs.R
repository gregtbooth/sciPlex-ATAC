basepath = "/home/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB13"))
setwd(paste0(out_dir, "results/NB13"))

suppressPackageStartupMessages({
  library(ArchR)
  #library(Seurat)
  #library(monocle3)
  library(devtools)
  load_all('/git/github/monocle3_dev') # latest developer branch of monocle3
  library(dplyr)
  library(tidymodels)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
  library(RColorBrewer)
  #library(ggpubr)
  library(viridis)
  library(snowfall)
  library(furrr)
  library(gt)
  plan(multicore)
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered_annotated"), 
                       showLogo = FALSE) 

# Load CDS with T-cell trajectory pseudotime information (then pull only allo-activated cells)
cds = readRDS(file = paste0(out_dir, "results/NB8/cds_NaiveActive_Tcells_PT"))
cds_f = cds[,colData(cds)$cellType_broad == "Activated_Tcell" & colData(cds)$Stimulator == "StimA"]

#######################
#######################

# perform DA analysis on activated T-cells 
# 1) Sites influenced by position in pseudotime (activated Tcells only). 

# Only need to run fit_models once: load the saved output in next section 

# subset to peaks open in > 0.5% cells
cds_a <- detect_genes(cds_f)
cutoff = ncol(cds_a)*.005
cds_a = cds_a[rowData(cds_a)$num_cells_expressed > cutoff, ]

# 1)
# fit regression models for each peak across with Pseudotime as predictor 
peak_fits_PT <- fit_models(
  cds = cds_a, 
  model_formula_str = "~Pseudotime + nFrags", 
  expression_family = "binomial", 
  cores=3)

# get coeficients from fit_models
peak_fit_coefs_PT = coefficient_table(peak_fits_PT)

peak_fit_coefs_PT_f = dplyr::select(peak_fit_coefs_PT, 
                                    peak, num_cells_expressed, status, term,
                                    estimate, std_err, test_val, p_value, 
                                    normalized_effect, model_component, q_value)
write.table(peak_fit_coefs_PT_f, 
            file =  "alloActTcell_Pseudotime_DApeak_coefs.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")

###############################
# Run DAG tests on Gene scores 
###############################

# load Gene Score Matrix into memory
gMat = getMatrixFromProject(ArchRProj = prj, 
                            useMatrix = "GeneScoreMatrix", 
                            binarize = FALSE)

# Create CDS from Gene Matrix (SummarizedExperiment)
g_rd = data.frame(rowData(gMat)) %>% 
  dplyr::mutate(gene_coords = paste0(seqnames, "_", start, "_", end)) %>% 
  dplyr::rename(gene_short_name = name) %>% 
  dplyr::select(gene_coords, gene_short_name, idx)
row.names(g_rd) <- g_rd$gene_short_name
row.names(gMat) <- g_rd$gene_short_name
cds_geneScores = monocle3::new_cell_data_set(assays(gMat)$GeneScoreMatrix, 
                                             cell_metadata = colData(gMat), 
                                             gene_metadata = g_rd)  

# filter for Allo-Activated Tcells
cds_gs_f = cds_geneScores[,colData(cds_geneScores)$cellType_broad == "Activated_Tcell" & colData(cds_geneScores)$Stimulator == "StimA"]
cds_gs = cds_gs_f
colData(cds_gs)$cell =
  row.names(colData(cds_gs))

# add PT data to geneScoreMatrix colData
cdat_pt = colData(cds_f) %>% data.frame() %>%
  select(cell, root_node, Pseudotime)

cdat_gs = colData(cds_gs) %>% data.frame() %>% 
  left_join(cdat_pt, by = "cell")
colData(cds_gs)$root_node = cdat_gs$root_node
colData(cds_gs)$Pseudotime = cdat_gs$Pseudotime

# subset to genes open in > 5% cells
cds_gs_a <- detect_genes(cds_gs)
cutoff = ncol(cds_gs_a)*.05
cds_gs_a = cds_gs_a[rowData(cds_gs_a)$num_cells_expressed > cutoff, ]

# Fit model to gene accessibility scores by Pseudotime position
gene_fits_PT <- fit_models(
  cds = cds_gs_a, 
  model_formula_str = "~Pseudotime + nFrags", 
  expression_family = "quasipoisson", 
  cores=3)

# get coeficients from fit_models
gene_fit_coefs_PT = coefficient_table(gene_fits_PT)

gene_fit_coefs_PT_f = dplyr::select(gene_fit_coefs_PT, 
                                    gene_short_name, num_cells_expressed, status, term,
                                    estimate, std_err, test_val, p_value, 
                                    normalized_effect, model_component, q_value)

write.table(gene_fit_coefs_PT_f, 
            file =  "alloActTcell_Pseudotime_DAgene_coefs.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")

#######################
#######################


# sites by Pseudotime Term
AlloUp_sites_PT = filter(peak_fit_coefs_PT_f,term == "Pseudotime", p_value <= 0.05, estimate > 0) %>% 
  mutate(chr = stringr::str_split_fixed(peak, "_", 3)[,1], 
         start = stringr::str_split_fixed(peak, "_", 3)[,2], 
         end = stringr::str_split_fixed(peak, "_", 3)[,3]) %>% 
  dplyr::select(chr, start, end, peak)

AlloDown_sites_PT = filter(peak_fit_coefs_PT_f, term == "Pseudotime", p_value <= 0.05, estimate < 0) %>% 
  mutate( chr = stringr::str_split_fixed(peak, "_", 3)[,1], 
         start = stringr::str_split_fixed(peak, "_", 3)[,2], 
         end = stringr::str_split_fixed(peak, "_", 3)[,3]) %>% 
  dplyr::select(chr, start, end, peak)

write.table(AlloUp_sites_PT, file = "AlloActT_Up_sites_PT.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(AlloDown_sites_PT, file = "AlloActT_Down_sites_PT.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)










