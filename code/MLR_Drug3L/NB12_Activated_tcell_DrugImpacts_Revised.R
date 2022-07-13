basepath = "/home/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB12"))
setwd(paste0(out_dir, "results/NB12"))

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

# Revise pseudotime bins so that they only factor in the allo-activated T-cells
colData(cds_f)$Pseudotime_bin_new <- kmeans(colData(cds_f)$Pseudotime, centers=10)$cluster
# rename Pseudodose_bin to reflect ordering of average pseudodose
x = as.data.frame(colData(cds_f)) %>% group_by(Pseudotime_bin_new) %>% 
  mutate(mean_pt = mean(Pseudotime)) %>% ungroup() %>% 
  mutate(Pseudotime_bin_ordered = dense_rank((mean_pt))) %>% 
  dplyr::select(cell, Pseudotime, Pseudotime_bin_new, Pseudotime_bin_ordered)
colData(cds_f)$Pseudotime_bin_new = x$Pseudotime_bin_ordered

#print number of cells in each pseudodose bin 
print(table(colData(cds_f)$Pseudotime_bin_new))

#############################################################
# Run pseudobulk DA tests on peaks 
#############################################################
cdat_pt = colData(cds_f) %>% data.frame() %>%
  group_by(Pseudotime_bin_new) %>% 
  mutate(mean_pt = mean(Pseudotime)) %>% ungroup() %>% 
  mutate(PT_bin_rep = paste0("PTbin_", Pseudotime_bin_new, "_", Replicate)) %>% 
  select(cell, root_node, Pseudotime, Pseudotime_bin_new, PT_bin_rep, mean_pt)

colData(cds_f)$PT_bin_rep = cdat_pt$PT_bin_rep
colData(cds_f)$mean_pt = cdat_pt$mean_pt


# prep a grouping DF (PT position)
gdf = as.data.frame(colData(cds_f)) %>% 
  select(cell, group_id = PT_bin_rep, cell_group = mean_pt)

#### Code from Cole for regressing accessibility at each site by pseudobulked datasets: 
agg_expr_mat = monocle3::aggregate_gene_expression(cds_f,
                                                   cell_group_df = gdf,
                                                   norm_method="size_only",
                                                   scale_agg_values = FALSE,
                                                   pseudocount=0)

agg_coldata = gdf %>%
  dplyr::group_by(group_id, cell_group) %>%
  dplyr::summarize(num_cells_in_group = n()) %>%
  as.data.frame
agg_expr_mat = agg_expr_mat[,agg_coldata$group_id]
row.names(agg_coldata) = colnames(agg_expr_mat)


pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(cds_f) %>% as.data.frame)
pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)

#fit spline model for gene accessibility by mean pseudotime 
pseudo_fit = fit_models(pseudobulk_cds,
                        "~ splines::ns(cell_group, df = 3)",
                        #"~cell_group",
                        weights=colData(pseudobulk_cds)$num_cells_in_group,
                        expression_family = "quasipoisson",
                        cores=4)

pseudo_fit_coefs <- coefficient_table(pseudo_fit)

# Correct for multiple testing
simple_model_test_res =
  pseudo_fit_coefs %>%
  group_by(term) %>%
  mutate(q_value = p.adjust(p_value))

# save model output
smtr = select(simple_model_test_res, -model, -model_summary) 

saveRDS(smtr, 
        file =  "sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DA_pseudotime_peaks", 
)    

#############################################################
# Run DEG tests on Gene scores 
#############################################################

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
  group_by(Pseudotime_bin_new) %>% 
  mutate(mean_pt = mean(Pseudotime)) %>% ungroup() %>% 
  select(cell, root_node, Pseudotime, Pseudotime_bin_new, mean_pt)
  

cdat_gs = colData(cds_gs) %>% data.frame() %>% 
  left_join(cdat_pt, by = "cell") %>% 
  mutate(PT_bin_rep = paste0("PTbin_", Pseudotime_bin_new, "_", Replicate))
colData(cds_gs)$root_node = cdat_gs$root_node
colData(cds_gs)$Pseudotime = cdat_gs$Pseudotime
colData(cds_gs)$Pseudotime_bin = cdat_gs$Pseudotime_bin
colData(cds_gs)$PT_bin_rep = cdat_gs$PT_bin_rep
colData(cds_gs)$mean_pt = cdat_gs$mean_pt

# prep a grouping DF (PT position)
gdf = as.data.frame(colData(cds_gs)) %>% 
  select(cell, group_id = PT_bin_rep, cell_group = mean_pt)

#### Code from Cole for regressing accessibility at each site by pseudobulked datasets: 
agg_expr_mat = monocle3::aggregate_gene_expression(cds_gs,
                                                   cell_group_df = gdf,
                                                   norm_method="size_only",
                                                   scale_agg_values = FALSE,
                                                   pseudocount=0)

agg_coldata = gdf %>%
  dplyr::group_by(group_id, cell_group) %>%
  dplyr::summarize(num_cells_in_group = n()) %>%
  as.data.frame
agg_expr_mat = agg_expr_mat[,agg_coldata$group_id]
row.names(agg_coldata) = colnames(agg_expr_mat)


pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(cds_gs) %>% as.data.frame)
pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)


# fit spline model for gene accessibility by mean pseudotime 
pseudo_fit = fit_models(pseudobulk_cds,
                        "~ splines::ns(cell_group, df = 3)",
                        #"~cell_group",
                        weights=colData(pseudobulk_cds)$num_cells_in_group,
                        expression_family = "quasipoisson",
                        cores=4)

pseudo_fit_coefs <- coefficient_table(pseudo_fit)

# Correct for multiple testing
simple_model_test_res =
  pseudo_fit_coefs %>%
  group_by(term) %>%
  mutate(q_value = p.adjust(p_value))

# save model output
smtr = select(simple_model_test_res, -model, -model_summary) 

saveRDS(smtr, 
        file =  "sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DA_pseudotime_genes", 
        )    

smtr_t = readRDS(file = paste0(out_dir, "results/NB12/sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DA_pseudotime_genes"))

smtr_t %>%
  select(gene_short_name, estimate, std_err, test_val, p_value, q_value, status) %>% 
  filter(grepl("cell_group", term), p_value < 0.01, status != "FAIL") %>% 
  arrange(p_value) %>% 
  distinct(gene_short_name) %>% 
  head(n = 50) %>% data.frame()
#summarize(n())

#######################################
# make bigwig coverage tracks for allo activated T-cells for each pseudotime bin
#######################################

# add PT bins to ArchR project filtered for Allo-activated T-cells
prj_f = prj[row.names(prj@cellColData) %in% colData(cds_f)$cell,]

cdat = data.frame(colData(cds_f)) %>% 
  select(cell,root_node, Pseudotime, Pseudotime_bin_new) %>% 
  mutate(pt_bin = paste0("PTbin_", Pseudotime_bin_new))

pcdat = data.frame(prj_f@cellColData) %>% 
  tibble::rownames_to_column(var = "cell") %>% 
  left_join(cdat, by = "cell")

prj_f <- addCellColData(ArchRProj = prj_f, data = pcdat$pt_bin,
                      cells = pcdat$cell, name = "pt_bin", 
                      force = TRUE)

# separate tracks for all cell types within each condition (139 total)
getGroupBW(
  ArchRProj = prj_f,
  groupBy = "pt_bin",
  normMethod = "ReadsInTSS", # creates a group scale factor = 10k/sum(reads in TSS) 
  #normMethod = "nFrags", # creates a group scale factor = 10k/sum(nFrags)
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)


##################################################
# Function below designed to test for DA drug affected sites or genes 
# It takes a cds of cells and pulls all cells affected by a specific drug 
# It then groups cells by dose and creates pseudobulk ATAC profiles. 
# Pseudobulk signals within sites are then modeled as a function of drug dose 
DA_test_pseudobulk = function(cds, drug = "SAHA"){
  cds_gs_drug = cds[, colData(cds)$Drug == drug]
  # prep a grouping DF
  gdf = as.data.frame(colData(cds_gs_drug)) %>% 
    mutate(rel_dose = as.numeric(Relative_Dose), 
           group_id = paste0("dose_",Relative_Dose)) %>% 
    select(cell, group_id, cell_group = rel_dose)
  
  #### Code from Cole for regressing accessibility at each site by pseudobulked datasets: 
  agg_expr_mat = monocle3::aggregate_gene_expression(cds_gs_drug,
                                                     #cell_group_df = tibble::rownames_to_column(cds_f_SAHA@metadata[["cell_group_assignments"]]),
                                                     cell_group_df = gdf,
                                                     norm_method="size_only",
                                                     scale_agg_values = FALSE,
                                                     pseudocount=0)
  
  agg_coldata = gdf %>%
    dplyr::group_by(group_id, cell_group) %>%
    dplyr::summarize(num_cells_in_group = n()) %>%
    as.data.frame
  agg_expr_mat = agg_expr_mat[,agg_coldata$group_id]
  row.names(agg_coldata) = colnames(agg_expr_mat)
  
  pseudobulk_cds = new_cell_data_set(agg_expr_mat, cell_metadata = agg_coldata, rowData(cds_gs_drug) %>% as.data.frame)
  pseudobulk_cds = estimate_size_factors(pseudobulk_cds, round_exprs = FALSE)
  
  # fit models for each site/gene by drug dose
  pseudo_fit = fit_models(pseudobulk_cds,
                          "~cell_group",
                          weights=colData(pseudobulk_cds)$num_cells_in_group,
                          cores=4)
  
  pseudo_fit_coefs <- coefficient_table(pseudo_fit)
  
  # Correct for multiple testing
  simple_model_test_res =
    pseudo_fit_coefs %>%
    group_by(term) %>%
    mutate(q_value = p.adjust(p_value))
  
  # save model output
  smtr = select(simple_model_test_res, -model, -model_summary) 
  
  write.table(smtr, 
              file =  paste0("AlloActivatedTcell_DA_genes", drug,".txt"), 
              row.names=F, 
              quote=F, 
              sep="\t")
  return(smtr)
}
##################################################
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

## Run the pseudobulk DA test on aggregated data for each drug 
#drug_list = data.frame(colData(cds_gs)) %>% 
#  select(Drug) %>% 
#  distinct()

#for(d in drug_list$Drug){
#  cat("\nfinding DA genes for ", d, "\n")
#  DAres = DA_test_pseudobulk(cds_gs, drug = d)
#}


#cds_gs_drug = cds_gs[, colData(cds_gs)$Drug == "CycA-Rapa"]


# prep a grouping DF (relative drug dose)
#gdf = as.data.frame(colData(cds_gs)) %>% 
#  mutate(rel_dose = as.numeric(Relative_Dose), 
#         group_id = paste0("dose_",Relative_Dose)) %>% 
#  select(cell, group_id, cell_group = rel_dose)

## Run the pseudobulk DA test on aggregated data for each drug 
drug_list = data.frame(colData(cds_gs)) %>% 
  select(Drug) %>% 
  distinct()

for(d in drug_list$Drug){
  cat("\nfinding DA genes for ", d, "\n")
  DAres = DA_test_pseudobulk(cds_gs, drug = d)
}


