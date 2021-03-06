basepath = "/home/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB9"))
setwd(paste0(out_dir, "results/NB9"))

suppressPackageStartupMessages({
  library(ArchR)
  #library(Seurat)
  library(monocle3)
  library(dplyr)
  library(tidymodels)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
  library(RColorBrewer)
  library(ggpubr)
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

#######################

#start from processed cds
cds = readRDS(file = paste0(out_dir, "results/NB8/cds_NaiveActive_Tcells_PT"))

# isolate only activated Tcells from allo stim 
cds_f = cds[,colData(cds)$cellType_broad == "Activated_Tcell" & colData(cds)$Stimulator == "StimA"]

#Plot UMAPs of activated Tcells by T-cell trajectory position: 
pdf("AlloActivated_Tcells_pseudotime_umap_byDrug.pdf", width = 3, height = 3)
colData(cds_f) %>%
  as.data.frame() %>%
  ggplot(aes( x = UMAP1, y =  UMAP2, color = Pseudotime))+
  geom_point_rast(size=2, stroke = 0.5) +
  facet_wrap(~Drug) +
  monocle3:::monocle_theme_opts() +
  #theme_void() +
  #theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
  theme(legend.position = "none") +
  scale_color_viridis(option = "inferno")
dev.off() 

#Plot UMAPs of activated Tcells by T-cell relative Dose: 
pdf("AlloActivated_Tcells_RelativeDose_umap_byDrug.pdf", width = 3, height = 3)
colData(cds_f) %>%
  as.data.frame() %>%
  ggplot(aes( x = UMAP1, y =  UMAP2, color = Relative_Dose))+
  geom_point_rast(size=2, stroke = 0.5) +
  facet_wrap(~Drug) +
  monocle3:::monocle_theme_opts() +
  #theme_void() +
  #theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
  theme(legend.position = "none") +
  scale_color_manual(values = c("0" = "#001889",
                                "1"="#74008D",
                                "10"="#C03180",
                                "100"="#ECBE00")) 
dev.off() 


# box plots comparing pseudotimes positions of allo stimulated Tcells by dose (facet drug) 
cdat = data.frame(colData(cds_f))

cdat %>%
  filter(is.finite(Pseudotime)) %>% 
  ggplot(aes(x = Relative_Dose, y = Pseudotime)) +
  geom_violin(aes(fill = Relative_Dose)) +
  geom_boxplot(width=0.1) +
  facet_wrap(~Drug,ncol = 3) +
  scale_fill_brewer(palette='Set1') +
  coord_flip() +
  theme_bw()+ 
  xlab("Relative Dose") +
  ylab("Pseudotime Position") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("0" = "#001889",
                                "1"="#74008D",
                                "10"="#C03180",
                                "100"="#ECBE00")) 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("BoxPlot_alloActivatedTcells_byDrugDose.pdf"),
       height = 5 , width = 5, unit = "in")
  
################
################
  
# perform DA analysis on activated T-cells 
# 1) Sites influenced by Drug Dose (Allo activated Tcells only).

# Only need to run fit_models once: load the saved output in next section 
  
# function for calling DE genes or sites
test_drug_dependence <- function(cds_subset, 
                                 cds, 
                                 reference_cells=NULL, 
                                 min_fraction_cells=0.005, 
                                 pseudodose=0.01, 
                                 residualModelFormulaTerms=NULL, 
                                 expFamily = "quasipoisson", 
                                 cores=1){
  cell_ids_to_test = as.character(cds_subset$cell)
  
  print(length(cell_ids_to_test))
  
  cds_subset = cds[,reference_cells]
  cds_subset = detect_genes(cds_subset)
  genes_in_reference_set = row.names(subset(fData(cds_subset), num_cells_expressed > ncol(cds_subset)*min_fraction_cells))
  
  cds_subset = cds[,cell_ids_to_test]
  cds_subset = detect_genes(cds_subset)
  genes_in_treated_set = row.names(subset(fData(cds_subset), num_cells_expressed > ncol(cds_subset)*min_fraction_cells))
  print(length(genes_in_treated_set))
  cds_subset = cds[,base::union(cell_ids_to_test, reference_cells)]
  cds_subset = cds_subset[base::union(genes_in_reference_set, genes_in_treated_set),]
  tmnt = unique(pData(cds_subset)$treatment)[1]
  ct = unique(pData(cds_subset)$cell_type)
  cell_treated = as.character(pData(cds_subset)$treatment) == as.character(tmnt)
  dose_of_drug = log((pData(cds_subset)$dose * cell_treated)  + pseudodose)
  drug_var_name = paste("dose_", tmnt, sep="")
  pData(cds_subset)[,drug_var_name] = dose_of_drug
  message(paste(ncol(cds_subset), ct, "cells treated with", tmnt))
  message(paste("Fitting drug/dose models to ", nrow(cds_subset), "genes"))
  modelFormulaStr = paste("~",  paste("dose_",tmnt, sep=""))
  
  Expression_Family = expFamily
  if (is.null(residualModelFormulaTerms) == FALSE & length(residualModelFormulaTerms) > 0){
    for (i in (1:length(residualModelFormulaTerms))){
      if (length(unique(pData(cds_subset)[,residualModelFormulaTerms[i]])) > 1){
        modelFormulaStr = paste(modelFormulaStr, "+", residualModelFormulaTerms[i])
      }
    }
  }
  if (Expression_Family %in% c("zipoisson", "zinegbinomial"))
    modelFormulaStr = paste(modelFormulaStr, "| log_hash_umis")
  message(paste("expression family : ", Expression_Family, sep = ""))
  message(paste("general linear model formula: ", modelFormulaStr, sep = ""))
  
  sci_plex_model_tbl = fit_models(cds_subset, expression_family = Expression_Family, model_formula_str = modelFormulaStr, cores=cores, reltol=1e-5, verbose=TRUE)
  sci_plex_model_coefs = coefficient_table(sci_plex_model_tbl)
  sci_plex_model_coefs
}

  
cds_atac = cds_f 

colData(cds_atac)$cell =
  row.names(colData(cds_atac))

colData(cds_atac)$log_nFrags =
  log(colData(cds_atac)$nFrags)

  
colData(cds_atac)$dose_character =
  factor(colData(cds_atac)$Relative_Dose,
         levels = c("0", "1", "10", "100"))
  
colData(cds_atac)$new_treatment_label =
  sapply(colData(cds_atac)$Drug, function(x) {
    if (grepl("Meth-CycA", x))
      return("Meth_CycA")
    if (grepl("CycA-Rapa", x))
      return("CycA_Rapa")
    if (grepl("BMS", x))
      return("BMS")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
    if (grepl("Rapa", x))
      return("Rapa")
    if (grepl("CycA", x))
      return("CycA")
    if (grepl("Meth", x))
      return("Meth")
  })

colData(cds_atac)$solvent =
  sapply(pData(cds_atac)$Drug, function(x) {
    if (grepl("Dex", x))
      return("Ethanol")
    else
      return("DMSO")
  })

colData(cds_atac)$vehicle =
  sapply(pData(cds_atac)$Relative_Dose, function(x) {
    if (x == 0)
      return(TRUE)
    else
      return(FALSE)
  })

colData(cds_atac)$treatment = colData(cds_atac)$new_treatment_label
colData(cds_atac)$cell_type = colData(cds_atac)$cellType_broad

# Use Relative doses rather than absolute concentrations for regressions. 
colData(cds_atac)$dose = as.numeric(colData(cds_atac)$Relative_Dose)


###########################################
# Run the above function for DAsite analysis (by drug)
sciPlex_cds = cds_atac

vehicle_cells = data.frame(pData(sciPlex_cds)) %>% filter(vehicle == TRUE)
vehicle_cells = as.character(vehicle_cells$cell)

reference_cells = sample_n(data.frame(pData(sciPlex_cds)), 1000)$cell

# run test_drug_dependence function on cds (group_by drug).
plan("multicore")
old_omp_num_threads = Sys.getenv("OMP_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 10)
simple_model_test_res =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  dplyr::mutate(
    test_res = purrr::map(
      data,
      .f = test_drug_dependence,
      sciPlex_cds,
      union(vehicle_cells, reference_cells),
      expFamily = "binomial",
      residualModelFormulaTerms = c("log_nFrags", "TSSEnrichment"),
      min_fraction_cells=0.001,
      cores = 1
    )
  )

simple_model_test_res = 
  simple_model_test_res %>% 
  dplyr::select(-data) %>% 
  unnest()

Sys.setenv(OMP_NUM_THREADS = old_omp_num_threads)


# Correct for multiple testing
simple_model_test_res =
  simple_model_test_res %>%
  group_by(term) %>%
  mutate(q_value = p.adjust(p_value))


# save model output
smtr = select(simple_model_test_res, -model, -model_summary)

write.table(smtr, 
            file =  "sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DApeaks.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")


# create upset plot of DA sites
simple_model_hit_table = smtr %>%
  filter(grepl("dose_", term)) %>%
  mutate(hit = ifelse(p_value < 0.01, 1, 0)) %>%
  dplyr::select(gene_id, term, hit) %>%
  spread(term, hit, fill = 0)


pdf(file =  "sciPlexATAC_alloActTcells_DAsites_upset.pdf",
    onefile = FALSE,
    height = 5, width = 10) # or other device
upset(
  as.data.frame(simple_model_hit_table),
  mainbar.y.label = "DA sites in intersection",
  sets.x.label = "Total DA sites", 
  nsets = 10
)
dev.off()


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

colData(cds_gs)$log_nFrags =
  log(colData(cds_gs)$nFrags)


colData(cds_gs)$dose_character =
  factor(colData(cds_gs)$Relative_Dose,
         levels = c("0", "1", "10", "100"))

colData(cds_gs)$new_treatment_label =
  sapply(colData(cds_gs)$Drug, function(x) {
    if (grepl("Meth-CycA", x))
      return("Meth_CycA")
    if (grepl("CycA-Rapa", x))
      return("CycA_Rapa")
    if (grepl("BMS", x))
      return("BMS")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
    if (grepl("Rapa", x))
      return("Rapa")
    if (grepl("CycA", x))
      return("CycA")
    if (grepl("Meth", x))
      return("Meth")
  })

colData(cds_gs)$solvent =
  sapply(pData(cds_gs)$Drug, function(x) {
    if (grepl("Dex", x))
      return("Ethanol")
    else
      return("DMSO")
  })

colData(cds_gs)$vehicle =
  sapply(pData(cds_gs)$Relative_Dose, function(x) {
    if (x == 0)
      return(TRUE)
    else
      return(FALSE)
  })

colData(cds_gs)$treatment = colData(cds_gs)$new_treatment_label
colData(cds_gs)$cell_type = colData(cds_gs)$cellType_broad

# Use Relative doses rather than absolute concentrations for regressions. 
colData(cds_gs)$dose = as.numeric(colData(cds_gs)$Relative_Dose)

sciPlex_cds = cds_gs

vehicle_cells = data.frame(colData(sciPlex_cds)) %>% 
  dplyr::filter(vehicle == TRUE)
vehicle_cells = as.character(vehicle_cells$cell)

# draw reference cells from only vehicle treated population. 
reference_cells = dplyr::sample_n(data.frame(pData(sciPlex_cds[,vehicle_cells])), 500)$cell

# run test_drug_dependence function on cds (group_by drug).
future::plan("multiprocess")
old_omp_num_threads = Sys.getenv("OMP_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 10)
simple_model_test_res =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  dplyr::group_by(cell_type, treatment) %>%
  tidyr::nest_legacy() %>%
  dplyr::mutate(
    test_res = purrr::map(
      data,
      .f = test_drug_dependence,
      sciPlex_cds,
      union(vehicle_cells, reference_cells),
      residualModelFormulaTerms=NULL, 
      expFamily = "quasipoisson", 
      min_fraction_cells=0.01,
      cores = 1
    )
  )

simple_model_test_res = 
  simple_model_test_res %>% 
  dplyr::select(-data) %>% 
  tidyr::unnest_legacy()

Sys.setenv(OMP_NUM_THREADS = old_omp_num_threads)


# Correct for multiple testing
simple_model_test_res =
  simple_model_test_res %>%
  group_by(term) %>%
  mutate(q_value = p.adjust(p_value))


# save model output
smtr = select(simple_model_test_res, -model, -model_summary)

write.table(smtr, 
            file =  "sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DAgenes_NEW.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")


smtr = read.table(file = paste0(out_dir,"results/NB10/sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DAgenes_NEW.txt"), header = TRUE)

# create upset plot of DA sites
simple_model_hit_table = smtr %>%
  filter(grepl("dose_", term)) %>%
  mutate(hit = ifelse(q_value < 0.05, 1, 0)) %>%
  dplyr::select(gene_id, term, hit) %>%
  spread(term, hit, fill = 0)


pdf(file =  "sciPlexATAC_alloActTcells_DAgenes_qval_upset_NEW.pdf",
    onefile = FALSE,
    height = 5, width = 10) # or other device
upset(
  as.data.frame(simple_model_hit_table),
  mainbar.y.label = "DA genes in intersection",
  sets.x.label = "Total DA genes", 
  nsets = 10
)
dev.off()


  















  ##########################
  #### Testing ###########
  ##########################
  
  # function for calling DE genes or sites
  test_drug_dependence <- function(cds_subset, 
                                   cds, 
                                   reference_cells=NULL, 
                                   min_fraction_cells=0.005, 
                                   pseudodose=0.01, 
                                   residualModelFormulaTerms=NULL, 
                                   expFamily = "quasipoisson", 
                                   cores=1){
    cell_ids_to_test = as.character(cds_subset$cell)
    
    print(length(cell_ids_to_test))
    
    cds_subset = cds[,reference_cells]
    cds_subset = detect_genes(cds_subset)
    genes_in_reference_set = row.names(subset(fData(cds_subset), num_cells_expressed > ncol(cds_subset)*min_fraction_cells))
    
    cds_subset = cds[,cell_ids_to_test]
    cds_subset = detect_genes(cds_subset)
    genes_in_treated_set = row.names(subset(fData(cds_subset), num_cells_expressed > ncol(cds_subset)*min_fraction_cells))
    print(length(genes_in_treated_set))
    cds_subset = cds[,base::union(cell_ids_to_test, reference_cells)]
    cds_subset = cds_subset[base::union(genes_in_reference_set, genes_in_treated_set),]
    tmnt = unique(pData(cds_subset)$treatment)[1]
    ct = unique(pData(cds_subset)$cell_type)
    cell_treated = as.character(pData(cds_subset)$treatment) == as.character(tmnt)
    dose_of_drug = log((pData(cds_subset)$dose * cell_treated)  + pseudodose)
    drug_var_name = paste("dose_", tmnt, sep="")
    pData(cds_subset)[,drug_var_name] = dose_of_drug
    message(paste(ncol(cds_subset), ct, "cells treated with", tmnt))
    message(paste("Fitting drug/dose models to ", nrow(cds_subset), "genes"))
    modelFormulaStr = paste("~",  paste("dose_",tmnt, sep=""))
    
    Expression_Family = expFamily
    if (is.null(residualModelFormulaTerms) == FALSE & length(residualModelFormulaTerms) > 0){
      for (i in (1:length(residualModelFormulaTerms))){
        if (length(unique(pData(cds_subset)[,residualModelFormulaTerms[i]])) > 1){
          modelFormulaStr = paste(modelFormulaStr, "+", residualModelFormulaTerms[i])
        }
      }
    }
    if (Expression_Family %in% c("zipoisson", "zinegbinomial"))
      modelFormulaStr = paste(modelFormulaStr, "| log_hash_umis")
    message(paste("expression family : ", Expression_Family, sep = ""))
    message(paste("general linear model formula: ", modelFormulaStr, sep = ""))
    
    sci_plex_model_tbl = fit_models(cds_subset, expression_family = Expression_Family, model_formula_str = modelFormulaStr, cores=cores, reltol=1e-5, verbose=TRUE)
    sci_plex_model_coefs = coefficient_table(sci_plex_model_tbl)
    sci_plex_model_coefs
  }
  
  
  
  

cds = cds_gs
colData(cds)$dose = as.numeric(colData(cds)$Relative_Dose)

cds_subset = cds_gs
colData(cds_subset)$dose = as.numeric(colData(cds_subset)$Relative_Dose)


# run through the test_drug_dependence function
reference_cells = dplyr::sample_n(data.frame(pData(cds)), 1000)$cell
min_fraction_cells=0.005
pseudodose=0.01
cell_ids_to_test = as.character(cds_subset$cell)
expFamily = "quasipoisson"
cores=1

print(length(cell_ids_to_test))

cds_subset_1 = cds[,reference_cells]
cds_subset_1 = detect_genes(cds_subset_1)
genes_in_reference_set = row.names(subset(fData(cds_subset_1), num_cells_expressed > ncol(cds_subset_1)*min_fraction_cells))

cds_subset_2 = cds[,cell_ids_to_test]
cds_subset_2 = detect_genes(cds_subset_2)

genes_in_treated_set = row.names(subset(fData(cds_subset_2), num_cells_expressed > ncol(cds_subset_2)*min_fraction_cells))
print(length(genes_in_treated_set))


cds_subset_3 = cds[,base::union(cell_ids_to_test, reference_cells)]
cds_subset_3 = cds_subset_3[base::union(genes_in_reference_set, genes_in_treated_set),]

tmnt = unique(pData(cds_subset_3)$treatment)[1]

ct = unique(pData(cds_subset_3)$cell_type)

## this line defines only cells with treatment of interest as treated. (All others will be defined as untreated)
cell_treated = as.character(pData(cds_subset_3)$treatment) == as.character(tmnt)

dose_of_drug = log((pData(cds_subset_3)$dose * cell_treated)  + pseudodose)

drug_var_name = paste("dose_", tmnt, sep="")

pData(cds_subset)[,drug_var_name] = dose_of_drug

message(paste(ncol(cds_subset), ct, "cells treated with", tmnt))

message(paste("Fitting drug/dose models to ", nrow(cds_subset), "genes"))

modelFormulaStr = paste("~",  paste("dose_",tmnt, sep=""))

sci_plex_model_tbl = fit_models(cds_subset_3, expression_family = Expression_Family, 
                                model_formula_str = modelFormulaStr, cores=cores, 
                                reltol=1e-5, verbose=TRUE)

sci_plex_model_coefs = coefficient_table(sci_plex_model_tbl)
sci_plex_model_coefs







