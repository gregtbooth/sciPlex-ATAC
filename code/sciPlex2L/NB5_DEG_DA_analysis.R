basepath = "github/"
bin_directory = paste0(basepath, "analysis/bin/")
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB5"))
setwd(paste0(out_dir, "results/NB5"))

suppressPackageStartupMessages({
  library(tidymodels)
  library(tidyr)
  library(ArchR)
  library(monocle3)
  library(cicero)
  library(ggplot2)
  library(ggridges)
  library(devtools)
  library(furrr)
  library(viridis)
  library(piano)
  library(UpSetR)
  library(snowfall)
  library(dplyr)
  library(ggrastr)
  library(ggpubr)
  plan(multicore)
})

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)
set.seed(1989) 


###########################################
# function for calling DE genes or sites
test_drug_dependence <- function(cds_subset, 
                                 cds, 
                                 reference_cells=NULL, 
                                 min_fraction_cells=0.05, 
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

##########################################
# load sciPlexRNA data
cds_rna <- readRDS(paste0(basepath, "analysis/archr_revised/scRNA/sciPlex2_cds.RDS"))

# filter cells 
cds_rna <-
  cds_rna[, !is.na(pData(cds_rna)$top_to_second_best_ratio) &
            pData(cds_rna)$top_to_second_best_ratio > 10 &
            pData(cds_rna)$qval < 0.01 &
            pData(cds_rna)$hash_umis > 30]

colData(cds_rna)$cell_type = "A549"

# format colData of RNA cds
meta_data = 
  stringr::str_split_fixed(colData(cds_rna)$top_oligo,
                           pattern = "_",
                           n = 3)

colData(cds_rna)$cell =
  colData(cds_rna)$Cell

colData(cds_rna)$dose =
  meta_data[, 2] %>%
  as.numeric()

colData(cds_rna)$treatment =
  meta_data[, 1] %>%
  as.character()

colData(cds_rna)$culture_well =
  meta_data[, 3] %>%
  stringr::str_sub(start = 2,
                   end = 4) %>%
  as.character()

colData(cds_rna)$well_oligo =
  meta_data[, 3] %>%
  as.character()


colData(cds_rna)$culture_plate =
  meta_data[, 3] %>%
  stringr::str_sub(start = 1,
                   end = 1) %>%
  as.character()

colData(cds_rna)$vehicle =
  colData(cds_rna)$dose  == 0

colData(cds_rna)$dose_character =
  factor(colData(cds_rna)$dose,
         levels = c("0", "0.1", "0.5", "1", "5", "10", "50", "100"))

colData(cds_rna)$new_treatment_label =
  sapply(colData(cds_rna)$treatment, function(x) {
    if (grepl("BMS", x))
      return("BMS345541")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
  })

colData(cds_rna)$solvent =
  sapply(pData(cds_rna)$treatment, function(x) {
    if (grepl("Dex", x))
      return("Ethanol")
    else
      return("DMSO")
  })

colData(cds_rna)$replicate <-
  unlist(sapply(colData(cds_rna)$culture_well, function(x) {
    if (stringr::str_sub(x, -2) %in% c("01", "04", "07", "10"))
      return("rep1")
    if (stringr::str_sub(x, -2) %in% c("02", "05", "08", "11"))
      return("rep2")
    if (stringr::str_sub(x, -2) %in% c("03", "06", "09", "12"))
      return("rep3")
    else
      return(NA)
  }))
colData(cds_rna)$treatment_RNAbased = paste(colData(cds_rna)$treatment, colData(cds_rna)$dose, sep = "_")
colData(cds_rna)$treatmentRNA_broad = ifelse(colData(cds_rna)$dose == 0, "vehicle", colData(cds_rna)$treatment)

# save processed/formatted cds object 
saveRDS(cds_rna, paste0(basepath, "analysis/archr_revised/scRNA/sciPlex2_cds_NB5processed.RDS"))


# DEGene/site analysis 
sciPlex_cds = cds_rna

vehicle_cells = data.frame(pData(sciPlex_cds)) %>% filter(vehicle == TRUE)
vehicle_cells = as.character(vehicle_cells$cell)

reference_cells = sample_n(data.frame(pData(sciPlex_cds)), 1000)$cell

# run test_drug_dependence function on cds (group_by drug).
plan("multiprocess")
old_omp_num_threads = Sys.getenv("OMP_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
simple_model_test_res =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = purrr::map(
      data,
      .f = test_drug_dependence,
      sciPlex_cds,
      union(vehicle_cells, reference_cells),
      residualModelFormulaTerms = NULL,
      min_fraction_cells=0.01,
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

smtr = select(simple_model_test_res, -model, -model_summary)

write.table(smtr, 
            file = "sciPlex_RNA_A549_small_screen_DEGs.txt", 
            row.names=FALSE, 
            quote=FALSE, 
            sep="\t")

simple_model_hit_table = simple_model_test_res %>%
  filter(grepl("dose_", term)) %>%
  mutate(hit = ifelse(q_value < 0.05, 1, 0)) %>%
  dplyr::select(id, term, hit) %>%
  spread(term, hit, fill = 0)

#colnames(simple_model_hit_table) = stringr::str_replace(colnames(adjusted_model_hit_table), "dose_", "")
pdf(file =  "sciPlexRNA_hits_upset.pdf",
    onefile = FALSE,
    height = 3) # or other device
upset(
  as.data.frame(simple_model_hit_table),
  mainbar.y.label = "DEGs in intersection",
  sets.x.label = "Total DEGs"
)
dev.off()

#######################################################
#######################################################
# test drug dependence of peaks from atac_cds

# load ArchR project 
prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered/"),
                       showLogo = FALSE)

# load PeakMatrix into memory
pMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"sc2_filtered/ArrowFiles/sc2.arrow"),
                          ArchRProj = prj, 
                          cellNames = prj$cellNames,
                          useMatrix = "PeakMatrix", 
                          binarize = TRUE)

# Create CDS from peak Matrix (SummarizedExperiment)
rd = data.frame(rowRanges(pMat)) %>% 
  rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end))

# peaks in this list are in a DIFFERENT ORDER 
# relative to pMat (rowRanges above is correct)
pd = data.frame(getPeakSet(prj)) %>% 
  rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end), 
         gene_short_name = ifelse(peakType == "Promoter" & distToGeneStart < 500, nearestGene, "NA")) %>% 
  dplyr::select(rowname, score, replicateScoreQuantile, groupScoreQuantile, 
                Reproducibility, GroupReplicate, distToGeneStart,
                nearestGene, peakType, distToTSS, nearestTSS,
                gene_short_name, GC, idx, peak)

a_rd = left_join(rd, pd, by = "peak") %>% 
  select(rowname = rowname.x, idx = idx.x, score, replicateScoreQuantile, groupScoreQuantile, 
         Reproducibility, GroupReplicate, distToGeneStart,
         nearestGene, peakType, distToTSS, nearestTSS,
         gene_short_name, GC, peak)

row.names(a_rd) <- a_rd$peak
row.names(pMat) <- a_rd$peak
cds_atac = monocle3::new_cell_data_set(assays(pMat)$PeakMatrix, 
                                       cell_metadata = colData(pMat), 
                                       gene_metadata = a_rd)
#format colData of atac_cds

colData(cds_atac)$cell =
  row.names(colData(cds_atac))

colData(cds_atac)$log_nFrags =
  log(colData(cds_atac)$nFrags)

colData(cds_atac)$dose_character =
  factor(colData(cds_atac)$Relative_dose,
         levels = c("0", "0.1", "0.5", "1", "5", "10", "50", "100"))

colData(cds_atac)$new_treatment_label =
  sapply(colData(cds_atac)$treatment, function(x) {
    if (grepl("BMS", x))
      return("BMS345541")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
  })

colData(cds_atac)$solvent =
  sapply(pData(cds_atac)$treatment, function(x) {
    if (grepl("Dex", x))
      return("Ethanol")
    else
      return("DMSO")
  })

###########################################
# DAsite analysis 
sciPlex_cds = cds_atac

vehicle_cells = data.frame(pData(sciPlex_cds)) %>% filter(vehicle == TRUE)
vehicle_cells = as.character(vehicle_cells$cell)

reference_cells = sample_n(data.frame(pData(sciPlex_cds)), 1000)$cell

# run test_drug_dependence function on cds (group_by drug).
plan("multiprocess")
old_omp_num_threads = Sys.getenv("OMP_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 10)
simple_model_test_res =
  colData(sciPlex_cds) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  group_by(cell_type, treatment) %>%
  nest() %>%
  mutate(
    test_res = purrr::map(
      data,
      .f = test_drug_dependence,
      sciPlex_cds,
      union(vehicle_cells, reference_cells),
      expFamily = "binomial",
      residualModelFormulaTerms = c("log_nFrags", "TSSEnrichment"),
      min_fraction_cells=0.01,
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

# Can be downloaded from website found under Supplemental Table 2
## Insert link here
smtr = select(simple_model_test_res, -model, -model_summary)

write.table(smtr, 
            file =  "sciPlex_ATAC_A549_small_screen_DApeaks.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")



simple_model_hit_table = simple_model_test_res %>%
  filter(grepl("dose_", term)) %>%
  mutate(hit = ifelse(q_value < 0.05, 1, 0)) %>%
  dplyr::select(peak, term, hit) %>%
  spread(term, hit, fill = 0)


#colnames(simple_model_hit_table) = stringr::str_replace(colnames(adjusted_model_hit_table), "dose_", "")
pdf(file = "sciPlexATAC_hits_upset.pdf",
    onefile = FALSE,
    height = 3) # or other device
upset(
  as.data.frame(simple_model_hit_table),
  mainbar.y.label = "DA sites in intersection",
  sets.x.label = "Total DA sites"
)
dev.off()



#######################################################
#######################################################
# use Gene Activety scores from ArchR to test drug dependence

# load ArchR project 
prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered/"),
                       showLogo = FALSE)


# load PeakMatrix into memory
gMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"sc2_filtered/ArrowFiles/sc2.arrow"),
                          ArchRProj = prj, 
                          cellNames= prj$cellNames, # restricts to filtered cells
                          useMatrix = "GeneScoreMatrix", 
                          binarize = FALSE)

# Create CDS from peak Matrix (SummarizedExperiment)
g_rd = data.frame(rowData(gMat)) %>% 
  dplyr::mutate(gene_coords = paste0(seqnames, "_", start, "_", end)) %>% 
  dplyr::rename(gene_short_name = name) %>% 
  dplyr::select(gene_coords, gene_short_name, idx)
row.names(g_rd) <- g_rd$gene_short_name
row.names(gMat) <- g_rd$gene_short_name
cds_geneScores = monocle3::new_cell_data_set(assays(gMat)$GeneScoreMatrix, 
                                       cell_metadata = colData(gMat), 
                                       gene_metadata = g_rd)
#format colData of atac_cds

colData(cds_geneScores)$cell =
  row.names(colData(cds_geneScores))

colData(cds_geneScores)$log_nFrags =
  log(colData(cds_geneScores)$nFrags)

colData(cds_geneScores)$dose_character =
  factor(colData(cds_geneScores)$Relative_dose,
         levels = c("0", "0.1", "0.5", "1", "5", "10", "50", "100"))

colData(cds_geneScores)$new_treatment_label =
  sapply(colData(cds_geneScores)$treatment, function(x) {
    if (grepl("BMS", x))
      return("BMS345541")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
  })

colData(cds_geneScores)$solvent =
  sapply(pData(cds_geneScores)$treatment, function(x) {
    if (grepl("Dex", x))
      return("Ethanol")
    else
      return("DMSO")
  })

###########################################
# DAsite analysis 
sciPlex_cds = cds_geneScores

vehicle_cells = data.frame(colData(sciPlex_cds)) %>% 
  dplyr::filter(vehicle == TRUE)
vehicle_cells = as.character(vehicle_cells$cell)

reference_cells = dplyr::sample_n(data.frame(pData(sciPlex_cds)), 1000)$cell

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

## Insert link here
#simple_model_test_res = readRDS("sciPlex_ATAC_A549_small_screen_DAgeneScores")

smtr = select(simple_model_test_res, -model_summary, -model)

write.table(smtr, 
            file = "sciPlex_ATAC_A549_small_screen_DAgeneScores.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")



simple_model_hit_table = simple_model_test_res %>%
  filter(grepl("dose_", term)) %>%
  mutate(hit = ifelse(q_value < 0.05, 1, 0)) %>%
  dplyr::select(gene_short_name, term, hit) %>%
  spread(term, hit, fill = 0)


#colnames(simple_model_hit_table) = stringr::str_replace(colnames(adjusted_model_hit_table), "dose_", "")
pdf(file = "sciPlexATAC_ATACGeneScore_hits_upset.pdf",
    onefile = FALSE,
    height = 3) # or other device
upset(
  as.data.frame(simple_model_hit_table),
  mainbar.y.label = "DA geneScores in intersection",
  sets.x.label = "Total DA geneScores"
)
dev.off()



#################################################
# compare DEG (RNA) coefs vs DAG (ATAC) coefs. 
gene_model_res = read.table(file = "sciPlex_RNA_A549_small_screen_DEGs.txt", head = TRUE)
genescore_model_res = read.table( file = "sciPlex_ATAC_A549_small_screen_DAgeneScores.txt", head = TRUE)

# get all sig genes for all drugs
DEgene_model_res = 
  gene_model_res %>% 
  filter(grepl("dose", term)) %>% 
  select(id, gene_short_name, num_cells_expressed, 
         term, estimate, p_value, q_value) %>% 
  mutate(gene_and_term = paste(gene_short_name, term, sep ="_")) %>% 
  filter(q_value < 0.05)

# isolate gene score estimates for joining to DEGene info
gs_model_res = 
  genescore_model_res %>% 
  filter(grepl("dose", term)) %>% 
  mutate(gene_and_term = paste(gene_short_name, term, sep ="_")) %>% 
  select(gene_and_term, estimate_gs = estimate, q_value_gs = q_value)# %>% 
  #rename(estimate_gs = estimate, q_value_gs = q_value)

RNAvsGenescore_estim = 
  dplyr::left_join(DEgene_model_res, gs_model_res, by = "gene_and_term")

## scatterplots of drug specific coefficients
ggplot(RNAvsGenescore_estim, aes(x =  estimate, y = estimate_gs)) +
  geom_point_rast(size = 2, alpha = 0.5) +
  facet_wrap(~term) +
  xlab("RNA Coefficient") +
  ylab("ATAC Coefficient") +
  stat_cor(method = "pearson", aes(label = ..r.label..),
           label.x = -0.5, label.y = 0.3) +
  geom_smooth(method=lm, se = TRUE) + 
  theme_bw() 
ggsave(filename = "CoeficientScatter_RNAvsGeneScore_DEGs.pdf",
       width = 4, height = 4)



#######################################################
# write DA peak bed file for each drug

DA_model_res = read.table(file = "sciPlex_ATAC_A549_small_screen_DApeaks.txt", head = TRUE)

for (d in c("SAHA", "Dex", "BMS", "Nutlin")){
  sites = DA_model_res %>% 
    filter(term == paste0("dose_", d)) %>% 
           #q_value < 0.05) %>% 
    mutate(chr = stringr::str_split_fixed(peak, "_", 3)[,1],
           start = stringr::str_split_fixed(peak, "_", 3)[,2],
           end = stringr::str_split_fixed(peak, "_", 3)[,3]) %>% 
    select(chr, start, end, term, estimate, q_value)
  write.table(sites, file = paste0(d, "_tested_sites.bed"), 
              col.names = TRUE, 
              row.names = FALSE,
              quote = FALSE, 
              sep = "\t")
  }






