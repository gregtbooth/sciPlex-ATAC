# This notebook repeats much of the same analysis from NB14 but now using the Gene Score Matrix
# Will generate pseudobulk aggregate gene scores and smoothed heatmaps across pseudotime bins 
# focuses only on Allo-activated T-cells

basepath = "/home/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB15"))
setwd(paste0(out_dir, "results/NB15"))

suppressPackageStartupMessages({
  library(ArchR)
  #library(Seurat)
  library(monocle3)
  library(cicero)
  library(devtools)
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
  library(piano)
  plan(multicore)
})

suppressPackageStartupMessages({
  source(paste0(out_dir, "external_resources/bin/GSA_helper_functions.R"))
  source(paste0(out_dir, "external_resources/bin/loadGSCSafe.R"))
})

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered_annotated"), 
                       showLogo = FALSE) 

# Load CDS with T-cell trajectory pseudotime information (then pull only allo-activated cells)
cds = readRDS(file = paste0(out_dir, "results/NB8/cds_NaiveActive_Tcells_PT"))
cds_f = cds[,colData(cds)$cellType_broad == "Activated_Tcell" & colData(cds)$Stimulator == "StimA"]
#cds_a = detect_genes(cds_f)

####################################################
####################################################
# fit models for peaks over pseudotime across Activated Tcells
# Only run to the next double hash once (ouput saved)! 
##########################
# Create a aggregate cds (columns = number of bins)"
# get Gene score matrix and create cds object
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
cds_a <- detect_genes(cds_gs)
cutoff = ncol(cds_a)*.05
cds_a = cds_a[rowData(cds_a)$num_cells_expressed > cutoff, ]

# Redefine Pseudotime_bins (based only on active T-cells, but same pseudotime values)
set.seed(2017)
colData(cds_a)$Pseudotime_bin <- kmeans(colData(cds_a)$Pseudotime, centers=10)$cluster
table(colData(cds_a)$Pseudotime_bin)

# break bins down into 50-cell aggregates
cds_a <- detect_genes(cds_a)
grouping <- plyr::ddply(as.data.frame(colData(cds_a)), plyr::.(Pseudotime_bin), function(x) {
  cells <- x$cell
  if(length(cells) < 100) {
    return(data.frame(cell = x$cell, 
                      num_genes_expressed = x$num_genes_expressed,
                      Pseudotime = x$Pseudotime, 
                      cell_bin = paste(x$Pseudotime_bin, 1, sep="_")))
  } else {
    num_cuts <- floor(length(cells)/50)
    return(data.frame(cell=x$cell, 
                      num_genes_expressed = x$num_genes_expressed,
                      Pseudotime = x$Pseudotime,
                      cell_bin = paste(x$Pseudotime_bin, ntile(x$Pseudotime, n = num_cuts), sep="_")))
  }
})

row.names(grouping) <- grouping$cell
grouping <- grouping[row.names(colData(cds_a)),]
colData(cds_a)$cell_bin <- grouping$cell_bin

table(colData(cds_a)$cell_bin)

# Aggregates cells by Pseudotime bin
binned_cds <-aggregate_by_cell_bin(cds_a, "cell_bin")

group_info <- plyr::ddply(as.data.frame(colData(cds_a)), plyr::.(cell_bin), function(x) {
  data.frame(mean_pseudo = mean(x$Pseudotime), 
             mean_num_genes = mean(x$num_genes_expressed),
             ncells = nrow(x))
})
row.names(group_info) <- group_info$cell_bin
group_info$group_name <- row.names(group_info)

colData(binned_cds) <- merge(colData(binned_cds), group_info)

row.names(colData(binned_cds)) <- colData(binned_cds)$cell_bin
colData(binned_cds) <- colData(binned_cds)[row.names(group_info),]

rowData(binned_cds) <- rowData(cds_gs)[row.names(binned_cds),]

binned_cds <- detect_genes(binned_cds, min_expr=0.1)
binned_cds <- estimate_size_factors(binned_cds)

rowData(binned_cds)$use_for_ordering <- FALSE

# save agrregate CDS with cells binned by Pseudotime bins")
saveRDS(binned_cds, file = "cds_AlloActTcell_GeneScores_pseudotime_agg")

######################
## NOTE: This Script is pretty memory-intensive (saves/reads full models for each peak)
## run on cluster

# run fit_models on aggregated ATAC data by pseudotime with spline

# here counts are estimated to be distributed according to 
# negative binomial b/c counts are no longer binary (cells are binned)
site_fits = fit_models(
  binned_cds, 
  expression_family = "negbinomial", 
  model_formula_str = "~ splines::ns(mean_pseudo, df = 3) + splines::ns(mean_num_genes, df = 3)")

site_fits$site_name = row.names(rowData(binned_cds))
fit_coefs = coefficient_table(site_fits)
pseudotime_terms = fit_coefs %>% 
  dplyr::filter(grepl("mean_pseudo", term)) %>% 
  dplyr::select(-model_summary, -model)

#writing pseudotime models and coefficients terms to file)
saveRDS(site_fits, file  = "fit_GeneScore_models_AlloActTcell_pseudotime_agg")
saveRDS(pseudotime_terms, file = "spline_GeneScore_coefs_AlloActTcell_pseudotime_agg")

###################################################
###################################################

# plot raw and smoothed heatmaps of binned data across pseudotime for significant sites 
# load binned cds
cds_b = readRDS(file = "cds_AlloActTcell_GeneScores_pseudotime_agg")
# load models fit over pseudotime (to binned data)
gene_fits_b = readRDS(file = "fit_GeneScore_models_AlloActTcell_pseudotime_agg")

gene_coefs = read.table(paste0(out_dir,"results/NB13/alloActTcell_Pseudotime_DAgene_coefs.txt"), head = TRUE) %>% 
  filter(term == "Pseudotime")

DA_genes = filter(gene_coefs, grepl("Pseudotime", x = term), p_value < 0.05) %>% 
  distinct(gene_short_name, .keep_all = TRUE)

cds_b_sig <- cds_b[DA_genes$gene_short_name,]

# create a set of pseudotime values to predict peak accessibility from 
bin_num = 20 # set the number of bins for smoothing
pseudoorder <- seq(min(colData(cds_b_sig)$mean_pseudo),
                   max(colData(cds_b_sig)$mean_pseudo),
                   length.out=bin_num) 

############################################
# prepare smoothed heatmaps of DA sites across PT bins
############################################
#get smoothed curves for each peak over pseudotime (calculated using output from previously determined model fits)")
gene_fits_b_sig = dplyr::filter(gene_fits_b, site_name %in% row.names(rowData(cds_b_sig)))

# curves are generated from models for each peak, imputing supplied variables from "new_data" (needs the same names as model input)
curves = model_predictions(gene_fits_b_sig, new_data = data.frame(mean_pseudo=pseudoorder, mean_num_genes = 10000))
rownames(curves) <- gene_fits_b_sig$site_name
rm(gene_fits_b_sig)

#split sites into transient, opening, and closing  
hm_mat <- curves
hm_mat <- hm_mat[!rowSums(!is.finite(hm_mat)),] # remove rows with NAs
hm_mat_scaled <- hm_mat - apply(hm_mat, 1, min)
hm_mat_scaled <- hm_mat_scaled / apply(hm_mat_scaled, 1, max)
hm_mat_scaled_z <- t(scale(t(hm_mat))) # want to scale by row (hence transpose)

# identify transiently opened or closed peaks
temp <- hm_mat_scaled_z < 0 | hm_mat_scaled_z == "NaN"
transient <- apply(temp, 1, function(x) {
  current <- x[1]
  count <- 0
  for(y in x) {
    if (y != current) {
      count <- count + 1
      current <- y
    }
  }
  count
})
trans_max <- apply(hm_mat_scaled, 1, which.max)
trans <- hm_mat_scaled[transient > 1 & !trans_max %in% c(1:ceiling(bin_num/5), ceiling(bin_num-(bin_num/5)):bin_num),]
trans_list <- row.names(trans)

# get genes that open and close over pseudo dose (non transient)
nt_f <- hm_mat_scaled[!row.names(hm_mat_scaled) %in% trans_list,]
nt = nt_f[!is.na(nt_f[,1]),]
up <- nt[nt[,1] < nt[,bin_num],]
down <- nt[nt[,1] > nt[,bin_num],]

#order sites by pseudodose at which half the max accessibility is reached
up_cp <- apply(up, 1, function(x) which.min(abs(x - 0.5)) ) 
up <- up[order(up_cp),]

down_cp <- apply(down, 1, function(x) which.min(abs(x - 0.5)) )
down <- down[order(down_cp),]

all <- rbind(up, down)

if(!is.null(trans_list)){
  trans_cp <- apply(trans, 1, function(x) which.min(abs(x - 0.5)))
  trans <- trans[order(trans_cp),]
  
  all <- rbind(all, trans)
} else {
  trans_cp <- NULL
}

# make a column color bar dataframe for heatmap based on fraction allo_stim in each bin
#bin_df = data.frame("frac_allo" = bin_frac_allo)
#rownames(bin_df) = colnames(all)

# set colors for fraction allo labels
#ann_colors = list(
#  frac_allo = c("white", "firebrick"))

# print smoothed heatmap to file 
options(bitmapType='cairo')

pheatmap::pheatmap(
  all,
  color = colorRampPalette(c("#3C1642",  "#1DD3B0", "#AFFC41"), space = "Lab")(100),
  fontsize = 6,
  scale = "none",
  width=3,
  height=3,
  gaps_row=c(nrow(up), (nrow(up) + nrow(down))),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  #annotation_col = bin_df,
  #annotation_colors = ann_colors,
  filename= "smoothed_PTallo_actTcells_FitGenes_heatmap_20_p05.jpg")

# Print unsmoothed heatmap to file

cell_order <- colData(cds_b_sig)[order(colData(cds_b_sig)$mean_pseudo),]$cell_bin
un_smoothed <- counts(cds_b_sig)[row.names(all), cell_order]
un_smoothed <- t(t(un_smoothed)/colData(cds_b_sig)[cell_order,]$ncells)

pheatmap::pheatmap(
  un_smoothed,
  color = colorRampPalette(c("#3C1642",  "#1DD3B0", "#AFFC41"), space = "Lab")(100),
  fontsize = 6,
  #scale = "row",
  width=2.8,
  height=3.3,
  gaps_row=c(nrow(up), (nrow(up) + nrow(down))),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  filename= "unsmoothed_PTallo_actTcells_FitGenes_heatmap_20_p05.jpg")


##############################
##########################
########## GSEA ##########
##########################

# load DA coefficients comparing Allo vs. Bead Stim. T-cells. 
# DA based on Stimulation type (ST)
# NOTE:  Here positive coefs mean open in Bead relative to Allo

coefs_ST = read.table(paste0(out_dir,"results/NB13/alloActTcell_Pseudotime_DAgene_coefs.txt"), head = TRUE) %>% 
  filter(term == "Pseudotime")
         #p_value < 0.01) # Used to perform without this DA gene filter (04/25/22)


## Load Gene Set Collection
hallmarks <-loadGSCSafe(file= paste0(out_dir,"external_resources/GeneSets/h.all.v6.0.symbols.gmt"))

# This function is designed to collect GSA stats on DAgenes 
get_gsa_res <- function(DAgene_test_res, gsc, min_genes_in_set=5, ncores=1){
  DAgene_test_res = DAgene_test_res %>% 
    dplyr::group_by(gene_short_name) %>%
    dplyr::filter(abs(test_val) == max(abs(test_val))) 
  directions=data.frame(row.names=as.character(DAgene_test_res$gene_short_name), stats=DAgene_test_res$test_val)
  geneLevelStats=data.frame(row.names=as.character(DAgene_test_res$gene_short_name), stats=DAgene_test_res$test_val)
  tryCatch({
    gsares <- runGSA(geneLevelStats=geneLevelStats, direction=geneLevelStats, geneSetStat="mean", gsc=gsc, nPerm=10000, ncpu=ncores, gsSizeLim=c(min_genes_in_set, Inf), verbose=T)
    return (gsares)
  }, error = function(e) { print(e); return (NA) } )
}
Sys.setenv(OMP_NUM_THREADS = 1)


# run GSEA using DA genes (comparing Allo vs Bead directly)
cat("running GSEA analysis for DAgenes across PT for allo-activated T-cells")
  
DAgene_gsa_res = filter(coefs_ST, term == "Pseudotime") %>% 
  nest(-term) %>% 
  dplyr::mutate(
    gsares = purrr::map(.f = get_gsa_res, .x = data, hallmarks, ncores=4)
  ) 
  
# Extract p values and GSA test statistics into a table
DAgene_gsa_summary = 
  DAgene_gsa_res %>%
  dplyr::mutate(gsa_summary = purrr::map(.f = GSAsummaryTable, .x = gsares)) %>% 
  unnest(gsa_summary)


## Make a heatmap showing the effects of all the drugs on one cell line:
simple_model_gsa_summary = DAgene_gsa_summary %>%
  dplyr::mutate(up_stat = `Stat (mix.dir.up)`, dn_stat = `Stat (mix.dir.dn)`) %>%
  dplyr::select(term,
                Name,
                 `Genes (tot)`,
                 `p (mix.dir.up)`,
                up_stat,
                `p (mix.dir.dn)`,
                dn_stat)
simple_model_gsa_summary = simple_model_gsa_summary %>% 
  dplyr::mutate(gsa_stat = ifelse((is.na(up_stat) == FALSE & `p (mix.dir.up)` < `p (mix.dir.dn)`),
                                  up_stat, ifelse(is.na(dn_stat) == FALSE, -dn_stat, 0)
  ))
# save summary output from GSEA 
write.table(simple_model_gsa_summary, file = paste0("AlloActivatedTcells_Pseudotime_Response_GSEAsummary_p01_DAgenesOnly.txt"),
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 
  
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
  dplyr::mutate(Name = stringr::str_split_fixed(Name, "%", 3)[, 1], 
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

ggsave(file = paste0("AlloActivatedTcells_Pseudotime_Response_GSEA_barplot_p01_DAgenesOnly.pdf"), width = 5, height = 3)

