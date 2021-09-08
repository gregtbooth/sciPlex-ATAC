## NOTE: This Script is pretty memory-intensive (saves/reads full models for each peak)
## run on cluster

basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB9"))
setwd(paste0(out_dir, "results/NB9"))

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
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)


###########################################################
# Create a aggregate cds (columns = number of bins)"
###########################################################
# only need to run once (output saved)

# load processed cds with processed naive to active T-cell trajectory 
cds = readRDS(file = paste0(out_dir, "results/NB7/cds_NaiveActive_Tcells_PT"))
cds <- detect_genes(cds)

# subset to peaks open in > 0.1% cells
cutoff = ncol(cds)*.001
cds = cds[rowData(cds)$num_cells_expressed > cutoff, ]


table(colData(cds)$Pseudotime_bin)

# break bins down in to 50-cell aggregates
grouping <- plyr::ddply(as.data.frame(colData(cds)), plyr::.(Pseudotime_bin), function(x) {
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
grouping <- grouping[row.names(colData(cds)),]
colData(cds)$cell_bin <- grouping$cell_bin

table(colData(cds)$cell_bin)

# Aggregates cells by Pseudotime bin
binned_cds <-aggregate_by_cell_bin(cds, "cell_bin")

group_info <- plyr::ddply(as.data.frame(colData(cds)), plyr::.(cell_bin), function(x) {
  data.frame(mean_pseudo = mean(x$Pseudotime), 
             mean_num_genes = mean(x$num_genes_expressed),
             ncells = nrow(x))
})
row.names(group_info) <- group_info$cell_bin
group_info$group_name <- row.names(group_info)

colData(binned_cds) <- merge(colData(binned_cds), group_info)

row.names(colData(binned_cds)) <- colData(binned_cds)$cell_bin
colData(binned_cds) <- colData(binned_cds)[row.names(group_info),]

rowData(binned_cds) <- rowData(cds)[row.names(binned_cds),]

binned_cds <- detect_genes(binned_cds, min_expr=0.1)
binned_cds <- estimate_size_factors(binned_cds)

rowData(binned_cds)$use_for_ordering <- FALSE

# save agrregate CDS with cells binned by Pseudotime bins")
saveRDS(binned_cds, file = "cds_Tcell_pseudotime_agg")

######################
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
saveRDS(site_fits, file  = "fit_models_Tcell_pseudotime_agg")
saveRDS(pseudotime_terms, file = "spline_coefs_Tcell_pseudotime_agg")

###################################################
###################################################

# plot raw and smoothed heatmaps of binned data across pseudotime for significant sites 
# load binned cds
cds_b = readRDS(file = "cds_Tcell_pseudotime_agg")
# load models fit over pseudotime (to binned data)
gene_fits_b = readRDS(file = "fit_models_Tcell_pseudotime_agg")

# load regression coefficients 
# note These are the DA sites called using individual cells over Pseudotime (with spline fit)
peak_coefs = read.table(paste0(out_dir, "results/NB8/Tcell_Activation_PT_DApeak_spline_coefs.txt"),
                        sep = "\t", head = TRUE, stringsAsFactors = FALSE)
  
DA_sites = filter(peak_coefs, grepl("Pseudotime", x = term), q_value < 0.05) %>% 
  distinct(peak, .keep_all = TRUE)
  
cds_b_sig <- cds_b[DA_sites$peak,]
pseudoorder <- seq(min(colData(cds_b_sig)$mean_pseudo),
                   max(colData(cds_b_sig)$mean_pseudo),
                   length.out=100)
  
#get smoothed curves for each peak over pseudotime (calculated using output from previously determined model fits)")
gene_fits_b_sig = dplyr::filter(gene_fits_b, site_name %in% row.names(rowData(cds_b_sig)))

# curves are generated from models for each peak, imputing supplied variables from "new_data" (needs the same names as model input)
curves = model_predictions(gene_fits_b_sig, new_data = data.frame(mean_pseudo=pseudoorder, mean_num_genes = 10000))
rownames(curves) <- gene_fits_b_sig$site_name
rm(gene_fits_b)
rm(gene_fits_b_sig)
  
#split sites into transient, opening, and closing  
hm_mat <- curves
hm_mat_scaled <- hm_mat - apply(hm_mat, 1, min)
hm_mat_scaled <- hm_mat_scaled / apply(hm_mat_scaled, 1, max)
hm_mat_scaled_z <- t(scale(t(hm_mat)))

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
trans <- hm_mat_scaled[transient > 1 & !trans_max %in% c(1:20, 80:100),]
trans_list <- row.names(trans)
  
# get genes that open and close over pseudo dose (non transient)
nt_f <- hm_mat_scaled[!row.names(hm_mat_scaled) %in% trans_list,]
nt = nt_f[!is.na(nt_f[,1]),]
up <- nt[nt[,1] < nt[,100],]
down <- nt[nt[,1] > nt[,100],]
  
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
  
  
# print smoothed heatmap to file 
options(bitmapType='cairo')
  
pheatmap::pheatmap(
  all,
  color = colorRampPalette(c("#3C1642",  "#1DD3B0", "#AFFC41"), space = "Lab")(100),
  fontsize = 6,
  scale = "none",
  width=2,
  height=3,
  gaps_row=c(nrow(up), (nrow(up) + nrow(down))),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_colnames = FALSE,
  show_rownames = FALSE,
  filename= "smoothed_PT_Fit_heatmap.jpg")
  
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
  filename= "unsmoothed_PT_heatmap.jpg")
  
# Add column to coeficient table with direction of DA changes
DA_sites$direction <- sapply(DA_sites$site_name, function(x) {
  if(x %in% names(up_cp)) return("Opening")
  if(x %in% names(down_cp)) return("Closing")
  if(x %in% names(trans_cp)) return("Transient")
  return("Static")
})
  
saveRDS(DA_sites, file = "spline_coefs_Tcell_PT_directions")








