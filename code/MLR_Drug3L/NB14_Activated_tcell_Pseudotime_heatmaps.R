basepath = "/home/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB14"))
setwd(paste0(out_dir, "results/NB14"))

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
  library(forcats)
  plan(multicore)
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
cds_a = detect_genes(cds_f)

####################################################
####################################################
# fit models for peaks over pseudotime across Activated Tcells
# Only run to the next double hash once (ouput saved)! 
##########################
# Create a aggregate cds (columns = number of bins)"

# Use CDS (defined above) only containing activated T-cells
cds_a = detect_genes(cds_f)
cutoff = ncol(cds_a)*.005
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

rowData(binned_cds) <- rowData(cds)[row.names(binned_cds),]

binned_cds <- detect_genes(binned_cds, min_expr=0.1)
binned_cds <- estimate_size_factors(binned_cds)

rowData(binned_cds)$use_for_ordering <- FALSE

# save agrregate CDS with cells binned by Pseudotime bins")
saveRDS(binned_cds, file = "cds_AlloActTcell_pseudotime_agg")

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
saveRDS(site_fits, file  = "fit_models_AlloActTcell_pseudotime_agg")
saveRDS(pseudotime_terms, file = "spline_coefs_AlloActTcell_pseudotime_agg")

###################################################
###################################################

# plot raw and smoothed heatmaps of binned data across pseudotime for significant sites 
# load binned cds
cds_b = readRDS(file = "cds_AlloActTcell_pseudotime_agg")
# load models fit over pseudotime (to binned data)
gene_fits_b = readRDS(file = "fit_models_AlloActTcell_pseudotime_agg")

# load regression coefficients 
# note These are the DA sites called (from above) using individual cells over Pseudotime
peak_coefs = read.table(paste0(out_dir,"results/NB13/alloActTcell_Pseudotime_DApeak_coefs.txt"), head = TRUE) %>% 
  filter(term == "Pseudotime")

DA_sites = filter(peak_coefs, grepl("Pseudotime", x = term), p_value < 0.01) %>% 
  distinct(peak, .keep_all = TRUE)

cds_b_sig <- cds_b[DA_sites$peak,]

# create a set of pseudotime values to predict peak accessibility from 
bin_num = 20 # set the number of bins for smoothing
pseudoorder <- seq(min(colData(cds_b_sig)$mean_pseudo),
                   max(colData(cds_b_sig)$mean_pseudo),
                   length.out=bin_num) 

####################################################################
# function retrieves the fraction of Drug Treated (any dose) cells in each PT bin
get_frac_treated = function(cds = cds_a, Pseudoorder = pseudoorder, drug= "SAHA"){
  cdat = data.frame(colData(cds_a)) %>% 
    filter(Drug == drug) %>%
    mutate(vehicle = ifelse(Relative_Dose > 0, FALSE, TRUE))
  res = c()
  for (bin in 1:(length(Pseudoorder))){
    if (bin == 1){
      bin_low = 0
      bin_high = Pseudoorder[bin]
    }
    else{
      bin_low = Pseudoorder[bin-1]
      bin_high = Pseudoorder[bin]
    }
    cd = filter(cdat, Pseudotime >= bin_low & Pseudotime < bin_high) %>% 
      group_by(vehicle) %>% 
      summarise(count = n())
    if(nrow(cd) > 0){
      Untreated_cells = cd[cd$vehicle == TRUE, 2]
      Untreated_cells = ifelse(nrow(Untreated_cells) > 0, Untreated_cells$count, 0)
      Treated_cells = cd[cd$vehicle == FALSE, 2]
      Treated_cells = ifelse(nrow(Treated_cells) > 0, Treated_cells$count, 0)
      frac_Treated = Treated_cells/(Untreated_cells+Treated_cells)}
    else{
      frac_Treated = NA
    }
    res = c(res, frac_Treated)
  }
  res
}

# create data frame of treated cell/fraction per bin (separate column for each drug)
df = data.frame(PTbin = 1:(bin_num))
# calculate fraction of drugged cells per PT bin (separately for each drug)
drug_list = data.frame(colData(cds_a)) %>% 
  distinct(Drug)

for(d in 1:length(drug_list$Drug)){ 
  DRUG <- drug_list$Drug[d]
  cat("counting fraction of treated cells per PT bin for :" ,DRUG, "\n")
  fracCellsTreated = get_frac_treated(cds = cds_a, Pseudoorder = pseudoorder, drug= DRUG)
  #names(fracCellsTreated) <- DRUG
  df[DRUG] = fracCellsTreated# <- cbind(df,fracCellsTreated)
}

# reformat (gather) for ggplot plotting: 
df_gather = df %>% 
  gather(Drug, fracTreated, -(PTbin))

# use the sum of fraction treated across PT bins to order Drugs.
order = df_gather %>% 
  mutate(weight = ifelse(is.na(fracTreated), 0, fracTreated)) %>% 
  group_by(Drug) %>% 
  summarise(imbalance = sum(weight)) %>% 
  arrange(imbalance)

df_gather = left_join(df_gather, order, by = "Drug") %>% 
  mutate(Drug = fct_reorder(Drug, imbalance), 
         Drug_PTbin = paste0(Drug, "_", PTbin))

# plot heatmap of Treated cell fractions across PT bins 
pdf(file = paste0("HeatMap_percent_treated_cells_Recovered_acrossPTbins_20_new.pdf"), 
    width = 7, height = 3.5)
print(  
  ggplot(df_gather, aes(x = PTbin, y = Drug, fill = fracTreated)) + 
    #facet_wrap(~Stimulator) +
    geom_tile(color = "gray") +
    scale_fill_gradient(low = "#625a94", high = "#f9ad2a") + 
    ylab("Drug") +
    xlab("PT bin")  
    #theme_void()
)
dev.off()

#############################################
# function retrieves the Mean DOSE (relative) of Drug Treated cells in each PT bin
get_meanDose_treated = function(cds = cds_a, Pseudoorder = pseudoorder, drug= "SAHA"){
  cdat = data.frame(colData(cds_a)) %>% 
    filter(Drug == drug) 
  res = c()
  for (bin in 1:(length(Pseudoorder))){
    if (bin == 1){
      bin_low = 0
      bin_high = Pseudoorder[bin]
    }
    else{
      bin_low = Pseudoorder[bin-1]
      bin_high = Pseudoorder[bin]
      }
    cd = filter(cdat, Pseudotime >= bin_low & Pseudotime < bin_high)
    if(nrow(cd) > 0){
      mean_rel_dose = summarise(cd, mean_rel_dose = mean(as.numeric(Relative_Dose)))
      mrd = as.vector(mean_rel_dose)
    }
    else{
      mean_rel_dose = NA
      mrd = as.vector(mean_rel_dose)
    }
    res = c(res, mrd)
  }
  res_v = unname(unlist(res))
  res_v
}

# create data frame of treated cell/fraction per bin (separate column for each drug)
df_mean = data.frame(PTbin = 1:(bin_num))
# calculate fraction of drugged cells per PT bin (separately for each drug)
drug_list = data.frame(colData(cds_a)) %>% 
  distinct(Drug)

for(d in 1:length(drug_list$Drug)){ 
  DRUG <- drug_list$Drug[d]
  cat("Getting mean dose of cells per PT bin for :" ,DRUG, "\n")
  PTbinMeanDose = get_meanDose_treated(cds = cds_a, Pseudoorder = pseudoorder, drug= DRUG)
  #names(fracCellsTreated) <- DRUG
  df_mean[DRUG] = PTbinMeanDose# <- cbind(df,fracCellsTreated)
}

# reformat (gather) for ggplot plotting: 
df_gather_md = df_mean %>% 
  gather(Drug, meanDose, -(PTbin)) %>% 
  mutate(Drug_PTbin = paste0(Drug, "_", PTbin)) %>% 
  select(Drug_PTbin, meanDose) %>% 
  left_join(df_gather , by = "Drug_PTbin") # add ordering info from fraction treated plot

# plot heatmap of Treated cell fractions across PT bins 
pdf(file = paste0("HeatMap_MeanDose_cellsRecovered_acrossPTbins_20_new.pdf"), 
    width = 7, height = 3.5)
print(  
  ggplot(df_gather_md, aes(x = PTbin, y = Drug, fill = meanDose)) + 
    #facet_wrap(~Stimulator) +
    geom_tile(color = "gray") +
    scale_fill_gradient(low = "#625a94", high = "#f9ad2a") + 
    ylab("Drug") +
    xlab("PT bin") 
)
dev.off()

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
  filename= "smoothed_PTallo_actTcells_Fit_heatmap_20_p01.jpg")

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
  filename= "unsmoothed_PTallo_actTcells_heatmap_20_p01.jpg")


##############################
# using direction assignments, 
# look for enriched motifs/ TFBS chip
#############################

# Add column to coeficient table with direction of DA changes
peak_coefs$direction <- sapply(peak_coefs$peak, function(x) {
  if(x %in% names(up_cp)) return("Opening")
  if(x %in% names(down_cp)) return("Closing")
  if(x %in% names(trans_cp)) return("Transient")
  return("Static")
})


# function to retrieve and format peak by motif matrix
pull_peak_by_motif_mat = function(archr_project){
  motif_dat = readRDS(paste0(getOutputDirectory(archr_project), "/Annotations/Motif-In-Peaks-Summary.rds"))
  motif_mat = assay(motif_dat$motifMatches)
  peaks = data.frame(rowRanges(motif_dat$motifMatches)) %>% # GB Corrected 11/17/2020
    dplyr::mutate(peak = paste(seqnames, start, end, sep = "_"))
  row.names(motif_mat) <- peaks$peak
  return(motif_mat)
}

# Look for TF motifs which are significantly associated with opening, closing, transient or static sites. 
find_sig_motifs = function(site_coef_directions, archr_project){
  cat("load motif matrix.\n") 
  motif_mat = pull_peak_by_motif_mat(archr_project)
  # filter for only used peaks
  motif_mat_f = motif_mat[site_coef_directions$site_name,]
  cat("correcting motif names")
  m = colnames(motif_mat_f)
  m_fix = stringr::str_replace_all(m, 
                                   pattern = "::",
                                   replacement =  "_") %>% 
    stringr::str_replace_all(pattern = "[-.]",
                             replacement = "_") %>% 
    stringr::str_replace_all(pattern = "[()]",
                             replacement = "")
  colnames(motif_mat_f) <- m_fix
  motifs = colnames(motif_mat_f)
  # only consider motifs present in 0.5% of used peaks
  #motif_mat_f = motif_mat_f[, Matrix::colSums(motif_mat_f) > 0.005*nrow(motif_mat_f)]
  
  # reformat from logical to numeric
  mmf = matrix(as.numeric(motif_mat_f), nrow = nrow(x=motif_mat_f))
  row.names(mmf) = row.names(motif_mat_f)
  colnames(mmf) = colnames(motif_mat_f)
  motif_df = as.data.frame(mmf)
  motif_df$site_name = row.names(motif_df)
  # add direction info from coef_table
  motif_df = inner_join(site_coef_directions, motif_df, by = "site_name")
  # add binary columns describing directionality of each site 
  motif_df = cbind(motif_df, as.data.frame(model.matrix(~ 0 + direction, motif_df)))
  # define sites as Distal or Promoter
  motif_df = mutate(motif_df, type = ifelse(distToTSS > 500, "Distal", "Promoter"))
  cat("fit linear regression predicting DA direction from motif presence. \nThis will take a few minutes...\n")
  # For every motif within DA peaks, fit a linear regression to determine 
  # if it is significantly associated with opening, closing, or static peaks
  get.peak.DA.logistic.reg.stats = function(response, features, data) {
    t(sapply(features, function(feature) {
      res = summary(glm(paste(response, feature, sep=" ~ "), family = "binomial", data = data))
      if (nrow(res$coefficients) > 1) {
        return(c(res$coefficients[2,1], res$coefficients[2,4]))
      } else {
        return(c(NA, NA))
      }
    }))
  }
  
  cat("running model on promoter proximal sites...\n") 
  peak.DA.stats = list()
  # check if there are any sites that actually meet criteria (avoids errors)
  if(filter(motif_df, type == "Promoter", direction == "Opening") %>% nrow > 0){
    peak.DA.stats$promoter.opening = get.peak.DA.logistic.reg.stats(
      "directionOpening", motifs, motif_df %>% filter(type == "Promoter"))
    #filter for only the significant motifs
    peak.DA.stats$promoter.opening = data.frame(
      motif = row.names(peak.DA.stats$promoter.opening),
      beta = peak.DA.stats$promoter.opening[,1],
      p.val = peak.DA.stats$promoter.opening[,2]) %>%
      dplyr::filter(!is.na(p.val), beta > 0) %>% dplyr::arrange(p.val)
  }
  if(filter(motif_df, type == "Promoter", direction == "Closing") %>% nrow > 0){
    peak.DA.stats$promoter.closing = get.peak.DA.logistic.reg.stats(
      "directionClosing", motifs, motif_df %>% filter(type == "Promoter"))
    #filter for only the significant motifs
    peak.DA.stats$promoter.closing = data.frame(
      motif = row.names(peak.DA.stats$promoter.closing),
      beta = peak.DA.stats$promoter.closing[,1],
      p.val = peak.DA.stats$promoter.closing[,2]) %>%
      dplyr::filter(!is.na(p.val), beta > 0) %>% dplyr::arrange(p.val)
  }
  if(filter(motif_df, type == "Promoter", direction == "Transient") %>% nrow > 0){
    peak.DA.stats$promoter.transient = get.peak.DA.logistic.reg.stats(
      "directionTransient", motifs, motif_df %>% filter(type == "Promoter"))
    #filter for only the significant motifs
    peak.DA.stats$promoter.transient = data.frame(
      motif = row.names(peak.DA.stats$promoter.transient),
      beta = peak.DA.stats$promoter.transient[,1],
      p.val = peak.DA.stats$promoter.transient[,2]) %>%
      dplyr::filter(!is.na(p.val), beta > 0) %>% dplyr::arrange(p.val)
  }
  if(filter(motif_df, type == "Promoter", direction == "Static") %>% nrow > 0){
    peak.DA.stats$promoter.static = get.peak.DA.logistic.reg.stats(
      "directionStatic", motifs, motif_df %>% filter(type == "Promoter"))
    #filter for only the significant motifs
    peak.DA.stats$promoter.static = data.frame(
      motif = row.names(peak.DA.stats$promoter.static),
      beta = peak.DA.stats$promoter.static[,1],
      p.val = peak.DA.stats$promoter.static[,2]) %>%
      dplyr::filter(!is.na(p.val), beta > 0) %>% dplyr::arrange(p.val)
  }
  cat("running model on distal sites...\n") 
  if(filter(motif_df, type == "Distal", direction == "Opening") %>% nrow > 0){
    peak.DA.stats$distal.opening = get.peak.DA.logistic.reg.stats(
      "directionOpening", motifs, motif_df %>% filter(type == "Distal"))
    #filter for only the significant motifs
    peak.DA.stats$distal.opening = data.frame(
      motif = row.names(peak.DA.stats$distal.opening),
      beta = peak.DA.stats$distal.opening[,1],
      p.val = peak.DA.stats$distal.opening[,2]) %>%
      dplyr::filter(!is.na(p.val), beta > 0) %>% dplyr::arrange(p.val)
  }
  if(filter(motif_df, type == "Distal", direction == "Closing") %>% nrow > 0){
    peak.DA.stats$distal.closing = get.peak.DA.logistic.reg.stats(
      "directionClosing", motifs, motif_df %>% filter(type == "Distal"))
    #filter for only the significant motifs
    peak.DA.stats$distal.closing = data.frame(
      motif = row.names(peak.DA.stats$distal.closing),
      beta = peak.DA.stats$distal.closing[,1],
      p.val = peak.DA.stats$distal.closing[,2]) %>%
      dplyr::filter(!is.na(p.val), beta > 0) %>% dplyr::arrange(p.val)
  }
  if(filter(motif_df, type == "Distal", direction == "Transient") %>% nrow > 0){
    peak.DA.stats$distal.transient = get.peak.DA.logistic.reg.stats(
      "directionTransient", motifs, motif_df %>% filter(type == "Distal"))
    #filter for only the significant motifs
    peak.DA.stats$distal.transient = data.frame(
      motif = row.names(peak.DA.stats$distal.transient),
      beta = peak.DA.stats$distal.transient[,1],
      p.val = peak.DA.stats$distal.transient[,2]) %>%
      dplyr::filter(!is.na(p.val), beta > 0) %>% dplyr::arrange(p.val)
  }
  if(filter(motif_df, type == "Distal", direction == "Static") %>% nrow > 0){
    peak.DA.stats$distal.static = get.peak.DA.logistic.reg.stats(
      "directionStatic", motifs, motif_df %>% filter(type == "Distal"))
    #filter for only the significant motifs
    peak.DA.stats$distal.static = data.frame(
      motif = row.names(peak.DA.stats$distal.static),
      beta = peak.DA.stats$distal.static[,1],
      p.val = peak.DA.stats$distal.static[,2]) %>%
      dplyr::filter(!is.na(p.val), beta > 0) %>% dplyr::arrange(p.val)
  }
  
  cat("returning motif regression data\n")
  return(peak.DA.stats)
}


plot_top_motifs = function(peak.DA.stats, Feat_type = "distal"){
  if(Feat_type == "distal"){
    D.peak.enrich = rbind(
      peak.DA.stats$distal.opening %>% mutate(facet = "Opening"),
      peak.DA.stats$distal.closing %>% mutate(facet = "Closing"),
      peak.DA.stats$distal.static %>% mutate(facet = "Static"))
  }
  if(Feat_type == "promoter"){
    D.peak.enrich = rbind(
      peak.DA.stats$promoter.opening %>% mutate(facet = "Opening"),
      peak.DA.stats$promoter.closing %>% mutate(facet = "Closing"),
      peak.DA.stats$promoter.static %>% mutate(facet = "Static"))
  }
  if(!(Feat_type %in% c("promoter", "distal"))){cat("error: Feat_type not promoter or distal")}
  D.peak.enrich = D.peak.enrich %>% 
    dplyr::group_by(facet) %>% 
    dplyr::arrange(facet, log10(p.val)) %>%
    dplyr::mutate(rank = rank(log10(p.val), ties.method = "first")) %>% 
    dplyr::filter(rank <= 6, p.val < .05) %>% #1e-2
    ungroup() %>%
    dplyr::mutate(motif = as.character(motif))
  
  D.peak.enrich = distinct(D.peak.enrich, motif, .keep_all = T) #FIXIT: removes motifs common to openning and closing. 
  D.peak.enrich$motif_id = factor(D.peak.enrich$motif, levels = rev(D.peak.enrich$motif))
  D.peak.enrich$facet = factor(D.peak.enrich$facet, levels = c("Opening", "Closing", "Static"))
  D.peak.enrich <- D.peak.enrich[order(D.peak.enrich$facet, rev(D.peak.enrich$beta)),]
  
  pdf(paste0(Feat_type,".peak.motif.enrichments_beta.pdf"), width = 2, height=2.4)
  print(
    ggplot(D.peak.enrich, aes(x = motif_id, y = (beta), fill = facet)) +
      geom_bar(stat="identity") +
      coord_flip() +
      scale_fill_manual(values=c("#EF5B5B", "#0FA3B1", "#8ABF69")) +
      xlab("Motif") +
      ylab("Regression Coefficient") +
      guides(fill = guide_legend(title = "Accessibility\nTrend")) +
      monocle3:::monocle_theme_opts() +
      theme(text = element_text(size=6)) +
      theme(legend.position="none") +
      theme(#axis.title.x = element_text(margin=margin(3, 0, 0, 0)),
        #axis.title.y = element_text(margin=margin(0, 3, 0, 0)),
        legend.key.width=unit(0.1, "in"),
        legend.key.height=unit(0.1, "in"))
    #legend.margin=margin(0, 0, 0, -8))
  )
  dev.off()
}

# run the functions above 
# load ArchR project
prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered_annotated"), 
                       showLogo = FALSE) 

peak_dat = getPeakSet(prj) %>% 
  data.frame() %>% 
  mutate(site_name = paste(seqnames, start, end, sep = "_")) %>% 
  dplyr::select(site_name, distToTSS, nearestGene, nearestTSS, GC)

DA.stats = peak_coefs %>% 
  mutate(site_name = peak) %>% 
  left_join(peak_dat, by = "site_name")

write.table(DA.stats, file = "alloActTcell_Pseudotime_DApeak_coefs_directions.txt", 
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
site_motif_stats = find_sig_motifs(site_coef_directions = DA.stats, archr_project = prj)
plot_top_motifs(peak.DA.stats = site_motif_stats, Feat_type = "promoter")
plot_top_motifs(peak.DA.stats = site_motif_stats, Feat_type = "distal")


####################################################
#Define early/late changes as sites where 50% max reached
# look for motif enrichments in early/late opening/close
####################################################
prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered_annotated"), 
                       showLogo = FALSE) 


# Add column to coeficient table with direction of DA changes
peak_coefs$direction <- sapply(peak_coefs$peak, function(x) {
  if(x %in% names(up_cp)) return("Opening")
  if(x %in% names(down_cp)) return("Closing")
  if(x %in% names(trans_cp)) return("Transient")
  return("Static")
})


#before/after middle pseudotime (bin 10)
early_up = up_cp[up_cp < 10] 
late_up = up_cp[up_cp >= 10] 
early_down = down_cp[down_cp < 10] 
late_down = down_cp[down_cp >= 10] 

peak_coefs$timing <- sapply(peak_coefs$peak, function(x) {
  if(x %in% names(early_up)) return("early")
  if(x %in% names(late_up)) return("late")
  if(x %in% names(early_down)) return("early")
  if(x %in% names(late_down)) return("late")
  if(x %in% names(trans_cp)) return("Transient")
  return("Static")
})

peak_coefs = mutate(peak_coefs, direction_timing = paste0(direction, "_", timing))


# pull peak motif data for enrichment analysis. 
# function to retrieve and format peak by motif matrix
pull_peak_by_motif_mat = function(archr_project){
  motif_dat = readRDS(paste0(getOutputDirectory(archr_project), "/Annotations/Motif-In-Peaks-Summary.rds"))
  motif_mat = assay(motif_dat$motifMatches)
  peaks = data.frame(rowRanges(motif_dat$motifMatches)) %>% # GB Corrected 11/17/2020
    dplyr::mutate(peak = paste(seqnames, start, end, sep = "_"))
  row.names(motif_mat) <- peaks$peak
  return(motif_mat)
}

# function regresses  motif presence within accessible site and opening or closing assignment
# For every motif within DA peaks, fit a linear regression to determine 
# if it is significantly associated with opening, closing, or static peaks
get.peak.DA.logistic.reg.stats = function(response, features, data) {
  t(sapply(features, function(feature) {
    res = summary(glm(paste(response, feature, sep=" ~ "), family = "binomial", data = data))
    if (nrow(res$coefficients) > 1) {
      return(c(res$coefficients[2,1], res$coefficients[2,4]))
    } else {
      return(c(NA, NA))
    }
  }))
}


# get peak x motif matrix 
motif_mat = pull_peak_by_motif_mat(prj)
  
# filter for only used peaks
motif_mat_f = motif_mat[peak_coefs$peak,]
  
# list of motif names
motifs = colnames(motif_mat_f)
  
# reformat from logical to numeric
mmf = matrix(as.numeric(motif_mat_f), nrow = nrow(x=motif_mat_f))
row.names(mmf) = row.names(motif_mat_f)
colnames(mmf) = colnames(motif_mat_f)
motif_df = as.data.frame(mmf)
motif_df$peak = row.names(motif_df)
# add direction info from coef_table
motif_df = dplyr::inner_join(peak_coefs, motif_df, by = "peak")
# add binary columns describing directionality of each site 
motif_df = cbind(motif_df, as.data.frame(model.matrix(~ 0 + direction_timing, motif_df)))

cat("fit linear regression predicting DA direction from motif presence. \nThis will take a few minutes...\n")

# set container
peak.DA.stats = list()

if(dplyr::filter(motif_df, direction_timing == "Opening_early") %>% nrow > 0){
  peak.DA.stats$promoter.opening.early = get.peak.DA.logistic.reg.stats(
    "direction_timingOpening_early", motifs, motif_df)
    #filter for only the significant motifs
    peak.DA.stats$promoter.opening.early = data.frame(
    motif = row.names(peak.DA.stats$promoter.opening.early),
    beta = peak.DA.stats$promoter.opening.early[,1],
    p.val = peak.DA.stats$promoter.opening.early[,2]) %>%
    mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
    dplyr::filter(!is.na(Padj), beta > 0) %>% 
    dplyr::arrange(Padj) %>% 
    dplyr::mutate(rank = rank(Padj), 
                  sig = ifelse(Padj <0.05, 1, 0), 
                  mlog10Padj = -log10(Padj)) 
    }

if(dplyr::filter(motif_df, direction_timing == "Opening_late") %>% nrow > 0){
  peak.DA.stats$promoter.opening.late = get.peak.DA.logistic.reg.stats(
    "direction_timingOpening_late", motifs, motif_df)
  #filter for only the significant motifs
  peak.DA.stats$promoter.opening.late = data.frame(
    motif = row.names(peak.DA.stats$promoter.opening.late),
    beta = peak.DA.stats$promoter.opening.late[,1],
    p.val = peak.DA.stats$promoter.opening.late[,2]) %>%
    mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
    dplyr::filter(!is.na(Padj), beta > 0) %>% 
    dplyr::arrange(Padj) %>% 
    dplyr::mutate(rank = rank(Padj), 
                  sig = ifelse(Padj <0.05, 1, 0), 
                  mlog10Padj = -log10(Padj)) 
  }

  
if(dplyr::filter(motif_df, direction_timing == "Closing_early") %>% nrow > 0){
    peak.DA.stats$promoter.closing.early = get.peak.DA.logistic.reg.stats(
    "direction_timingClosing_early", motifs, motif_df)
    #filter for only the significant motifs
    peak.DA.stats$promoter.closing.early = data.frame(
      motif = row.names(peak.DA.stats$promoter.closing.early),
      beta = peak.DA.stats$promoter.closing.early[,1],
      p.val = peak.DA.stats$promoter.closing.early[,2]) %>%
      mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
      dplyr::filter(!is.na(Padj), beta > 0) %>% 
      dplyr::arrange(Padj) %>% 
      dplyr::mutate(rank = rank(Padj), 
                    sig = ifelse(Padj <0.05, 1, 0), 
                    mlog10Padj = -log10(Padj)) 
    }

if(dplyr::filter(motif_df, direction_timing == "Closing_late") %>% nrow > 0){
  peak.DA.stats$promoter.closing.late = get.peak.DA.logistic.reg.stats(
    "direction_timingClosing_late", motifs, motif_df)
  #filter for only the significant motifs
  peak.DA.stats$promoter.closing.late = data.frame(
    motif = row.names(peak.DA.stats$promoter.closing.late),
    beta = peak.DA.stats$promoter.closing.late[,1],
    p.val = peak.DA.stats$promoter.closing.late[,2]) %>%
    mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
    dplyr::filter(!is.na(Padj), beta > 0) %>% 
    dplyr::arrange(Padj) %>% 
    dplyr::mutate(rank = rank(Padj), 
                  sig = ifelse(Padj <0.05, 1, 0), 
                  mlog10Padj = -log10(Padj)) 
}
  
# plot motifs ranked by Pvals 
ggUp_early <- ggplot(peak.DA.stats$promoter.opening.early, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point_rast(size = 1) +
  ggrepel::geom_label_repel(
    #data = peak.DA.stats$promoter.opening[peak.DA.stats$promoter.opening$sig ==1, ],
    data = peak.DA.stats$promoter.opening.early[1:5, ],
    aes(x = rank, y = mlog10Padj, label = motif), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + 
  theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp_late <- ggplot(peak.DA.stats$promoter.opening.late, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point_rast(size = 1) +
  ggrepel::geom_label_repel(
    #data = peak.DA.stats$promoter.opening[peak.DA.stats$promoter.opening$sig ==1, ],
    data = peak.DA.stats$promoter.opening.late[1:5, ],
    aes(x = rank, y = mlog10Padj, label = motif), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + 
  theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
  
ggDo_early <- ggplot(peak.DA.stats$promoter.closing.early, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point_rast(size = 1) +
  ggrepel::geom_label_repel(
    #data = peak.DA.stats$promoter.closing[peak.DA.stats$promoter.closing$sig ==1, ],
    data = peak.DA.stats$promoter.closing.early[1:5, ],
    aes(x = rank, y = mlog10Padj, label = motif), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + 
  theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
ggDo_late <- ggplot(peak.DA.stats$promoter.closing.late, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point_rast(size = 1) +
  ggrepel::geom_label_repel(
    #data = peak.DA.stats$promoter.closing[peak.DA.stats$promoter.closing$sig ==1, ],
    data = peak.DA.stats$promoter.closing.late[1:5, ],
    aes(x = rank, y = mlog10Padj, label = motif), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + 
  theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
ggsave(plot = ggUp_early, 
       filename =  "DAmotifs_PT_OpeningEarlySites_rast_Padj.pdf",
       width = 3.75, 
       height = 3.75)

ggsave(plot = ggUp_late, 
       filename = "DAmotifs_PT_OpeningLateSites_rast_Padj.pdf",
       width = 3.75, 
       height = 3.75)
  
ggsave(plot = ggDo_early, 
       filename = "DAmotifs_PT_ClosingEarlySites_rast_Padj.pdf",
       width = 3.75, 
       height = 3.75)

ggsave(plot = ggDo_late, 
       filename = "DAmotifs_PT_ClosingLateSites_rast_Padj.pdf",
       width = 3.75, 
       height = 3.75)


#####################
# quick check as to the fraction of promoter and distal sites in early vs late categories
data.frame(getPeakSet(prj)) %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end), 
         PT = ifelse(distToTSS < 500, "Promoter", "Distal")) %>% 
  left_join(peak_coefs, by = "peak") %>% 
  filter(p_value < 0.05, direction == "Opening") %>% 
  group_by(direction, PT) %>%
  summarise(n())



