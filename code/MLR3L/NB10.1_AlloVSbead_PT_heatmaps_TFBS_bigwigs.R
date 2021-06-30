basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB10_1"))
setwd(paste0(out_dir, "results/NB10_1"))

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
})

options(dplyr.summarise.inform = FALSE) 

set.seed(2017)

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

# load processed cds (with pseudotime info)
cds = readRDS(file = paste0(out_dir,"results/NB7/cds_NaiveActive_Tcells_PT"))

#broadly define allo and bead stimulations
cdat = colData(cds) %>% 
  data.frame() %>% 
  mutate(stim_type = 
           ifelse(Responder == "stimAlone", "stimAlone", 
                  ifelse(Stimulator == "Bead", "bead", 
                         ifelse(Stimulator == "noStim", "noStim",
                                ifelse(stringr::str_sub(Stimulator, -1) == stringr::str_sub(Responder, -1), "auto", "allo")))))

colData(cds)$stim_type = cdat$stim_type


# count Activated Tcells from each broad stim_type
# including no stim etc. 

T_act_counts = filter(cdat, cellType_broad == "Activated_Tcell") %>% 
  group_by(stim_type) %>% 
  summarise(Tcells_act = n()) 

TcellCounts = group_by(cdat, stim_type) %>% 
  summarise(Tcells_total = n()) %>% 
  left_join(T_act_counts, by = "stim_type")

# Isolate only the bead or allo-Activated Tcells 
idx_cells = 
  (colData(cds)$cellType_broad == "Activated_Tcell" & 
     (colData(cds)$stim_type == "allo" | colData(cds)$stim_type == "bead"))

cds_a = cds[,idx_cells]


####################################################
####################################################
# fit models for peaks over pseudotime across Activated Tcells
# Only run to the next double hash once (ouput saved)! 
##########################
# Create a aggregate cds (columns = number of bins)"

# Use CDS (defined above) only containing activated T-cells
cutoff = ncol(cds_a)*.001
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
saveRDS(binned_cds, file = "cds_actTcell_pseudotime_agg")

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
saveRDS(site_fits, file  = "fit_models_actTcell_pseudotime_agg")
saveRDS(pseudotime_terms, file = "spline_coefs_actTcell_pseudotime_agg")

###################################################
###################################################

# plot raw and smoothed heatmaps of binned data across pseudotime for significant sites 
# load binned cds
cds_b = readRDS(file = "cds_actTcell_pseudotime_agg")
# load models fit over pseudotime (to binned data)
gene_fits_b = readRDS(file = "fit_models_actTcell_pseudotime_agg")

# load regression coefficients 
# note These are the DA sites called (from above) using individual cells over Pseudotime
peak_coefs = read.table(paste0(out_dir,"results/NB10/Tcell_Pseudotime_DApeak_coefs.txt"), head = TRUE) %>% 
  filter(term == "Pseudotime")

DA_sites = filter(peak_coefs, grepl("Pseudotime", x = term), p_value < 0.05) %>% 
  distinct(peak, .keep_all = TRUE)

cds_b_sig <- cds_b[DA_sites$peak,]

# create a set of pseudotime values to predict peak accessibility from 
bin_num = 100 # set the number of bins for smoothing
pseudoorder <- seq(min(colData(cds_b_sig)$mean_pseudo),
                   max(colData(cds_b_sig)$mean_pseudo),
                   length.out=bin_num) 


# function retrieves the fraction of allo to bead stimulated cells in each 
get_frac_allo = function(cds = cds_a, Pseudoorder = pseudoorder){
  cdat = data.frame(colData(cds))
  res = c(NA)
  for (bin in 1:(length(Pseudoorder)-1)){
    bin_low = Pseudoorder[bin]
    bin_high = Pseudoorder[bin + 1]
    cd = filter(cdat, Pseudotime >= bin_low & Pseudotime < bin_high) %>% 
      group_by(stim_type) %>% 
      summarise(count = n())
    if(nrow(cd) > 0){
      allo_cells = cd[cd$stim_type == "allo",2]
      allo_cells = ifelse(nrow(allo_cells) > 0, allo_cells$count, 0)
      bead_cells = cd[cd$stim_type == "bead",2]
      bead_cells = ifelse(nrow(bead_cells) > 0, bead_cells$count, 0)
      frac_allo = allo_cells/(bead_cells+allo_cells)}
    else{
      frac_allo = NA
    }
    res = c(res, frac_allo)
  }
  res
}

bin_frac_allo = get_frac_allo(cds = cds_a, Pseudoorder = pseudoorder)

#get smoothed curves for each peak over pseudotime (calculated using output from previously determined model fits)")
gene_fits_b_sig = dplyr::filter(gene_fits_b, site_name %in% row.names(rowData(cds_b_sig)))

# curves are generated from models for each peak, imputing supplied variables from "new_data" (needs the same names as model input)
curves = model_predictions(gene_fits_b_sig, new_data = data.frame(mean_pseudo=pseudoorder, mean_num_genes = 10000))
rownames(curves) <- gene_fits_b_sig$site_name
#rm(gene_fits_b)
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
bin_df = data.frame("frac_allo" = bin_frac_allo)
rownames(bin_df) = colnames(all)
# set colors for fraction allo labels
ann_colors = list(
  frac_allo = c("white", "firebrick"))

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
  annotation_col = bin_df,
  annotation_colors = ann_colors,
  filename= "smoothed_PTactTcells_Fit_heatmap_2.jpg")

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
  filename= "unsmoothed_PTactTcells_heatmap.jpg")

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
prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE)  

peak_dat = getPeakSet(prj) %>% 
  data.frame() %>% 
  mutate(site_name = paste(seqnames, start, end, sep = "_")) %>% 
  dplyr::select(site_name, distToTSS, nearestGene, nearestTSS, GC)

DA.stats = peak_coefs %>% 
  mutate(site_name = peak) %>% 
  left_join(peak_dat, by = "site_name")

write.table(DA.stats, file = "Tcell_Pseudotime_DApeak_coefs_directions.txt", 
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
site_motif_stats = find_sig_motifs(site_coef_directions = DA.stats, archr_project = prj)
plot_top_motifs(peak.DA.stats = site_motif_stats, Feat_type = "promoter")
plot_top_motifs(peak.DA.stats = site_motif_stats, Feat_type = "distal")



#############################
# Make Bigwig tracks
#############################

# list of Activated T-cells 
ActTcells = colData(cds) %>% 
  data.frame() %>% 
  mutate(stim_type = 
           ifelse(Responder == "stimAlone", "stimAlone", 
                  ifelse(Stimulator == "Bead", "bead", 
                         ifelse(Stimulator == "noStim", "noStim",
                                ifelse(stringr::str_sub(Stimulator, -1) == stringr::str_sub(Responder, -1), "auto", "allo"))))) %>% 
  filter(cellType_broad == "Activated_Tcell", 
         stim_type %in% c("allo", "bead"))

# re-aggregate cells based on pseudotime bins (for bulk bw tracks)
set.seed(2017)
ActTcells$Pseudotime_bin <- kmeans(ActTcells$Pseudotime, centers=5)$cluster 

bin_ranks = group_by(ActTcells, Pseudotime_bin) %>% 
  summarise(mean_pseudo = mean(Pseudotime)) %>% 
  mutate(PT_bin_rank =  rank(mean_pseudo, ties.method = 'first'), 
         PT_bin_rank_char = paste0("PT_bin_", PT_bin_rank)) %>% 
  dplyr::select(Pseudotime_bin, PT_bin_rank_char)

ActTcells = left_join(ActTcells, bin_ranks, by = "Pseudotime_bin")
  
table(ActTcells$PT_bin_rank_char)
# generate bigwig tracks for activated T-cells only 
prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 

prj = prj[prj$cellNames %in% ActTcells$cell,]

#Add info project metadata (temporarily for making bigwig tracks)
prj <- addCellColData(ArchRProj = prj, data = ActTcells$PT_bin_rank_char,
                      cells = ActTcells$cell, name = "PT_bin_rank_char", 
                      force = TRUE)

prj <- addCellColData(ArchRProj = prj, data = ActTcells$stim_type,
                      cells = ActTcells$cell, name = "ActTcell_stim_type", 
                      force = TRUE)


# Make bigwig tracks

# by pseudotime bin
getGroupBW(
  ArchRProj = prj,
  groupBy = "PT_bin_rank_char",
  normMethod = "ReadsInTSS", # creates a group scale factor = 10k/sum(reads in TSS) 
  #normMethod = "nFrags", # creates a group scale factor = 10k/sum(nFrags)
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

# by Allo vs Bead 
getGroupBW(
  ArchRProj = prj,
  groupBy = "ActTcell_stim_type",
  normMethod = "ReadsInTSS", # creates a group scale factor = 10k/sum(reads in TSS) 
  #normMethod = "nFrags", # creates a group scale factor = 10k/sum(nFrags)
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

###########################
# plot browser track
###########################

###
# get sig DA genes (just for browsing)
coefs_PT = read.table(paste0(out_dir, "results/NB13/Tcell_Pseudotime_DAgene_coefs.txt"), head = TRUE) %>% 
  filter(term == "Pseudotime", 
         p_value < 0.05) %>% 
  arrange(p_value)
###

# Use GVIZ for plotting pseudobulk data from bigwig

plot_pb_browsertracks = function(archrPrj = prj, bw_list = actTcell_bw_list,
                                 chrom = "chr5", start_window = 95284800, 
                                 end_window = 95316600, locus = "locus", Ylim = c(0,20))
{
  dir.create("browser_tracks/")
  ## Load BW files to data tracks
  PT_1 <- DataTrack(bw_list$PT_1, type = "histogram", fill = "#1D1147FF",
                            col.histogram = "#1D1147FF", name = "PT_1", 
                            ylim = Ylim)
  PT_2 <- DataTrack(bw_list$PT_2, type = "histogram", fill = "#822681FF", 
                              col.histogram = "#822681FF", name = "PT_2", 
                              ylim = Ylim)
  PT_3 <- DataTrack(bw_list$PT_3, type = "histogram", fill = "#B63679FF",
                              col.histogram = "#B63679FF", name = "PT_3", 
                              ylim = Ylim)
  PT_4 <- DataTrack(bw_list$PT_4, type = "histogram", fill = "#FB8861FF",
                            col.histogram = "#FB8861FF", name = "PT_4", 
                            ylim = Ylim)
  PT_5 <- DataTrack(bw_list$PT_5, type = "histogram", fill = "#FEC287FF", 
                            col.histogram = "#FEC287FF", name = "PT_5", 
                            ylim = Ylim)
  bead <- DataTrack(bw_list$bead, type = "histogram", fill = "darkblue", 
                             col.histogram = "darkblue", name = "Bead", 
                             ylim = Ylim)
  allo <- DataTrack(bw_list$allo, type = "histogram", fill = "firebrick", 
                             col.histogram = "firebrick", name = "Allo", 
                             ylim = Ylim)
  
  # prepare extra feature tracks to display
  bm <- useMart(host = "grch37.ensembl.org", 
                biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl")
  
  bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chrom,
                                start = start_window, end = end_window, 
                                filter = list(with_ox_refseq_mrna = TRUE),
                                stacking = "dense")
  
  ideoTrack <- IdeogramTrack(genome="hg19", chromosome= chrom)
  axisTrack <- GenomeAxisTrack()
  
  allpeaks = getPeakSet(archrPrj)
  # must restrict GRanges object to specific chromosome
  peaks = allpeaks[allpeaks@seqnames == chrom,]
  pTrack <- AnnotationTrack(peaks, genome= "hg19", name = "peaks", 
                            col = "darkblue", fill = "darkblue")
  
  #plot all tracks together
  pdf(paste0("browser_tracks/gviz_", locus, ".pdf"), width = 3, height = 6)
  print(
    plotTracks(c(ideoTrack, 
                 axisTrack, 
                 PT_1,
                 PT_2, 
                 PT_3, 
                 PT_4, 
                 PT_5, 
                 bead, 
                 allo, 
                 pTrack,
                 bmt),
               sizes = c(0.5, 1,1,1,1,1,1,1,1, 0.5, 0.5), 
               from = start_window,
               to = end_window, 
               chromosome = chrom, 
               type = "histogram")
  )
  dev.off()
}

# SAHA tracks
actTcell_bw_list = list(
  PT_1 = paste0(out_dir,"mlr_filtered_annotated/GroupBigWigs/PT_bin_rank_char/PT_bin_1-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  PT_2 = paste0(out_dir,"mlr_filtered_annotated/GroupBigWigs/PT_bin_rank_char/PT_bin_2-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  PT_3 = paste0(out_dir,"mlr_filtered_annotated/GroupBigWigs/PT_bin_rank_char/PT_bin_3-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  PT_4 = paste0(out_dir,"mlr_filtered_annotated/GroupBigWigs/PT_bin_rank_char/PT_bin_4-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  PT_5 = paste0(out_dir,"mlr_filtered_annotated/GroupBigWigs/PT_bin_rank_char/PT_bin_5-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  bead = paste0(out_dir,"mlr_filtered_annotated/GroupBigWigs/ActTcell_stim_type/bead-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  allo = paste0(out_dir,"mlr_filtered_annotated/GroupBigWigs/ActTcell_stim_type/allo-TileSize-100-normMethod-ReadsInTSS-ArchR.bw")
)


# browser track (ELL2)
plot_pb_browsertracks(archrPrj = prj, bw_list = actTcell_bw_list,
                      chrom = "chr5", start_window = 95283800, 
                      end_window = 95317600, locus = "ELL2", Ylim = c(0,30))





###########################

# Not finished (probably not necessary)

#Define early/late changes as sites where 50% max reached
#before/after middle pseudotime (bin 50)
early_up = up_cp[up_cp < 50] 
late_up = up_cp[up_cp >= 50] 
early_down = down_cp[down_cp < 50] 
late_down = down_cp[down_cp >= 50] 

peak_coefs$timing <- sapply(peak_coefs$peak, function(x) {
  if(x %in% names(early_up)) return("early")
  if(x %in% names(late_up)) return("late")
  if(x %in% names(early_down)) return("early")
  if(x %in% names(late_down)) return("late")
  if(x %in% names(trans_cp)) return("Transient")
  return("Static")
})













