basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190521_scichem2_AllPlates/"
bin_directory = paste0(basepath, "analysis/bin/")
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB8"))
setwd(paste0(out_dir, "results/NB8/"))

PD_directory = paste0(out_dir, "results/NB7/")
DA_directory = paste0(out_dir, "results/NB5/")

suppressPackageStartupMessages({
  library(ArchR)
  library(ggridges)
  library(tidymodels)
  library(devtools)
  library(monocle3)
  library(furrr)
  library(viridis)
  library(piano)
  library(UpSetR)
  library(snowfall)
  library(plyr)
  library(dplyr)
  library(cicero)
  plan(multicore)
})

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)
set.seed(1989) 

prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered"),
                       showLogo = FALSE)

# function prints raw and smoothed heatmaps of binned data across pseudodose for significant sites 
## also returns the input coeficient table with 
plot_pd_heatmap = function(d = "SAHA"){
  cat("loading ", d, " data for making heatmaps\n")
  dir.create(d)
  cds_b = readRDS(paste0(PD_directory, d, "/ATAC_", d, "_agg_pseudodose_cds"))
  rowData(cds_b)$site_name = row.names(rowData(cds_b))
  gene_fits_b = readRDS(paste0(PD_directory, d, "/ATAC_", d, "_binned_spline_models_by_site"))
  gene_fits_b$site_name = gene_fits_b$peak
  
  # Note if you want to use sites called via different methods, change the file below
  # terms have altered nomenclature for drugs
  drug_simple = function(drug){
    res = ifelse(drug == "BMS345541", "BMS",
            ifelse(drug == "Nutlin3A", "Nutlin", drug))
    return(res)
  }
  term_drug = drug_simple(d)
  dose_coefs_all = read.table(paste0(DA_directory, "sciPlex_ATAC_A549_small_screen_DApeaks.txt"), head = TRUE) %>% 
    dplyr::filter(term == paste0("dose_", term_drug)) %>% 
    mutate(site_name = as.character(peak))
  
  DA_sites = dose_coefs_all %>% 
    dplyr::filter(q_value < 0.05) 
  
  cds_b_sig <- cds_b[DA_sites$site_name,]
  pseudoorder <- seq(min(colData(cds_b_sig)$mean_pseudo), max(colData(cds_b_sig)$mean_pseudo), length.out=100)
  
  cat("get smoothed curves for each peak over pseudotime (calculated using output from preiously determined model fits)\n")
  gene_fits_b_sig = dplyr::filter(gene_fits_b, site_name %in% row.names(rowData(cds_b_sig)))
  # curves are generated from models for each peak, imputing supplied variables from "new_data" (needs the same names as model input)
  curves = model_predictions(gene_fits_b_sig, new_data = data.frame(mean_pseudo=pseudoorder, mean_num_genes = 10000))
  rownames(curves) <- gene_fits_b_sig$site_name
  rm(gene_fits_b)
  rm(gene_fits_b_sig)
  
  cat("splitting sites into transient, opening, and closing...\n") 
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
  
  cat("ordering sites by pseudodose at which half the max accessibility is reached...\n")
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
  
  
  cat("printing smoothed heatmap to file: ", d, "_smoothed_pdFit_heatmap.jpg...\n") 
  options(bitmapType='cairo')
  
  pheatmap::pheatmap(all,
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
                     filename= paste0(d, "/",  d, "_smoothed_pdFit_heatmap.jpg"))
  
  cat("printing unsmoothed heatmap to file: ",  d, "_smoothed_pdFit_heatmap.jpg...\n")
  
  cell_order <- colData(cds_b_sig)[order(colData(cds_b_sig)$mean_pseudo),]$cell_bin
  un_smoothed <- counts(cds_b_sig)[row.names(all), cell_order]
  un_smoothed <- t(t(un_smoothed)/colData(cds_b_sig)[cell_order,]$ncells)
  
  pheatmap::pheatmap(un_smoothed,
                     #annotation_colors = list(Type = c(Proximal = "#241023", `Other Distal` = "#47A025", `Distal Enhancer` = "#735CDD")),
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
                     filename= paste0(d, "/",  d, "_heatmap_unsmoothed.jpg"))
  
  # Add column to coeficient table with direction of DA changes
  dose_coefs_all$direction <- sapply(dose_coefs_all$site_name, function(x) {
    if(x %in% names(up_cp)) return("Opening")
    if(x %in% names(down_cp)) return("Closing")
    if(x %in% names(trans_cp)) return("Transient")
    return("Static")
  })
  return(dose_coefs_all)
}

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

#############################################
# Isolate the top 6 ranked motifs for each catagory from distal or promoter sites
plot_top_motifs = function(peak.DA.stats, Feat_type = "distal", d = "SAHA"){
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
  
  pdf(paste0(d, "/", d, "_", Feat_type,".peak.motif.enrichments_beta.pdf"), width = 2, height=2.4)
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


################################################################################
# run the functions above for each drug. 

#for(D in c("SAHA", "Dex", "BMS345541", "Nutlin3A")){
for(D in c("Nutlin3A")){
  DA.stats = plot_pd_heatmap(d = D)
  write.table(DA.stats, file = paste0(D, "/", D, "_", "peak_acc_coefs_PDdirections.txt"), 
              quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  site_motif_stats = find_sig_motifs(site_coef_directions = DA.stats, archr_project = prj)
  plot_top_motifs(peak.DA.stats = site_motif_stats, Feat_type = "promoter", d = D)
  plot_top_motifs(peak.DA.stats = site_motif_stats, Feat_type = "distal", d = D)
}


##########################################################################




