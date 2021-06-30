# This script was developed to assess the predictive power of sequence features (motifs),
# both within the gene's promoter and within cicero-connected distal regulatory elements,
# on gene expression. To achieve this, this script does the following
# 1) identify cicero connection on the binary peak by cell count matrix (for all conditions).
# 2) retreive motif sequences present within each peak (based on ArchR database)
# 3) for each gene, identify the peak which most closely represents the promoter.

## Subsequent steps and modeling are performed separately for each drug in screen

# 4) for each promoter peak, retrieve and summarize all cicero-connected peaks and drug 
# effects on RNA expression. 
# 5) build a promter peak x motif matrix where each motif is represented twice, once for whether
# it is present in the promoter, and again for it's presence within connected peaks.
# Note that a motif presence within cicero-connected sites is wieghted (mulitplied) by 
# each connected peak's "connection score".
# 6) using pre-existing Expression data for each gene from the same conditions, we then
# run linear models to predict a given genes expression (or drug-dependent expression 
# coefficient) based on each motif presence solely within the promoter or within the
# promoter and connected distal sites. 
# 7) ultimately, we hope these models reveal sequence elements within promoters and connected
# distal sites which predict how genes will react to drugs. 

basepath = "github/"
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB10"))
setwd(paste0(out_dir, "results/NB10/"))


suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(cicero)
  library(ggplot2)
  library(ggridges)
  library(tidymodels)
  library(ggrastr)
  library(ggpubr)
  library(glmnet)
})

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)
set.seed(1989) 

# load processed ArchR project
prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered"),
                       showLogo = FALSE)

##############################################
# 1) Get Cicero connections
##############################################
# only need to run once (output saved)
# make a peak x cell cds object for running cicero. 

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

cds_a = monocle3::new_cell_data_set(assays(pMat)$PeakMatrix, 
                                       cell_metadata = colData(pMat), 
                                       gene_metadata = a_rd)

#find connections from ALL cells 
# preprocess cds for cicero: 
message("Getting Cicero connections using All cells\n")

cds_archr = cds_a

set.seed(2017)
cds_archr <- detect_genes(cds_archr)
cds_archr = cds_pl[rowData(cds_archr)$num_cells_expressed > 0,]
cds_archr <- estimate_size_factors(cds_archr)
cds_archr = preprocess_cds(cds_archr, method = "LSI", num_dimensions=50)
cds_archr = reduce_dimension(cds_archr, reduction_method = 'UMAP', preprocess_method = "LSI")

umap_coords <- reducedDims(cds_archr)$UMAP
cicero_cds <- make_cicero_cds(cds_archr, reduced_coordinates = umap_coords)

# run cicero on preprocessed data 
data("human.hg19.genome")
conns <- run_cicero(cicero_cds, 
                    human.hg19.genome, 
                    window = 5e+05,
                    silent = FALSE,
                    sample_num = 100)

filter(conns, coaccess > 0) %>% 
  write.table(file = "cicero_peak_conns.txt", 
              sep = "\t",
              quote = FALSE,
              row.names = FALSE, 
              col.names = TRUE)

#############

## recalculate cicero connections from subsetted drug treatments
drugs = c("Dex", "SAHA", "BMS", "Nutlin")
for (d in drugs){
  
  dir.create(d)
  message("Getting Cicero connections using", d, "-treated cells only\n")
  
  # restrict cds to cells from each treatment
  cds_archr = cds_a[,colData(cds_a)$treatment == d]
  
  # preprocess cds for cicero: 
  set.seed(2017)
  cds_archr <- detect_genes(cds_archr)
  #cds_archr = cds_archr[rowData(cds_archr)$num_cells_expressed > 0,]
  cds_archr <- estimate_size_factors(cds_archr)
  cds_archr = preprocess_cds(cds_archr, method = "LSI", num_dimensions=50)
  cds_archr = reduce_dimension(cds_archr, reduction_method = 'UMAP', preprocess_method = "LSI")

  umap_coords <- reducedDims(cds_archr)$UMAP
  cicero_cds <- make_cicero_cds(cds_archr, reduced_coordinates = umap_coords)

  # run cicero on preprocessed data 
  data("human.hg19.genome")
  conns <- run_cicero(cicero_cds, 
                      human.hg19.genome, 
                      window = 5e+05,
                      silent = FALSE,
                      sample_num = 100)

  filter(conns, coaccess > 0) %>% 
    write.table(file = paste0(d, "/", d, "_cicero_peak_conns.txt"), 
                sep = "\t",
                quote = FALSE,
                row.names = FALSE, 
                col.names = TRUE)
}
##############################################
# 2) retrieve motif sequences present within each accessible peak 
# (based on ArchR database)
##############################################

# load peak by motif matrix: 
pull_peak_by_motif_mat = function(archr_project){
  motif_dat = readRDS(paste0(getOutputDirectory(archr_project), "/Annotations/Motif-In-Peaks-Summary.rds"))
  motif_mat = assay(motif_dat$motifMatches)
  peaks = data.frame(rowRanges(motif_dat$motifMatches)) %>% # GB Corrected 11/17/2020
    dplyr::mutate(peak = paste(seqnames, start, end, sep = "_"))
  row.names(motif_mat) <- peaks$peak
  return(motif_mat)
}

motif_mat = pull_peak_by_motif_mat(prj)

# reformat motif_mat from logical to numeric
mmf = matrix(as.numeric(motif_mat), nrow = nrow(x=motif_mat))
row.names(mmf) = row.names(motif_mat)
colnames(mmf) = colnames(motif_mat)
motif_df = as.data.frame(mmf)
motif_df$site_name = row.names(motif_df)

##############################################
# 3) for each gene, identify the peak which most closely represents the promoter.
##############################################

# Identify Motifs Present within Promoters and 
# Motifs present in linked distal sites for each promoter
# load and format various peak and gene information

DA.res = read.table(paste0(out_dir, "results/NB5/sciPlex_ATAC_A549_small_screen_DApeaks.txt"), head = T)
RNA.res = read.table(file = paste0(out_dir, "results/NB5/sciPlex_RNA_A549_small_screen_DEGs.txt"), head = T)
motif.df = motif_df

cons.info = 
  read.table("cicero_peak_conns.txt", head = TRUE) %>% 
  mutate(p1_center = (as.numeric(stringr::str_split_fixed(Peak1, pattern = "_", 3)[,2]) + 
                        as.numeric(stringr::str_split_fixed(Peak1, pattern = "_", 3)[,3]))/2, 
         p2_center = (as.numeric(stringr::str_split_fixed(Peak2, pattern = "_", 3)[,2]) + 
                        as.numeric(stringr::str_split_fixed(Peak2, pattern = "_", 3)[,3]))/2, 
         dist = abs(p2_center - p1_center))


######################################################
# 4) for each promoter peak, retrieve and summarize all cicero-connected peaks and drug 
# effects on RNA expression. 
######################################################

# function returns matrix where each row is a promoter site 
# columns contain summaries of cicero connected sites, 
# drug dependent effects on gene expression 
# and weighted values for each motifs presence within 
# promoters and distal sites
get_TSS_info = function(d, DA_res = DA.res, motif_df = motif.df, 
                        cons_info = cons.info, RNA_res = RNA.res, 
                        cicero_threshold = 0.25){
  dir.create(d)
  
  # format peaks into promoters and distal sites. Keep info about drug dependent changes
  da_res = 
    filter(DA_res, term == paste0("dose_", d)) %>% 
    rename(site_name = peak) %>% 
    mutate(gene_short_name = ifelse(peakType == "Promoter" & distToTSS < 500, 
                                    nearestGene, 'NA'), 
           peak.type = ifelse(peakType == "Promoter" & distToTSS < 500, 
                              "Promoter", 'Distal'), 
           promoter_gene = gene_short_name, 
          direction = ifelse(q_value < 0.05 & estimate > 0, "Opening", 
                              ifelse(q_value < 0.05 & estimate < 0, "Closing", "Static"))) %>% 
    left_join(RNA_res, by = "gene_short_name") %>% 
    select( c(nearestGene, peak.type, distToTSS, 
              nearestTSS, GC, site_name, gene_short_name, 
              promoter_gene, id, direction)) %>% 
    rename(tss.id = id) %>% 
    distinct(site_name, .keep_all = TRUE)

  # select and rename only the necessary columns
  peak.info = arrange(da_res, site_name) %>% 
    dplyr::select(peak = site_name, peak.type = peak.type,
                  promoter.gene = promoter_gene, tss.id = tss.id, acc.trend = direction)

  row.names(peak.info) = NULL

  # load the cicero connection info
  glasso.edges = cons_info[, c("Peak1", "Peak2", "dist", "coaccess")]
  
  # set thresholds and define cicero-connections (graphical lasso)
  glasso.edges = glasso.edges[!is.na(glasso.edges$dist) & glasso.edges$dist >= 1000 & 
                                !is.na(glasso.edges$coaccess) & glasso.edges$coaccess > cicero_threshold,]

  # not sure what this part does (it creates a very redundant matrix with each promoter represented many times)
  promoter.info = do.call(rbind, lapply((peak.info %>% filter(peak.type == "Promoter"))$tss.id,
                                        function(x) {
                                          all.tss = strsplit(x, " ")[[1]]
                                          data.frame(tss.id.melt = all.tss, tss.id = x, stringsAsFactors = F)
                                        })) %>% inner_join(peak.info %>% filter(peak.type == "Promoter"), by = "tss.id")
  #save(promoter.info, file=paste0(d, "/", d, "_promoter.info"))

  # get list of motif names
  motifs = names(motif_df)[1:(ncol(motif_df)-1)]

  # add motif presence (0 or 1) to the promoter.info peak list  
  promoter.elements = data.frame(left_join(promoter.info %>%
                                             dplyr::select(peak, tss.id = tss.id.melt, p.acc.trend = acc.trend), 
                                           motif_df, by = c("peak" = "site_name")))

  # add "p." prefix to motifs (for presence in promoter peak)
  for (motif in motifs) {
    promoter.elements[, paste("p.", motif, sep="")] =
      ifelse(is.na(promoter.elements[, motif]), 0, promoter.elements[, motif])
    promoter.elements[, motif] = NULL
  }
  
  # reduces the redundancy to one row per promoter peak
  # I believe this also looks a promoters with multiple peaks and pools their motifs
  tss.promoter.motifs = reshape2::dcast(
    reshape2::melt(promoter.elements, id.vars = c("peak", "tss.id", "p.acc.trend"),
         variable.name = "feature", value.name = "value"),
    tss.id ~ feature, value.var = "value", fill = 0, fun.aggregate = max)

  # create binary motif matrix for peaks defined as "distal"
  distal.elements = filter(peak.info, peak.type == "Distal") %>%
    dplyr::select(peak, acc.trend) %>%
    left_join(motif_df, by = c("peak" = "site_name"))
  
  # add "d." prefix to motifs (for presence/absence in distal peaks)
  for (motif in motifs) {
    distal.elements[, paste("d.", motif, sep="")] =
      ifelse(is.na(distal.elements[, motif]), 0, distal.elements[, motif])
    distal.elements[, motif] = NULL
  }

  # pair distal sites with their connected promoters
  # Note that gene promoter with multiple promoter peaks can be connected to the 
  # same distal peak through each promoter peak
  # identify maximum co-accessibility score for each TSS - distal peak connection. 
  tss.distal.peaks = data.frame(promoter.elements[, c("tss.id", "peak")] %>%
                                  inner_join(glasso.edges, by = c("peak" = "Peak1")) %>%
                                  group_by(tss.id, Peak2) %>% 
                                  summarize(coaccess = max(coaccess)) %>%
                                  inner_join(distal.elements, by = c("Peak2" = "peak")))

  # multiply motif presense (1 or 0) by maximum co-accessibility of any connected peak containing the motif
  for (motif in motifs) {
    tss.distal.peaks[, paste("d.", motif, sep="")] =
      tss.distal.peaks[, "coaccess"] * tss.distal.peaks[, paste("d.", motif, sep="")]
  }

  # for each promoter, summarize stats for distal peak connections 
  tss.distal.peak.stats = group_by(tss.distal.peaks[, 1:4], tss.id) %>%
    summarize(
      n.dp = n(),
      n.o.dp = sum(acc.trend == "Opening"),
      n.c.dp = sum(acc.trend == "Closing"),
      best.o.dp = max(ifelse(acc.trend == "Opening", coaccess, 0)),
      best.c.dp = max(ifelse(acc.trend == "Closing", coaccess, 0)))

  # These next two steps take the sum of motif values from all connected peaks to each TSS
  # melt the tss by distal motif matrix (long form)
  tss.distal.motifs = reshape2::melt(tss.distal.peaks[, c(1, 5:ncol(tss.distal.peaks))],
                           id.vars = c("tss.id"), variable.name = "feature", value.name = "value")

  # dcast into matrix with one row for each TSS, with motif scores now equal to the 
  # sum of scores for each distally connected peak to that TSS
  tss.distal.motifs = dcast(tss.distal.motifs, tss.id ~ feature, value.var = "value", fill = 0,
                            fun.aggregate = function(x) sum(x))
  
  message("number of promoters with connected distal sites = ", dim(tss.distal.motifs)[1], "\n")
  
  # now that we have one distal motif score for each TSS.id, join that matrix to the promoter motif scores 
  tss.motifs = data.frame(right_join(tss.distal.peak.stats, tss.promoter.motifs, by = "tss.id")) %>%
    left_join(tss.distal.motifs, by = "tss.id")

  # change any NAs to 0
  for (i in 2:ncol(tss.motifs)) {
    tss.motifs[, i] = ifelse(is.na(tss.motifs[, i]), 0, tss.motifs[, i])
  }
  
  message("total number of promoters in motif analysis = ", dim(tss.motifs)[1], "\n")
  
  return(tss.motifs)
}


##############################################
# 6 and 7) Create glm to predict gene expression coefficient from 
# promoter motifs and promoter + linked distal motifs
#############################################

# function to return cv.glmnet model res for sequence motif predictive effect 
# on drug dependent changes in gene expression 
# (***only predicting for DEGs***)

run_motif_glm_model = function(d, motif_df, tss.motifs, gene_model_res, 
                               DEG_only = TRUE, has_conn_only = FALSE) {
  # Add promoter and distal motif matrix to Gene expression coefficients
  motifs = names(motif_df)[1:(ncol(motif_df)-1)]
  tss.RNA.df = 
    gene_model_res %>% 
    dplyr::filter(term == paste0("dose_", d)) %>% 
    dplyr::select(tss.id = id, gene = gene_short_name, term, expr.estimate = estimate, pval = p_value, 
                  qval = q_value, num.cells.expr = num_cells_expressed) %>% 
    dplyr::inner_join(tss.motifs, by = "tss.id")
  
  # **** Restrict analysis to DEGs *****
  if(DEG_only){
  tss.RNA.df = filter(tss.RNA.df, qval < 0.05)}
  
  # **** Restrict analysis to DEG promoters with > 0 distal connections ******
  if(has_conn_only){
  tss.RNA.df = filter(tss.RNA.df, n.dp > 0)}
  
  message("final number of genes modeled = ", dim(tss.RNA.df)[1], "\n")

  models = list()

  tmp = filter(tss.RNA.df, !is.na(expr.estimate))

  # GLM for promoter motifs vs expression coefficients
  models[[paste("p", "Expression", sep=".")]] = cv.glmnet(
    as.matrix(tmp[, paste("p.", motifs, sep="")]),
    as.vector(tmp$expr.estimate),
    family = "gaussian", alpha = 0.5)

  # GLM for promoter and connected distal motifs vs expression coefficients
  models[[paste("pd", "Expression", sep=".")]] = cv.glmnet(
    as.matrix(tmp[, c(paste("p.", motifs, sep=""), paste("d.", motifs, sep=""))]),
    as.vector(tmp$expr.estimate),
    family = "gaussian", alpha = 0.5)
  
  return(models)
}



# Retrieve model R2 values for promoter model, and promoter and distal model for plotting
make_r2_df = function(models) {
  # quick function to get r2 values
  get.model.R2 = function(model) {
    model$glmnet.fit$dev.ratio[which(model$lambda == model$lambda.min)]
  }

  model.r2.df = data.frame(
    target = "Expression",
    p.r2 = get.model.R2(models[["p.Expression"]]),
    pd.r2 = get.model.R2(models[["pd.Expression"]])) %>%
    mutate(ratio = pd.r2 / p.r2, plot.label = round(ratio, 2))

  rownames(model.r2.df) = NULL
  model.r2.df
  
  # plot change in R2 value for promoters and promoters + distal site motifs 
  plot.df = reshape2::melt(model.r2.df %>% filter(target == "Expression") %>% dplyr::select(target, plot.label, p.r2, pd.r2),
                 id.vars = c("target", "plot.label"), variable.name = "model", value.name = "r2")
  
  return(plot.df)
}


# Function for getting the top motifs and coefficients for plotting
extract_top_motif_coefs = function(models){
  
  model.coef = list()
  
  # extract model coefs when lambda achieves error within 1 se of minimum
  model.coef$Expression = coef(models$pd.Expression,
                               s = "lambda.1se")[-1, 1]

  # reformat motif coefficients (for predicting expression) into df
  coef.mat = do.call(cbind, model.coef)
  coef.df = reshape2::melt(coef.mat)
  names(coef.df) = c("feature", "target", "value")
  coef.df = coef.df %>% 
    group_by(feature) %>%
    mutate(best.value = max(value)) %>% 
    ungroup() %>%
    filter(abs(best.value) >= 0.005) %>%
    inner_join(coef.df %>% 
                 filter(target == "Expression") %>%
                 dplyr::select(feature, Expression.value = value), by = "feature") %>%
    arrange(-Expression.value)

  coef.df$feature = factor(coef.df$feature, levels = rev(unique(coef.df$feature)))
  coef.df$distal.feature = grepl("^d[.]", coef.df$feature)

  name.map  = list(
    "d.MEF2.FAM" = "MEF2 family",
    "d.MEIS1" = "MEIS1",
    "p.MEIS1" = "MEIS1",
    "d.MRF.E.BOX" = "MRF E-box",
    "p.TEAD.FAM" = "TEAD family"
  )

  coef.df$feature.label = sapply(as.character(coef.df$feature), function(x) {
    if (x %in% names(name.map)) {
      name.map[[x]]
    } else { x }
  })
  
 return(coef.df) 
}

###############################################
# run functions for each drug 
# saha 
saha.tss.motifs = get_TSS_info(d = "SAHA", cicero_threshold = 0.1)
save(saha.tss.motifs, file = "SAHA/SAHA_tss.motifs")
saha.models = run_motif_glm_model(d = "SAHA", 
                                 motif_df = motif.df,
                                 tss.motifs = saha.tss.motifs, 
                                 gene_model_res = RNA.res, 
                                 DEG_only = TRUE,
                                 has_conn_only = FALSE)
saha.model.r2.df = make_r2_df(models = saha.models)
saha.model.coefs = extract_top_motif_coefs(models = saha.models)

# dex 
dex.tss.motifs = get_TSS_info(d = "Dex", cicero_threshold = 0.1)
save(dex.tss.motifs, file = "Dex/Dex_tss.motifs")
dex.models = run_motif_glm_model(d = "Dex", 
                                 motif_df = motif.df,
                                 tss.motifs = dex.tss.motifs, 
                                 gene_model_res = RNA.res,
                                 DEG_only = TRUE,
                                 has_conn_only = FALSE)
dex.model.r2.df = make_r2_df(models = dex.models)
dex.model.coefs = extract_top_motif_coefs(models = dex.models)

# BMS
bms.tss.motifs = get_TSS_info(d = "BMS", cicero_threshold = 0.1)
save(bms.tss.motifs, file = "BMS/BMS_tss.motifs")
bms.models = run_motif_glm_model(d = "BMS", 
                                 motif_df = motif.df,
                                 tss.motifs = bms.tss.motifs, 
                                 gene_model_res = RNA.res,
                                 DEG_only = TRUE,
                                 has_conn_only = FALSE)
bms.model.r2.df = make_r2_df(models = bms.models)
bms.model.coefs = extract_top_motif_coefs(models = bms.models)

# Nutlin
nutlin.tss.motifs = get_TSS_info(d = "Nutlin", cicero_threshold = 0.1)
save(nutlin.tss.motifs, file = "Nutlin/Nutlin_tss.motifs")
nutlin.models = run_motif_glm_model(d = "Nutlin", 
                                 motif_df = motif.df,
                                 tss.motifs = nutlin.tss.motifs, 
                                 gene_model_res = RNA.res,
                                 DEG_only = TRUE,
                                 has_conn_only = FALSE)
nutlin.model.r2.df = make_r2_df(models = nutlin.models)
nutlin.model.coefs = extract_top_motif_coefs(models = nutlin.models)


plot_model_res = function(d, plot.df, coef.df){
  pdf(paste0(d,"/p_vs_pd_Motif_expr.model.r2.pdf"), width = 0.8, height = 2.0)
  print(
  ggplot(plot.df, aes(x = model, y = r2)) +
    geom_bar(aes(fill = model), stat = "identity", position = "dodge") +
    geom_text(aes(y = r2 + 0.03, label = plot.label, group = model, alpha = model), size = 2.2) +
    xlab("Model") +
    ylab("Cross-validated R^2") +
    scale_x_discrete(labels = c("Promoter motifs", "Promoter and\ndistal motifs")) +
    scale_y_continuous(labels = percent) +
    scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
    scale_alpha_manual(values = c(0, 1)) +
    guides(fill = F, alpha = F) +
    monocle3:::monocle_theme_opts() +
    theme(text = element_text(size=6)) +
    theme(axis.line = element_line(),
          axis.title.x = element_text(margin = margin(5, 0, 0, 0)),
          axis.title.y = element_text(margin = margin(0, 4, 0, 0)),
          axis.text.x = element_text(angle = 90, vjust = 0.5))
  )
  dev.off()

  pdf(paste0(d,"/Expression.model.Motif.coef.pdf"),  width = 2.5, height = 2.0)
  print(
  ggplot(coef.df, aes(x = feature, y = value, fill = distal.feature)) +
    facet_wrap(~ target) +
    coord_flip() +
    scale_x_discrete(labels = rev(coef.df$feature.label[seq(1, nrow(coef.df), 3)])) +
    #scale_y_continuous(breaks = c(-0.01, 0, 0.01), labels = c(-0.01, 0, 0.01)) +
    scale_fill_manual(values = c("grey70", "firebrick"), labels = c("Promoter", "Distal")) +
    geom_bar(stat = "identity") +
    xlab("Motif family") +
    ylab("Regression coefficient") +
    guides(fill = guide_legend(title = "Element type")) +
    theme_bw(base_size = 6) +
    monocle3:::monocle_theme_opts() +
    theme(axis.line = element_line(),
          axis.title.x = element_text(margin = margin(4, 0, 0, 0)),
          axis.title.y = element_text(margin = margin(0, 4, 0, 0)),
          legend.key.size = unit(0.3, "cm"),
          legend.margin = margin(10, 0, 10, -5))
  )
  dev.off()
}



plot_model_res(d = "SAHA", plot.df = saha.model.r2.df, coef.df = saha.model.coefs)
plot_model_res(d = "Dex", plot.df = dex.model.r2.df, coef.df = dex.model.coefs)
plot_model_res(d = "BMS", plot.df = bms.model.r2.df , coef.df = bms.model.coefs)
plot_model_res(d = "Nutlin", plot.df = nutlin.model.r2.df, coef.df = nutlin.model.coefs)



##########################################################################
##########################################################################
# quick check cicero conns. 
cic_conns = read.table("cicero_peak_conns.txt", head = TRUE)

# download and unzip gtf 
temp <- tempfile()
download.file("ftp://ftp.ensembl.org/pub/release-65/gtf/homo_sapiens/Homo_sapiens.GRCh37.65.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name


# look at conns around  COPS7A locus (SAHA-affected)
pdf("test_COPS7A_conns_browser.pdf", width = 5, height=3)
print(
plot_connections(cic_conns, "chr12", 6800304, 6865007,
                 gene_model = gene_anno, 
                 coaccess_cutoff = .1, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )
)
dev.off()


# look at conns around MT1X (Dex-affected)
pdf("MT1X_conns_browser.pdf", width = 5, height=3)
print(
  plot_connections(cic_conns, "chr16", 56697397, 56739939,
                   gene_model = gene_anno, 
                   coaccess_cutoff = .1, 
                   connection_width = .5, 
                   collapseTranscripts = "longest" )
)
dev.off()



