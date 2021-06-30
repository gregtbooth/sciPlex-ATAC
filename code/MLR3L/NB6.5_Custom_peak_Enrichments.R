# Greg 02/23/2021

# This Code is intended to provide a framework for testing for enrichment of any 
# set of sequence features within peaks of accessibility called using ArchR.

# for example: Testing weather individual TFBS are "enriched" in ATAC peaks specific to 
# each cluster 

# This framework is heavily built around ArchR containers. 
# This works by first taking custom elements as genomic coordinates, then building a 
# peak x element  matrix (using ArchR). Then running a regression for each element 
# to test whether it is significantly associated with cluster specific peaks 
# (a binary response, separate for each cluster or cell group)

basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB6"))
setwd(paste0(out_dir, "results/NB6"))

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(monocle3)
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

prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 

###########################################
# Generalizable functions after "annotating" 
# Archr project with data.
###########################################
# Function for pulling peak x element matrix from saved ArchR project
pull_peak_by_feature_mat = function(archr_project, file = "Motif-Matches-In-Peaks.rds"){
  Ann = readRDS(paste0(getOutputDirectory(archr_project), "/Annotations/", file))
  pf_mat = assay(Ann)
  peaks = data.frame(rowRanges(Ann)) %>% # GB Corrected 11/17/2020
    dplyr::mutate(peak = paste(seqnames, start, end, sep = "_"))
  row.names(pf_mat) <- peaks$peak
  return(pf_mat)
}


# Function:
# fit a linear regression to determine if features predict binary peak behavior.  
# each peak is an observation, features are predictors and the response variable is binary classification 
# for each peak (i.e. sig opening vs not). This gives a coef estimate for each feature (i.e. motif)
get.feature.logistic.reg.stats = function(response, features, data) {
  t(sapply(features, function(feature) {
    res = summary(glm(paste(response, feature, sep=" ~ "), family = "binomial", data = data))
    if (nrow(res$coefficients) > 1) {
      return(c(res$coefficients[2,1], res$coefficients[2,4]))
    } else {
      return(c(NA, NA))
    }
  }))
}

# function tests weather DA peaks are enriched for features 
feature_enrichment = function(
  archr_project,
  marker_test_results, 
  feature_mat,
  grouping = "Monocyte", 
  pVal_cutoff = 0.05)
{
  # Deafine peaks as DA or unchanged for cell type 
  peak_coefs = 
    marker_test_results %>%
    dplyr::filter(cell_group == grouping) %>% 
    dplyr::mutate(direction = ifelse(marker_test_p_value < pVal_cutoff, "Opening",  "Unchanged"))
 
  # filter peak x feature matrix for only used peaks (in DA analysis)
  feature_mat_f = feature_mat[peak_coefs$gene_id,]
  
  # list of feature names
  features = colnames(feature_mat_f)
  
  # replace any characters that will screw up formulas
  features = stringr::str_replace_all(features, "[+-.*]", "_")
  
  colnames(feature_mat_f) <- features
  
  # reformat from logical to numeric
  fmf = matrix(as.numeric(feature_mat_f), nrow = nrow(x=feature_mat_f))
  row.names(fmf) = row.names(feature_mat_f)
  colnames(fmf) = colnames(feature_mat_f)
  feature_df = as.data.frame(fmf)
  feature_df$gene_id = row.names(feature_df)
  # add direction info from coef_table
  feature_df = dplyr::inner_join(peak_coefs, feature_df, by = "gene_id")
  # add binary columns (dummy variables) describing directionality of each site 
  feature_df = cbind(feature_df, as.data.frame(model.matrix(~ 0 + direction, feature_df)))
  
  cat("fit linear regression predicting DA direction from motif presence. \nThis will take a few minutes...\n")
  
  # set container
  peak.DA.stats = list()
  
  if(dplyr::filter(feature_df, direction == "Opening") %>% 
     nrow > 0){
    peak.DA.stats$opening = get.feature.logistic.reg.stats(
      "directionOpening", features, data = feature_df)
    
    peak.DA.stats$opening = data.frame(
      feature = row.names(peak.DA.stats$opening),
      beta = peak.DA.stats$opening[,1],
      p.val = peak.DA.stats$opening[,2]) %>%
      mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
      dplyr::filter(!is.na(Padj), beta > 0) %>% 
      dplyr::arrange(Padj) %>% 
      dplyr::mutate(rank = rank(Padj), 
                    sig = ifelse(Padj <0.05, 1, 0), 
                    mlog10Padj = -log10(Padj)) 
  }
  return(peak.DA.stats$opening) 
}

# Automate the enrichment analysis across multiple groups of cells
# outputs list of feature enrichment tables with one table for each cell group
enrichment_test = function(
  Archr_Project = prj,
  Marker_test_results = marker_test_res, # i.e. fit_models or find_marker_genes results from monocle3
  Feature_mat = Encode_mat, # which features to run regressions on 
  PVal_cutoff = 0.05, 
  samples = c("Monocyte", "Activated_Tcell") # groups used to find DA sites in results from fit_models or find_marker_genes
)
{
  res = list()
  for (s in samples){
    res[[s]] = feature_enrichment(
      archr_project = Archr_Project, 
      marker_test_results = Marker_test_results,
      feature_mat = Feature_mat, # which features to run regressions on 
      grouping = s,  # which cell type specific peaks 
      pVal_cutoff = PVal_cutoff
    )
  }
  res
}

#####
# function for making heatmap matrix of feature enrichments for each cell group
# similar to ArchRs plotEnrichHeatmap function

make_heatmap_mat = function(enrichment_res, top_n = 10, cutoff = 0.05){
  # prep container for results 
  res = data.frame(feature = character(), beta = numeric(), p.val = numeric(), 
                   Padj = numeric(), rank = numeric(), sig = numeric(), 
                   cellgroup = character())
  for (i in names(enrichment_res)){
    enrichment_res[[i]][["cellgroup"]] = i
    res = rbind(res, enrichment_res[[i]])
  }
  # isolate features enriched in at least one group
  res_f = filter(res, rank <= top_n, Padj <= cutoff)
  # keep values for important features across groups
  res_ff = filter(res, feature %in% res_f$feature)
  # reformat as matrix with -log10(pAdj) 
  res_mat = acast(res_ff, cellgroup~feature, value.var = "mlog10Padj", fill = 0)
  res_mat 
}

#####
# Function for creatign a ComplexHeatmap object for plotting
create_heatmap = function(heatmap_matrix){
  ComplexHeatmap::Heatmap(
    matrix = heatmap_matrix, 
    name = "-log10(adj pval)", 
    col = paletteContinuous(set = "comet", n = 100),
    heatmap_legend_param = list(
      # at = c(0, 1), 
      #labels = c(round(min(limits), 2), round(max(limits), 2)),
      #labels = c(round(0, 2), round(100, 2)),
      # color_bar = "continuous",
      legend_direction = "horizontal", 
      legend_width = unit(3, "cm")), 
    rect_gp = gpar(col = NA), 
    show_column_names = TRUE, 
    column_names_side = "top",
    cluster_columns = TRUE, 
    show_column_dend = FALSE, 
    clustering_method_columns = "ward.D2", 
    column_names_gp = gpar(fontsize = 6), 
    column_names_max_height = unit(100, "mm"), 
    show_row_names = TRUE, 
    row_names_gp = gpar(fontsize = 8), 
    cluster_rows = TRUE, 
    show_row_dend = FALSE,
    clustering_method_rows = "ward.D2", 
    split = NULL, 
    #top_annotation = ht1Anno, 
    use_raster = TRUE, 
    raster_device = "png",
    raster_quality = 5)
}

####################################
# Encode TFBS (chipSeq) enrichment
###################################
# Load required data 

# read in pre-calculated marker_test_res (from monocle3 top_markers)
marker_test_res = read.table(file = "marker_test_res_Top25Kpeaks.txt", head = TRUE)

# "Annotate" my peaks (find overlaps) with ArchR's stored Encode ChIP-seq peaks 
# only need to do once 
#prj <- addArchRAnnotations(ArchRProj = prj, collection = "EncodeTFBS")

# pull peak x encode TFBS matrix 
Encode_mat = pull_peak_by_feature_mat(archr_project = prj, file = "EncodeTFBS-Matches-In-Peaks.rds")
# correct_formatting of encode colnames
cn = paste0(stringr::str_split_fixed(colnames(Encode_mat), pattern = "\\.", n = 3)[,2],
            "_",
            stringr::str_split_fixed(colnames(Encode_mat), pattern = "\\.", n = 3)[,1])
colnames(Encode_mat)<- cn


#####
# Which ChIPseq profiles explain Top Marker peaks
Encode_TFBS_enrichments = enrichment_test(
  Archr_Project = prj,
  Marker_test_results = marker_test_res,
  Feature_mat = Encode_mat, # which features to run regressions on 
  PVal_cutoff = 0.05, 
  samples = c("Monocyte", "Activated_Tcell", "Bcell", "Naive_Tcell", "NK")
)

# create pValue enrichment matrix
encode_heatmap_mat = make_heatmap_mat(
  enrichment_res = Encode_TFBS_enrichments, 
  top_n = 10, 
  cutoff = 0.05)

# creat heatmap object for plotting
encode_heatmap = create_heatmap(encode_heatmap_mat)

# plot
pdf(file = "heatmap_Encode_TFBS_enrichment.pdf", height = 3, width = 4)
ComplexHeatmap::draw(encode_heatmap, heatmap_legend_side = "bottom")
dev.off()


####################################
# motif enrichment
###################################

# pull peak x motif matrix 
motif_mat = pull_peak_by_feature_mat(archr_project = prj, file = "Motif-Matches-In-Peaks.rds")

#####
# Which ChIPseq profiles explain Top Marker peaks
motif_enrichments = enrichment_test(
  Archr_Project = prj,
  Marker_test_results = marker_test_res,
  Feature_mat = motif_mat, # which features to run regressions on 
  PVal_cutoff = 0.05, 
  samples = c("Monocyte", "Activated_Tcell", "Bcell", "Naive_Tcell", "NK")
)

# create pValue enrichment matrix
motif_heatmap_mat = make_heatmap_mat(
  enrichment_res = motif_enrichments, 
  top_n = 10, 
  cutoff = 0.05)

# creat heatmap object for plotting
motif_heatmap = create_heatmap(motif_heatmap_mat)

# plot
pdf(file = "heatmap_motif_enrichments.pdf", height = 3, width = 4)
ComplexHeatmap::draw(motif_heatmap, heatmap_legend_side = "bottom")
dev.off()


####################################
# Bulk Track (FACS sorted ATAC) enrichment
###################################
bulk_mat = pull_peak_by_feature_mat(archr_project = prj, file = "ATAC-Matches-In-Peaks.rds")

#####
# Which ChIPseq profiles explain Top Marker peaks
bulk_enrichments = enrichment_test(
  Archr_Project = prj,
  Marker_test_results = marker_test_res,
  Feature_mat = bulk_mat, # which features to run regressions on 
  PVal_cutoff = 0.05, 
  samples = c("Monocyte", "Activated_Tcell", "Bcell", "Naive_Tcell", "NK")
)

# create pValue enrichment matrix
bulk_heatmap_mat = make_heatmap_mat(
  enrichment_res = bulk_enrichments, 
  top_n = 10, 
  cutoff = 0.05)

# creat heatmap object for plotting
bulk_heatmap = create_heatmap(bulk_heatmap_mat)

# plot
pdf(file = "heatmap_BulkATAC_enrichments.pdf", height = 3, width = 4)
ComplexHeatmap::draw(motif_heatmap, heatmap_legend_side = "bottom")
dev.off()














##################################################
# function for plotting  motifs ranked by Pvals
plot_motif_enrichment = function(motif_test_results, cell_type = "Monocyte")
{
  ggUp <- ggplot(motif_test_results, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point_rast(size = 1) +
    ggrepel::geom_label_repel(
      data = motif_test_results[1:5, ],
      aes(x = rank, y = mlog10Padj, label = motif), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + 
    theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
  ggsave(plot = ggUp, 
         filename = paste0(cell_type, "_motifs_openingSites_rast_Padj.pdf"),
         width = 3.75, 
         height = 3.75)
  
}













# run above functions for each cell type
for(ct in c("Activated_Tcell", "Bcell", "Naive_Tcell", "NK", "Monocyte")){
  # run motif enrichment regressions
  monocyte_motif_test_res = 
    celltype_motif_enrichment(
      archr_project = prj,
      marker_test_results = marker_test_res, 
      cell_type = ct,
      pVal_cutoff = 0.005)
  # plot ranked motif enrichments
  plot_motif_enrichment(monocyte_motif_test_res, cell_type = ct)
}





