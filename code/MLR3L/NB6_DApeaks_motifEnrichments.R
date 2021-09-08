basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
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

#################################
# print bed file of peaks (for checking in IGV): 
peaks = data.frame(getPeakSet(prj)) %>% 
  dplyr::select( -width, -strand, -score, -replicateScoreQuantile, -groupScoreQuantile, 
                 -Reproducibility, -GroupReplicate)

write.table(peaks, file="peaks.bed", sep = "\t", quote=FALSE, row.names = FALSE)


#################################

# get marker Peaks for defined cell types
markersPeaks <- getMarkerFeatures(
  ArchRProj = prj, 
  useMatrix = "PeakMatrix", 
  groupBy = "cellType_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# sig sites 
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList

# make heat map.
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks,
        name = "Peak-Marker-Heatmap", 
        width = 8, 
        height = 6, 
        addDOC = FALSE)

##################################
# Try looking specifically for DA sites 
# between Naive and Activated T-cells

markerTest <- getMarkerFeatures(
  ArchRProj = prj, 
  useMatrix = "PeakMatrix",
  groupBy = "cellType_broad",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Activated_Tcell",
  bgdGroups = "Naive_Tcell"
)

pma <- plotMarkers(
  seMarker = markerTest,
  name = "Activated_Tcell", 
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5",
  plotAs = "MA")

pv <- markerPlot(
  seMarker = markerTest, 
  name = "Activated_Tcell",
  cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", 
  plotAs = "Volcano")

plotPDF(pma, pv, 
        name = "Activated-vs-Naive_Tcells-Markers-MA-Volcano",
        width = 5, 
        height = 5, 
        addDOC = FALSE)


##################################
# port peak x cell data back into monocle
##################################

# load PeakMatrix into memory
pMat = getMatrixFromArrow(
  ArrowFile = paste0(out_dir, "mlr_filtered_annotated/ArrowFiles/MLR.arrow"),
  ArchRProj = prj, 
  cellNames = prj$cellNames,
  useMatrix = "PeakMatrix",
  binarize = TRUE)

# Create CDS from peak Matrix (SummarizedExperiment)
rd = data.frame(rowRanges(pMat)) %>% 
  tibble::rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end))

pd = data.frame(getPeakSet(prj)) %>% 
  tibble::rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end), 
         gene_short_name = ifelse(peakType == "Promoter" & distToGeneStart < 500, nearestGene, "NA")) %>% 
  dplyr::select(rowname, score, replicateScoreQuantile, groupScoreQuantile, 
                Reproducibility, GroupReplicate, distToGeneStart,
                nearestGene, peakType, distToTSS, nearestTSS,
                gene_short_name, GC, idx, peak)

a_rd = left_join(rd, pd, by = "peak") %>% 
  dplyr::select(idx = idx.x, score, replicateScoreQuantile, groupScoreQuantile, 
         Reproducibility, GroupReplicate, distToGeneStart,
         nearestGene, peakType, distToTSS, nearestTSS,
         gene_short_name, GC, peak)

row.names(a_rd) <- a_rd$peak
row.names(pMat) <- a_rd$peak
cds_atac = monocle3::new_cell_data_set(assays(pMat)$PeakMatrix, 
                                       cell_metadata = colData(pMat), 
                                       gene_metadata = a_rd)

# preprocess cell x peak cds object: 
# reduce dimensions
set.seed(2017)
cds_pl <- detect_genes(cds_atac)
#ncells = ncol(cds_b)*0.002
#cds_pl = cds_pl[rowData(cds_pl)$num_cells_expressed > ncells,]
cds_pl <- estimate_size_factors(cds_pl)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
#cds_pl = align_cds(cds_pl, preprocess_method = "LSI", residual_model_formula_str = "~log(nFrags)")
reducedDim(cds_pl) <- reducedDim(cds_pl)[,2:50]
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# save cds object
saveRDS(cds_pl, file = "cds_process.rds")

###############################################
# start with saved processed cds
cds_pl = readRDS(file = "cds_process.rds")

# prepare relevant custom UMAPs 
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))

# by monocle cluster
pdf("Monocle3_peak-based_UMAP_cellTypes.pdf", width = 2.25, height = 1.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = cellType_broad)) +
  geom_point_rast(size=0.4, stroke = 0) + 
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()


############
# Using processed peak x cell matrix 
# run top_markers

marker_test_res <- top_markers(
  cds_pl,
  group_cells_by="cellType_broad",
  genes_to_test_per_group = 25000, 
  reduction_method = "UMAP",
  marker_sig_test = TRUE,
  reference_cells = 1000,
  speedglm.maxiter = 25,
  cores = 32,
  verbose = FALSE
)

write.table(
  marker_test_res, 
  file = "marker_test_res_Top25Kpeaks.txt",
  quote = FALSE, 
  sep = "\t", 
  row.names = FALSE, 
  col.names = TRUE)

#############################
#############################

# read in pre-calculated marker_test_res
marker_test_res = 
  read.table(file = "marker_test_res_Top25Kpeaks.txt", 
             head = TRUE)

# add relaxed q_values (based on 25K tests per cell type)
marker_test_res =
  marker_test_res %>% 
  group_by(cell_group) %>% 
  mutate(marker_test_q_value_relaxed = p.adjust(marker_test_p_value, n=25000))

# plot number of sig marker peaks per cell type 
# (qvalue < 0.05)
sig_peaks = marker_test_res %>% 
  filter(marker_test_q_value_relaxed < 0.05) %>% 
  group_by(cell_group) %>% 
  summarise(peaks = n())

ggplot(sig_peaks, aes(x = cell_group, y = peaks)) +
  geom_col() + 
  geom_text(aes(label = peaks), vjust = -0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 

ggsave(filename = "sig_peaks_perCelltype.pdf",
       width = 2.5, 
       height = 3)

# plot number of sig marker peaks per cell type
# (qvalue < 0.05)
sig_peaks_p = marker_test_res %>% 
  filter(marker_test_p_value < 0.005) %>% 
  group_by(cell_group) %>% 
  summarise(peaks = n())

ggplot(sig_peaks_p, aes(x = cell_group, y = peaks)) +
  geom_col() + 
  geom_text(aes(label = peaks), vjust = -0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) 

ggsave(filename = "sig_peaks_perCelltype_pvalue_005.pdf",
       width = 2.5, 
       height = 3)

###########################################

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

  
celltype_motif_enrichment = function(
  archr_project,
  marker_test_results, 
  cell_type = "Monocyte", 
  pVal_cutoff = 0.05)
{
  # get all coefficients for cell type
  peak_coefs = 
    marker_test_results %>%
    dplyr::filter(cell_group == cell_type) %>% 
    dplyr::mutate(direction = ifelse(marker_test_p_value < pVal_cutoff, "Opening",  "Unchanged"))
  
  # get peak x motif matrix 
  motif_mat = pull_peak_by_motif_mat(archr_project)
  
  # filter for only used peaks
  motif_mat_f = motif_mat[peak_coefs$gene_id,]
  
  # list of motif names
  motifs = colnames(motif_mat_f)
  
  # reformat from logical to numeric
  mmf = matrix(as.numeric(motif_mat_f), nrow = nrow(x=motif_mat_f))
  row.names(mmf) = row.names(motif_mat_f)
  colnames(mmf) = colnames(motif_mat_f)
  motif_df = as.data.frame(mmf)
  motif_df$gene_id = row.names(motif_df)
  # add direction info from coef_table
  motif_df = dplyr::inner_join(peak_coefs, motif_df, by = "gene_id")
  # add binary columns describing directionality of each site 
  motif_df = cbind(motif_df, as.data.frame(model.matrix(~ 0 + direction, motif_df)))
  
  cat("fit linear regression predicting DA direction from motif presence. \nThis will take a few minutes...\n")
  
  # set container
  peak.DA.stats = list()
  
  if(dplyr::filter(motif_df, direction == "Opening") %>% 
     nrow > 0){
    peak.DA.stats$promoter.opening = get.peak.DA.logistic.reg.stats(
      "directionOpening", motifs, motif_df)

    peak.DA.stats$promoter.opening = data.frame(
      motif = row.names(peak.DA.stats$promoter.opening),
      beta = peak.DA.stats$promoter.opening[,1],
      p.val = peak.DA.stats$promoter.opening[,2]) %>%
      mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
      dplyr::filter(!is.na(Padj), beta > 0) %>% 
      dplyr::arrange(Padj) %>% 
      dplyr::mutate(rank = rank(Padj), 
                    sig = ifelse(Padj <0.05, 1, 0), 
                    mlog10Padj = -log10(Padj)) 
  }
  return(peak.DA.stats$promoter.opening) 
}
  

# funciton for plotting  motifs ranked by Pvals
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


###################################
###################################
# try the same analysis as above, but using peak x cell-type bulk-profle matrix
# Want to assess how each cluster's accessibility relates to cell type defined 
# bulk profiles

# add default ArchR peak x bulk ATAC matrix (FACS sorted PBMC and BMMC cell types)
prj <- addArchRAnnotations(
  ArchRProj = prj, 
  collection = "ATAC")



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



# compare annotations with FACS sorted bulk ATAC profiles (from archr)
# Using the pre-filtered cells with annotations


# use ArchR to get cell-type specific marker peaks
markersPeaks <- getMarkerFeatures(
  ArchRProj = prj, 
  useMatrix = "PeakMatrix", 
  groupBy = "cellType_broad",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Add monocle (top_markers) data to peaks 



# filter peaks based on column criteria
markerList <- getMarkers(
  markersPeaks, 
  cutOff = "MeanDiff >= 0.25 & FDR <= 0.25"
)

markerList



# look for overlap of cell-type bulk tracks and cluster specific peaks
enrichATAC <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = prj,
  peakAnnotation = "ATAC",
  cutOff = "FDR <= 0.25 & Log2FC >= 0.5"
)


# make heatmap 
heatmapATAC <- plotEnrichHeatmap(
  enrichATAC,
  n = 7, 
  transpose = TRUE)

# plot heatmap
plotPDF(heatmapATAC, 
        name = "bulkATAC-Enriched-Marker-Heatmap", 
        width = 8, 
        height = 6, 
        addDOC = FALSE)

