basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB4"))
setwd(paste0(out_dir, "results/NB4"))

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
  library(RColorBrewer)
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered"), 
                       showLogo = FALSE) 

##############################
# run iterative LSI 
prj <- addIterativeLSI(ArchRProj = prj, 
                       useMatrix = "TileMatrix", 
                       name = "IterativeLSI", 
                       iterations = 2,
                       dimsToUse = 2:50,
                       varFeatures = 100000,
                       force = TRUE)

# Add clusters
prj <- addClusters(input = prj, 
                   name = "Clusters",
                   reducedDims = "IterativeLSI", 
                   method = "Seurat",
                   corCutOff = 0.75,
                   knnAssign = 10,
                   force = TRUE)

# Add UMAP embedding
prj <- addUMAP(ArchRProj = prj, 
               reducedDims = "IterativeLSI", 
               name = "UMAP",
               nNeighbors = 40,
               minDist = 0.4,
               metric = "cosine",
               force = TRUE)

# plot ArchR embedding
p <- plotEmbedding(ArchRProj = prj,
                   colorBy = "cellColData",
                   name = "Clusters", 
                   embedding = "UMAP")
plotPDF(p, 
        name = "Plot-archrUMAP-clusters.pdf", 
        addDOC = FALSE, width = 5, height = 5)


##############################
# plot imputed marker gene scores over 
# monocle embedding

prj <- addImputeWeights(ArchRProj = prj)

saveArchRProject(ArchRProj = prj, 
                 load = FALSE)


# plot imputed gene scores in monocles UMAP embedding. 
markerGenes  <- c(
  "IL7R", "CCR7", #Naive CD4+ T
  "S100A4", #Memory CD4+ T
  "CD14", "LYZ",  #CD14 + Mono
  "CEBPB", "MPO","IRF8", #Monocytes
  "FCGR3A", "MS4A7", #FCGR3A+ Mono
  "CD8A", #CD8+ TCells
  "CD3D",  "TBX21", #TCells
  "MS4A1", # Bcells
  "GNLY", "NKG7", #NK
  "FCER1A", "CST3", #FCER1A + monocytes
  "PPBP", #platelet (A-nuclear, will not show up in ATAC)
  "CCNB1", "CCNB2","CDK1" #proliferation markers
)

p <- plotEmbedding(
  ArchRProj = prj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP_mon3",
  imputeWeights = getImputeWeights(prj)
)

plotPDF(p, 
        name = "Plot-monocle3_UMAP-Marker-Genes-W-Imputation.pdf", 
        addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(
  ArchRProj = prj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP_mon3",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

plotPDF(p1, 
        name = "Plot-monocle3_UMAP-Marker-Genes-RAW.pdf", 
        addDOC = FALSE, width = 5, height = 5)

#######################################################
# find "marker genes" which differ across monocle clusters

markersGS <- getMarkerFeatures(
  ArchRProj = prj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "clusters_mon3",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1.25")

# plot heatmap of maker fold change (gene scores) by cluster:
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

plotPDF(heatmapGS,
        name = "ArchR_GeneScores-Marker-Heatmap", 
        width = 8,
        height = 6, 
        addDOC = FALSE)


####################################################
# count the number and proportion of cells 
# within each cluster from each condition 
cdat = data.frame(prj@cellColData) %>% 
  tibble::rownames_to_column(var = "cellNames") 

# count treatment cells in clusters
cell_counts_summary = group_by(cdat, Responder, Stimulator, clusters_mon3) %>% 
  summarise(n_cells = n()) 

# stacked bar of proportion stimulators per cluster
cell_counts_summary %>%
  dplyr::select(Responder,Stimulator, clusters_mon3, n_cells) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = clusters_mon3, y =n_cells, fill = Stimulator), 
           color = "black", size = .25, stat = "identity", 
           position = "fill") +
  facet_wrap(~Responder) +
  #scale_fill_manual("Dose (uM)", values = clr_values) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  #scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Cluster") +
  ylab("Proportion") +
ggsave("stackedbarplot_stimulatorPerclusters.pdf",
        width = 3, height = 2.25, unit = "in", dpi = 1200)

# stacked bar of proportion cluster per stimulator
cell_counts_summary %>% 
  #filter(Responder == "stimAlone") %>% 
  dplyr::select(Responder,Stimulator, clusters_mon3, n_cells) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = Stimulator, y =n_cells, fill = clusters_mon3), 
           color = "black", size = .25, stat = "identity", 
           position = "fill") +
  facet_wrap(~Responder) +
  #scale_fill_manual("Dose (uM)", values = clr_values) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  #scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Stimulator") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("stackedbarplot_clustersPerStimulator.pdf",
        width = 3, height = 2.5, unit = "in", dpi = 1200)
#ggsave("stackedbarplot_clusters_StimAlone.pdf",
#       width = 2, height = 1.25, unit = "in", dpi = 1200)



####################################################
# integrate ATAC with annotated scRNA-seq data (10x)
####################################################
# set up project for integration 
prj_int = prj

# load 10X PBMC data (Downloaded from Seurat page)
pbmc.rna <- readRDS("~/shared/nobackup/Seurat/10K_PBMC_scRNAseq_10X/pbmc_10k_v3.rds")

# unconstrained integration of scRNA with scATAC
# not adding to Arrow, but just for assessing assignments
prj_int <- addGeneIntegrationMatrix(
  ArchRProj = prj_int, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = pbmc.rna,
  addToArrow = FALSE,
  groupRNA = "celltype",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un", 
  force = TRUE
)

cM <- as.matrix(confusionMatrix(prj_int$clusters_mon3, prj_int$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

# transform to fraction of rowSums (total cells/condition)
cM_frac <- apply(cM, 2, function(i) i/sum(i))
cM_melt = melt(cM_frac)
colnames(cM_melt) = c("Cluster", "predictionBased", "fracPredicted")

# plot confusion matrix
ggplot(cM_melt, aes(x = Cluster, y = predictionBased, fill = fracPredicted)) + 
  geom_tile(color = "gray") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  ylab("cell-type prediction") +
  xlab("mon3_cluster") + 
  labs(fill = "frac. predicted") 
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "Confusion_predCellTypeVsCluster_unconstrained.pdf",
       width = 6, height = 3.5)

# report the median prediction scores for each cluster
data.frame(prj_int@cellColData) %>% 
  group_by(clusters_mon3) %>% 
  summarise(med_pscore = median(predictedScore_Un))

# create a heat map of prediction scores for each pairing
pscores = data.frame(prj_int@cellColData) %>% 
  group_by(predictedGroup_Un, clusters_mon3) %>% 
  summarise(med_pscore = median(predictedScore_Un))
colnames(pscores) = c("predCellType", "Cluster", "med_predScore")

ggplot(pscores, aes(x = Cluster, y = predCellType, fill = med_predScore)) + 
  geom_tile(color = "gray") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  ylab("cell-type prediction") +
  xlab("mon3_cluster") + 
  labs(fill = "Med. Prediction Score") 
#theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "Confusion_CellTypeVsCluster_pred_scores.pdf",
       width = 6, height = 3.5)


##############################

# plot umap colored by predictedcell-type annotations 
p1 <- plotEmbedding(
  ArchRProj = prj_int, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  embedding = "UMAP_mon3",
)

plotPDF(p1, 
        name = "Plot-monocle3_UMAP-predictedCellTypes.pdf", 
        addDOC = FALSE, width = 5, height = 5)


# add metadata column that only 
# labels cells with high prediction scores 
cdat = data.frame(prj_int@cellColData) %>% 
  tibble::rownames_to_column(var = "cellNames") %>% 
  mutate(cellType_pred_highscore = ifelse(predictedScore_Un > 0.5, 
                                          predictedGroup_Un, "low_conf"))

prj_int <- addCellColData(
  ArchRProj = prj_int, 
  data = cdat$cellType_pred_highscore,
  cells = cdat$cellNames, 
  name = "cellType_pred_passQC",
  force = TRUE)

# replot UMAP with high-confidence celltype calls
p2 <- plotEmbedding(
  ArchRProj = prj_int, 
  colorBy = "cellColData", 
  name = "cellType_pred_passQC", 
  embedding = "UMAP_mon3",
)

plotPDF(p2, 
        name = "Plot-monocle3_UMAP-predictedCellTypes_passQC.pdf", 
        addDOC = FALSE, width = 5, height = 5)


###########################
# use cell type predictions and clusters 
# to establish broad cell annotations for 
# downstream analysis. 

cdat = data.frame(prj_int@cellColData) %>% 
  tibble::rownames_to_column(var = "cellNames") %>% 
  mutate(celltype_broad = 
           recode(clusters_mon3, 
                  "1" = "Naive_Tcell", 
                  "2" = "NK",
                  "3" = "Dead",
                  "4" = "Monocyte",
                  "5" = "Activated_Tcell",
                  "6" = "Dead", 
                  "7" = "Bcell"))

prj <- addCellColData(
  ArchRProj = prj, 
  data = cdat$celltype_broad,
  cells = cdat$cellNames, 
  name = "cellType_broad",
  force = TRUE)

# Save annotation info to project
saveArchRProject(ArchRProj= prj,
                 load = FALSE)


############################################
############################################
# start from saved project with info. 
prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered"), 
                       showLogo = FALSE) 

# plot UMAP with final annotations
p3 <- plotEmbedding(
  ArchRProj = prj, 
  colorBy = "cellColData", 
  name = "cellType_broad", 
  embedding = "UMAP_mon3",
)

plotPDF(p3, 
        name = "Plot-monocle3_UMAP-annotatedCellTypes_final.pdf", 
        addDOC = FALSE, width = 5, height = 5)

# custom UMAP plots with annotations 
cdat = data.frame(prj@cellColData) %>% 
  tibble::rownames_to_column(var = "Cell")
# retrieve stored monocle3 UMAP coords
umap_coords = data.frame(prj@embeddings$UMAP_mon3[[1]]) %>% 
  tibble::rownames_to_column(var = "Cell") %>% 
  rename(UMAP1 = LSI.UMAP_dim1, UMAP2 = LSI.UMAP_dim2)

cdat = left_join(cdat, umap_coords, by = "Cell")

# plot by broad annotations
pdf("UMAP_cellTypes_noLegend.pdf", width = 1.75, height = 1.75)
ggplot(cdat, aes(x = UMAP1, y = UMAP2, color = cellType_broad)) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  #theme(legend.position = "right", text = element_text(size = 6),  
  #      legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  theme(legend.position = "none") +
  #scale_color_brewer(palette = "Dark2") +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()



# Remake stacked barplots with annotations
cell_counts_summary = data.frame(prj@cellColData) %>% 
  group_by( Responder, Stimulator, cellType_broad) %>% 
  summarise(n_cells = n()) 

cell_counts_summary %>%
  dplyr::select(Responder,Stimulator, cellType_broad, n_cells) %>%
  distinct() %>%
  filter(Responder != "stimAlone") %>% 
  ggplot() +
  geom_bar(aes(x = Stimulator, y =n_cells, fill = cellType_broad), 
           color = "black", size = .25, stat = "identity", 
           position = "fill") +
  facet_wrap(~Responder) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  xlab("Stimulator") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("stackedbarplot_cellTypesPerStimulator_noStimAlone.pdf",
       width = 2.5, height = 2.2, unit = "in", dpi = 1200)


# now just plot the stim alone condition (for supplementary)
cell_counts_summary %>%
  dplyr::select(Responder,Stimulator, cellType_broad, n_cells) %>%
  distinct() %>%
  filter(Responder == "stimAlone") %>% 
  ggplot() +
  geom_bar(aes(x = Stimulator, y =n_cells, fill = cellType_broad), 
           color = "black", size = .25, stat = "identity", 
           position = "fill") +
  facet_wrap(~Responder) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  xlab("Stimulator") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("stackedbarplot_cellTypesPerStimulator_onlyStimAlone.pdf",
       width = 1.75, height = 1.25, unit = "in", dpi = 1200)



###################################### 
# Using Replicate information, plot cell type
# proportions for each condition. 
# want to look for significant changes 
# (i.e. activated Tcells after bead stimulation)
cell_counts_condition = data.frame(prj@cellColData) %>% 
  group_by(Replicate, Responder, Stimulator) %>% 
  summarise(total_cells = n()) %>% 
  mutate(condition = paste(Responder, Stimulator, Replicate, sep = "_"))

cell_counts_summary = data.frame(prj@cellColData) %>% 
  group_by(Replicate, Responder, Stimulator, cellType_broad) %>% 
  summarise(n_cells = n()) %>% 
  mutate(condition = paste(Responder, Stimulator, Replicate, sep = "_")) %>% 
  left_join(cell_counts_condition, by = "condition") %>% 
  select(
    Replicate = Replicate.x,
    Responder = Responder.x, 
    Stimulator = Stimulator.x, 
    cellType_broad, 
    n_cells, 
    total_cells) %>% 
  mutate(percent_cellType = (n_cells/total_cells)*100)

# make box plots of cell proportion faceted by responder. 
# prepare separate plots for each cell type 

celltypes = cell_counts_summary %>% 
  distinct(cellType_broad)

for (ct in celltypes$cellType_broad){
  pdf(file = paste0("BoxPlot_PercentRecoveredCellType_byStimulator_", ct,".pdf"), 
      width = 4, height = 3)
  print(  
  cell_counts_summary %>% 
    filter(cellType_broad == ct) %>% 
    ggplot() +
      geom_boxplot(aes(x =  Stimulator, y = percent_cellType)) +
      facet_wrap("~Responder") +
      scale_color_brewer(palette='Set1') +
      xlab("Stimulator") +
      ylab("cells") +
      ggtitle(paste0("Percent recovered ", ct)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    #theme_bw()
  )
  dev.off()
}


# filter out dead cells and make new project: 
prj_celltypes = prj[prj$cellType_broad != "Dead",]

saveArchRProject(ArchRProj= prj_celltypes,
                 outputDirectory = paste0(out_dir, "mlr_filtered_annotated"),
                 dropCells = TRUE,
                 load = FALSE)


