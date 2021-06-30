# This script takes the cell x tile matrix from the filtered 
# ArchR package and ports it to a monocle cds object for 
# dimensionality reduction and UMAP embedding. 
# set working environment in R

basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB3"))
setwd(paste0(out_dir, "results/NB3"))

suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered"), 
                           showLogo = FALSE) 
###############################################
# create monocle UMAP from bin x cell matrix
###############################################

getAvailableMatrices(prj)

# load PeakMatrix into memory
bMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"mlr_filtered/ArrowFiles/MLR.arrow"),
                          ArchRProj = prj, 
                          cellNames = prj$cellNames,
                          useMatrix = "TileMatrix", 
                          binarize = TRUE)

# format rowData 
rd = data.frame(rowData(bMat)) %>% 
  dplyr::mutate(stop = start + 499,
                bin = paste(seqnames, start, stop, sep = "_")) %>% 
  dplyr::select(bin)
row.names(rd) <- rd$bin
row.names(bMat) <- rd$bin

# Create CDS from peak Matrix (SummarizedExperiment)
cds_b = monocle3::new_cell_data_set(assays(bMat)$TileMatrix, 
                                    cell_metadata = colData(bMat),
                                    gene_metadata = rd)

###############################################
# preprocess monocle3 cds
###############################################

# reduce dimensions
set.seed(2017)
cds_pl <- detect_genes(cds_b)
ncells = ncol(cds_b)*0.002
cds_pl = cds_pl[rowData(cds_pl)$num_cells_expressed > ncells,]
cds_pl <- estimate_size_factors(cds_pl)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
#cds_pl = align_cds(cds_pl, preprocess_method = "LSI", residual_model_formula_str = "~log(nFrags)")
reducedDim(cds_pl) <- reducedDim(cds_pl)[,2:50]
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# prepare relevant custom UMAPs 
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))

# by monocle cluster
pdf("Monocle3_bMat_UMAP_clusters.pdf", width = 2, height = 1.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()

# by Stimulator
pdf("Monocle3_bMat_UMAP_stimulator.pdf", width = 2, height = 1.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = Stimulator)) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()

# by Responder
pdf("Monocle3_bMat_UMAP_responder.pdf", width = 2, height = 1.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = Responder)) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()

# add monocle UMAP embedding to ArchR project
UMAP_mon3 = data.frame(reducedDims(cds_pl)$UMAP)
colnames(UMAP_mon3) <- c("LSI#UMAP_dim1", "LSI#UMAP_dim2")
prj@embeddings$UMAP_mon3[[1]] <- UMAP_mon3

# add monocle3 clusters to cellColData
prj <- addCellColData(ArchRProj = prj, 
                      data = as.character(TCdat_pr$cluster),
                      cells = row.names(TCdat_pr), 
                      name = "clusters_mon3", 
                      force = TRUE)

# test plotting with ArchR
p <- plotEmbedding(ArchRProj = prj,
                   colorBy = "cellColData",
                   name = "clusters_mon3",
                   embedding = "UMAP_mon3")

plotPDF(p, 
        name = "ArchRPlot_Mon3_UMAP_bins_cluster.pdf", 
        addDOC = FALSE,
        width = 5,
        height = 5)

saveArchRProject(ArchRProj = prj, 
                 load = FALSE)











