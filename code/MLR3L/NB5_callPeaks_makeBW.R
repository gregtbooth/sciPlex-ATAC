## on cluster, initiate qlogin session with 1 core. 
# qlogin -q trapnell-login.q -l mfree=50G -pe serial 1
# 
#load the following modules 
#module load MACS/2.2.7.1 
#
# start R (4.0.0)

#set working environment in R
basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB5"))
setwd(paste0(out_dir, "results/NB5"))

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

###############################
#Only need to run this part once
###############################
# rerun iterative LSI on filtered annotated cells
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

prj <- addImputeWeights(ArchRProj = prj)

# plot ArchR embedding
p <- plotEmbedding(ArchRProj = prj,
                   colorBy = "cellColData",
                   name = "Clusters", 
                   embedding = "UMAP")
plotPDF(p, 
        name = "Plot-archrUMAP-clusters.pdf", 
        addDOC = FALSE, width = 5, height = 5)



###############################################
# re process filtered cells with monocle3 
# create monocle UMAP from bin x cell matrix
###############################################

getAvailableMatrices(prj)

# load PeakMatrix into memory
bMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"mlr_filtered_annotated/ArrowFiles/MLR.arrow"),
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
pdf("Monocle3_bMat_UMAP_cellTypes.pdf", width = 2, height = 1.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = cellType_broad)) +
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

# test plotting with ArchR
p <- plotEmbedding(ArchRProj = prj,
                   colorBy = "cellColData",
                   name = "cellType_broad",
                   embedding = "UMAP_mon3")

plotPDF(p, 
        name = "ArchRPlot_Mon3_UMAP_bins_cellType_broad.pdf", 
        addDOC = FALSE,
        width = 5,
        height = 5)

saveArchRProject(ArchRProj = prj, 
                 load = FALSE)

##################################################
##################################################

#Once above has been run and saved,
#start peak calling here

###############################################
# Call Peaks
###############################################
prj = loadArchRProject(
  path = paste0(out_dir, "mlr_filtered_annotated"), 
  showLogo = FALSE) 

# make pseduo bulk replicates for monocle3-defined clusters
prj = addGroupCoverages(
  ArchRProj = prj,
  groupBy=  "cellType_broad",
  useLabels = TRUE,
  minCells= 50,
  maxCells = 500,
  minReplicates = 2,
  maxReplicates = 5,
  sampleRatio = 0.8,  
  threads =  getArchRThreads(),
  verbose =  TRUE
)

# call peaks. 
prj <- addReproduciblePeakSet(
  ArchRProj = prj, 
  groupBy = "cellType_broad", 
  peakMethod = "Macs2",
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 25,
  excludeChr = c("chrM", "chrY"), 
  extendSummits = 250,
  promoterRegion = c(2000, 100),
  genomeAnnotation = getGenomeAnnotation(prj),
  geneAnnotation = getGeneAnnotation(prj),
  plot = TRUE,
  threads = getArchRThreads(),
  force = TRUE
)

# add peak matrix to arrow file
prj <- addPeakMatrix(
  prj, 
  binarize = FALSE, # don't binarize here (can be binarized when retrieved)
  force = TRUE)

# add peak motif information
prj <- addMotifAnnotations(
  ArchRProj = prj, 
  motifSet = "cisbp", 
  name = "Motif")

# save output
saveArchRProject(
  ArchRProj = prj, 
  load = FALSE)


###################################
# create pseudobulk bigwig files for 
# all celltypes from all treatment conditions

prj$celltype_condition = paste0(prj$Responder, "_", prj$Stimulator, "_", prj$cellType_broad)

# separate tracks for all cell types within each condition (139 total)
getGroupBW(
  ArchRProj = prj,
  groupBy = "celltype_condition",
  normMethod = "ReadsInTSS", # creates a group scale factor = 10k/sum(reads in TSS) 
  #normMethod = "nFrags", # creates a group scale factor = 10k/sum(nFrags)
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

# separate tracks for each cell type (combined conditions)
getGroupBW(
  ArchRProj = prj,
  groupBy = "cellType_broad",
  normMethod = "ReadsInTSS", # creates a group scale factor = 10k/sum(reads in TSS) 
  #normMethod = "nFrags", # creates a group scale factor = 10k/sum(nFrags)
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)








