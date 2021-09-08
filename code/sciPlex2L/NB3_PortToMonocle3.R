# This script takes the cell x tile matrix from the filtered 
# ArchR package and ports it to a monocle cds object for 
# dimensionality reduction and UMAP embedding. 
# set working environment in R

basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190521_scichem2_AllPlates/"
out_dir =paste0(basepath, "analysis/archr_revised/")
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

addArchRGenome("hg19")
prj = loadArchRProject(path = paste0(out_dir,"sc2_filtered/"),
                       showLogo = FALSE)


###############################################
# create monocle UMAP from bin x cell matrix
###############################################

getAvailableMatrices(prj)

# load PeakMatrix into memory
bMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"sc2_filtered/ArrowFiles/sc2.arrow"),
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
# filter out unused cells 
cds_b = cds_b[,!is.na(colData(cds_b)$treatment)]

cd = data.frame(colData(cds_b)) %>% 
  mutate(treatment_f = ifelse(vehicle, "Vehicle", treatment))
colData(cds_b)$treatment_f = cd$treatment_f

# preprocess bin cds
# reduce dimensions
set.seed(2017)
cds_pl <- detect_genes(cds_b)
# restrict to bins accessible in >0.1% of cells
ncells = ncol(cds_b)*0.005
cds_pl = cds_pl[rowData(cds_pl)$num_cells_expressed > ncells,]
cds_pl <- estimate_size_factors(cds_pl)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
#cds_pl = align_cds(cds_pl, preprocess_method = "LSI", residual_model_formula_str = "~log(nFrags)+TSSEnrichment")
#cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "Aligned")
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# custom plot
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))
TCdat_pr$dose_character <- as.character(TCdat_pr$Relative_dose)
TCdat_pr$dose_character <- factor(TCdat_pr$dose_character,
                                  levels = c("0", "0.1", "0.5", "1", "5","10", "50", "100")) 


#pdf by clusters
pdf("Monocle3_UMAP_bins_cluster.pdf", width = 2, height = 1.25)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  #geom_point(size = 0.2, stroke = 0) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts()
dev.off()

#pdf by drug
pdf("Monocle3_UMAP_bins_drug2.pdf", width = 2, height = 1.25)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = treatment_f)) +
  #geom_point(size = 0.2, stroke = 0) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  scale_color_manual("Treatment",
                     labels = c("BMS" = "BMS345541", "Dex" = "Dex",
                                "Nutlin" = "Nutlin3A", "SAHA" = "SAHA", "Vehicle" = "Vehicle"),
                     values = c("BMS" = "darkmagenta", "Dex" = "deepskyblue3",
                                "Nutlin" = "firebrick3", "SAHA" = "springgreen4", "Vehicle" = "dimgrey")) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts()
dev.off()

# faceted pdf by drug and colored by dose
pdf("Monocle3_UMAP_bins_byDose_facet.pdf", width = 2, height = 1.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = dose_character)) +
  #geom_point(data = TCdat_pr, color = "gray", size = 0.05, stroke = 0) +
  geom_point_rast(data = TCdat_pr, size = 0.75, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  facet_wrap(~treatment, ncol = 2) +
  theme(legend.position = "right", text = element_text(size = 6), 
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  scale_color_manual("Dose",values = c("0"="gray", "0.1"="#1D1147FF", "0.5"="#51127CFF", "1"="#822681FF", "5"="#B63679FF",
                                       "10"="#E65164FF", "50" = "#FB8861FF", "100"="#FEC287FF")) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()


# faceted pdf by replicate and colored by drug
pdf("Monocle3_UMAP_bins_byDrug_facetReplicates2.pdf", width = 2, height = 1.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = treatment_f)) +
  #geom_point(data = TCdat_pr, color = "gray", size = 0.05, stroke = 0) +
  geom_point_rast(data = TCdat_pr, size = 0.75, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  facet_wrap(~replicate, ncol = 2) +
  theme(legend.position = "right", text = element_text(size = 6), 
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  scale_color_manual("Treatment",
                     labels = c("BMS" = "BMS345541", "Dex" = "Dex",
                                "Nutlin" = "Nutlin3A", "SAHA" = "SAHA", "Vehicle" = "Vehicle"),
                     values = c("BMS" = "darkmagenta", "Dex" = "deepskyblue3",
                                "Nutlin" = "firebrick3", "SAHA" = "springgreen4", "Vehicle" = "dimgrey")) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()



#######################################
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


saveArchRProject(ArchRProj = prj, 
                 load = FALSE)



# test plotting with ArchR
p <- plotEmbedding(ArchRProj = prj,
                   colorBy = "cellColData",
                   name = "clusters_mon3",
                   embedding = "UMAP_mon3")

plotPDF(p, 
        name = "ArchRPlot_Mon3_UMAP_bins_cluster.pdf", 
        addDOC = FALSE, width = 5, height = 5)



