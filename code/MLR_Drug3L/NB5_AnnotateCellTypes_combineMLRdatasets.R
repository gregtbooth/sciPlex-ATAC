basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/211103_3Level_scichem_MLR_Drugs/"
basepath_oldData = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
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
  library(RColorBrewer)
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

# Load new (Drug) data as archr project
prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered"), 
                       showLogo = FALSE) 

# Load old (MLR) data as archr project
prj_o = loadArchRProject(path = paste0(basepath_oldData, "analysis/archr/mlr_filtered"), 
                       showLogo = FALSE) 


#############################################
# get bMat from new data
# load PeakMatrix into memory
bMat = getMatrixFromProject( ArchRProj = prj, 
                             useMatrix = "TileMatrix", 
                             binarize = TRUE)

# format rowData 
rd = data.frame(rowData(bMat)) %>% 
  dplyr::mutate(stop = start + 499,
                bin = paste(seqnames, start, stop, sep = "_")) %>% 
  dplyr::select(bin)
row.names(rd) <- rd$bin
row.names(bMat) <- rd$bin

# format colData 
cdat = data.frame(colData(bMat)) %>% 
  tibble::rownames_to_column(var = "cellname") %>% 
  mutate(cellType_broad = "Unknown")
row.names(cdat) = cdat$cellname 
cdat = select(cdat, Sample, cellType_broad, Responder, 
              Stimulator, Drug, Dose, Dose_Unit, Replicate,
              Relative_Dose, cellname, cm3 = clusters_mon3) 


# Create CDS from peak Matrix (SummarizedExperiment)
cds_b = monocle3::new_cell_data_set(assays(bMat)$TileMatrix, 
                                    cell_metadata = cdat,
                                    gene_metadata = rd)

# get bMat from old data
bMat_o = getMatrixFromProject( ArchRProj = prj_o, 
                             useMatrix = "TileMatrix", 
                             binarize = TRUE)

# format rowData 
rd_o= data.frame(rowData(bMat_o)) %>% 
  dplyr::mutate(stop = start + 499,
                bin = paste(seqnames, start, stop, sep = "_")) %>% 
  dplyr::select(bin)
row.names(rd_o) <- rd_o$bin
row.names(bMat_o) <- rd_o$bin

# format colData 
cdat_o = data.frame(colData(bMat_o)) %>% 
  tibble::rownames_to_column(var = "cellname") %>% 
  mutate(Drug = "NA",
         Dose = "NA",
         Dose_Unit = "NA",
         Relative_Dose = "NA", 
         cm3 = "NA")
row.names(cdat_o) = cdat_o$cellname 
cdat_o = select(cdat_o, Sample, cellType_broad, Responder, 
              Stimulator, Drug, Dose, Dose_Unit, Replicate,
              Relative_Dose, cellname, cm3) 

# Create CDS from peak Matrix (SummarizedExperiment)
cds_b_o = monocle3::new_cell_data_set(assays(bMat_o)$TileMatrix, 
                                    cell_metadata = cdat_o,
                                    gene_metadata = rd_o)

##########################################################
# remove drugged cells from newer data 
cds_bf = cds_b[,colData(cds_b)$Relative_Dose == 0]
cdat_f = data.frame(colData(cds_bf)) %>% 
  select(-Size_Factor)

# Manually combine the two cds objects

nd = counts(cds_bf)
od = counts(cds_b_o)
combined_counts = cbind(nd, od)
cdat_combined = rbind(cdat_f, cdat_o)

cds_combined = monocle3::new_cell_data_set(combined_counts, 
                                           cell_metadata = cdat_combined,
                                           gene_metadata = rd)

#########################################################

# plot both datasets together 
set.seed(2017)
cds_pl <- detect_genes(cds_combined)
ncells = ncol(cds_combined)*0.002
cds_pl = cds_pl[rowData(cds_pl)$num_cells_expressed > ncells,]
cds_pl <- estimate_size_factors(cds_pl)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
cds_pl = align_cds(cds_pl, preprocess_method = "LSI", alignment_group = "Sample")
reducedDim(cds_pl) <- reducedDim(cds_pl)[,2:50]
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "Aligned")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# prepare relevant custom UMAPs 
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))

# by monocle cluster
pdf("Monocle3_CombinedData_VehicleOnly_originalClusters_newData_aligned.pdf", width = 2, height = 1.75)
TCdat_pr %>% filter(cellType_broad == "Unknown") %>% 
ggplot( aes(x = UMAP_1, y = UMAP_2, color = cm3)) +
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


###########################
# use coembeded vehicle treated cells with previously annotated MLR data  
# to establish broad cell annotations for  
# downstream analysis. 

cd = data.frame(prj@cellColData) %>% 
  tibble::rownames_to_column(var = "cellNames") %>% 
  mutate(cellType_broad = 
           recode(clusters_mon3, 
                  "1" = "Dead", 
                  "2" = "Naive_Tcell",
                  "3" = "Monocyte",
                  "4" = "Naive_Tcell",
                  "5" = "NK",
                  "6" = "Naive_Tcell", 
                  "7" = "NK",
                  "8" = "Monocyte",
                  "9" = "Activated_Tcell", 
                  "10" = "Bcell"))


group_by(cd, cellType_broad) %>% 
  summarise(count = n())
  
#cd_o = data.frame(prj_o@cellColData) %>% 
#  group_by(cellType_broad) %>% 
#  summarise(count = n())
  
prj <- addCellColData(
  ArchRProj = prj, 
  data = cd$cellType_broad,
  cells = cd$cellNames, 
  name = "cellType_broad",
  force = TRUE)

# Save annotation info to project
saveArchRProject(ArchRProj= prj,
                 load = FALSE)


############################################
############################################
# start from saved project with info. 
prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered"), 
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
pdf("UMAP_cellTypes_with_Legend.pdf", width = 2.5, height = 1.75)
ggplot(cdat, aes(x = UMAP1, y = UMAP2, color = cellType_broad)) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  #theme(legend.position = "none") +
  #scale_color_brewer(palette = "Dark2") +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()


# filter out dead cells and make new project: 
prj_celltypes = prj[prj$cellType_broad != "Dead",]

saveArchRProject(ArchRProj= prj_celltypes,
                 outputDirectory = paste0(out_dir, "MLR_drug_filtered_annotated"),
                 dropCells = TRUE,
                 load = FALSE)







