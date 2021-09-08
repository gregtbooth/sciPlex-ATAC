# Note: I didn't pursue this analysis further, because I thin the two 
# clusters of "monocytes" represent disticnct cell types and not 
# allogeneic stimulation dependendent chromatin states. 

basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB11"))
setwd(paste0(out_dir, "results/NB11"))

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

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 

###################
###################

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

###################
###################

set.seed(2017)

# simplify stimulation labels
cdat = colData(cds_atac) %>% 
  data.frame() %>% 
  mutate(stim_type = 
           ifelse(Responder == "stimAlone", "stimAlone", 
                  ifelse(Stimulator == "Bead", "bead", 
                         ifelse(Stimulator == "noStim", "noStim",
                                ifelse(stringr::str_sub(Stimulator, -1) == stringr::str_sub(Responder, -1), "noStim", "allo")))))

colData(cds_atac)$stim_type = cdat$stim_type

# target cells 
cts = c("Monocyte")
stims = c("noStim", "allo")

# isolate monoctyes and desired stimulation
cds = cds_atac[,colData(cds_atac)$cellType_broad %in% cts & colData(cds_atac)$stim_type %in% stims]
colData(cds)$cell = row.names(colData(cds))

cds <- preprocess_cds(
  cds,
  method = 'LSI',
  num_dim = 50,
  verbose = T)

reducedDim(cds) <- reducedDim(cds)[,2:50]

cds = reduce_dimension(
  cds,
  max_components = 2,
  preprocess_method = "LSI",
  reduction_method = 'UMAP',
  umap.metric = 'cosine',
  umap.n_neighbors = 10,
  umap.min_dist = 0.1,
  verbose = TRUE)

cds <- cluster_cells(cds,reduction_method = "LSI", k = 20, resolution=1e-2)
colData(cds)$Cluster = clusters(cds, reduction_method="LSI")

cds <- cluster_cells(cds,reduction_method = "UMAP")
colData(cds)$louvain_component = cds@clusters[["UMAP"]]$partitions

graph_parameters = list()
graph_parameters[["minimal_branch_len"]] = 10
graph_parameters[["ncenter"]] = 750

cds <- learn_graph(
  cds, 
  learn_graph_control= graph_parameters,
  use_partition = FALSE, 
  close_loop = FALSE)

colData(cds)$UMAP1 = reducedDims(cds)[["UMAP"]][,1]
colData(cds)$UMAP2 = reducedDims(cds)[["UMAP"]][,2]

# Get the closest vertex for every cell
colData(cds)$closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex[,1]

ordering_summary =
  colData(cds) %>%
  as.data.frame() %>%
  dplyr::group_by(closest_vertex) %>%
  dplyr::count(stim_type) %>%
  dplyr::mutate(total_cells = sum(n), fraction_cellType = n / total_cells)

root_nodes =  
  ordering_summary %>% 
  filter(stim_type == "noStim" & fraction_cellType > 0.5)

root_nodes = root_nodes$closest_vertex

colData(cds)$root_node = colData(cds)$closest_vertex %in% root_nodes

root_cells  = 
  colData(cds) %>% 
  as.data.frame() %>% 
  filter(root_node) %>% 
  pull(cell) %>% 
  as.character()

cds <- order_cells(cds, root_cells=root_cells)

colData(cds)$Pseudotime = cds@principal_graph_aux[["UMAP"]]$pseudotime

# edit: rather than just binning, partition cells into Pseudotime groups using clustering. 

colData(cds)$Pseudotime_bin <- kmeans(colData(cds)$Pseudotime, centers=10)$cluster

# rename Pseudodose_bin to reflect ordering of average Pseudotime
x = as.data.frame(colData(cds)) %>% group_by(Pseudotime_bin) %>% 
  mutate(mean_pt = mean(Pseudotime)) %>% ungroup() %>% 
  mutate(Pseudotime_bin_ordered = dense_rank((mean_pt))) %>% 
  dplyr::select(cellType_broad, Pseudotime, Pseudotime_bin, Pseudotime_bin_ordered)
colData(cds)$Pseudotime_bin = x$Pseudotime_bin_ordered

cat("number of cells in each Pseudotime bin = \n")
print(table(colData(cds)$Pseudotime_bin))

cat("saving processed CDS with pseudodose_bins\n")
saveRDS(cds, file = "cds_monocytes_PT")

###################
###################

cds = readRDS(file = "cds_monocytes_PT")

# plot UMAPs with trajectories
cat("plotting UMAP by cluster ...\n") 
colData(cds) %>%
  as.data.frame() %>%
  ggplot()+
  geom_point(aes( x = UMAP1, y =  UMAP2, color = Cluster), stroke = 0, size = 0.5)  +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
  #scale_color_manual("cell type") +
  #theme_void() +
  #theme(legend.position = "right", text = element_text(size = 6)) +
  ggsave("UMAP_Monocyte_clusters.png",
         width = 2.5, height = 1.5, unit = "in", dpi = 900)

# plot UMAP by stim_type
cat("plotting UMAP by Stimulation type ...\n") 
colData(cds) %>%
  as.data.frame() %>%
  ggplot()+
  geom_point(aes( x = UMAP1, y =  UMAP2, color = stim_type), stroke = 0, size = 0.5)  +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
  ggsave("UMAP_Monocyte_stimType.png",
         width = 2.5, height = 1.5, unit = "in", dpi = 900)

# plot UMAP by Pseudotime bin
cat("plotting UMAP by Pseudotime bin ...\n") 
colData(cds) %>%
  as.data.frame() %>%
  ggplot()+
  geom_point(aes( x = UMAP1, y =  UMAP2, color = Pseudotime_bin), stroke = 0, size = 0.5)  +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
  ggsave("UMAP_Monocyte_Pseudotime_bin.png",
         width = 2.5, height = 1.5, unit = "in", dpi = 900)



  
  
  
  
  
  
  


