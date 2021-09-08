basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB7"))
setwd(paste0(out_dir, "results/NB7"))

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
  library(viridis) 
  library(ggpubr)
})

set.seed(2017) # ensures reproducibility of previous random number generation

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

# target cells 
cts = c("Naive_Tcell", "Activated_Tcell")

# isolate cells from desired celltypes 
cds = cds_atac[,colData(cds_atac)$cellType_broad %in% cts]
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

cds <- cluster_cells(cds,reduction_method = "LSI")
colData(cds)$Cluster = clusters(cds, reduction_method="LSI")

cds <- cluster_cells(cds,reduction_method = "UMAP")
colData(cds)$louvain_component = cds@clusters[["UMAP"]]$partitions

# find trajectory
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
  dplyr::count(cellType_broad) %>%
  dplyr::mutate(total_cells = sum(n), fraction_cellType = n / total_cells)

root_nodes =  
  ordering_summary %>% 
  filter(cellType_broad == "Naive_Tcell" & fraction_cellType > 0.5)

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
saveRDS(cds, file = "cds_NaiveActive_Tcells_PT")

#######################
#######################

#start from processed cds
cds = readRDS(file = "cds_NaiveActive_Tcells_PT")

# plot UMAPs with trajectories
cat("plotting UMAP by cellType ...\n") 
colData(cds) %>%
  as.data.frame() %>%
  ggplot()+
  geom_point(aes( x = UMAP1, y =  UMAP2, color = cellType_broad), stroke = 0, size = 0.5)  +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
  #scale_color_manual("cell type") +
  #theme_void() +
  #theme(legend.position = "right", text = element_text(size = 6)) +
  ggsave("Naive_to_Activated_Tcells_umap.png",
         width = 2.5, height = 1.5, unit = "in", dpi = 900)

################
cat("plotting UMAP by pseudotime bin ...\n") 

pdf("Naive_to_Activated_Tcells_pseudotime_umap2.pdf", width = 1.75, height = 1.75)
colData(cds) %>%
  as.data.frame() %>%
  ggplot(aes( x = UMAP1, y =  UMAP2, color = Pseudotime))+
  geom_point_rast(size=0.4, stroke = 0) +
  monocle3:::monocle_theme_opts() +
  theme_void() +
  #theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
  theme(legend.position = "none") +
  scale_color_viridis(option = "inferno") 
dev.off() 
   

###############################################
cat( "plotting umap with root nodes emphasized...\n")
ica_space_df <- 
  t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
  as.data.frame() %>% dplyr::select(prin_graph_dim_1 = V1,
                                    prin_graph_dim_2 = V2) %>% 
  dplyr::mutate(sample_name = rownames(.),sample_state = rownames(.))

edge_df <- 
  cds@principal_graph[["UMAP"]] %>% 
  igraph::as_data_frame() %>% 
  dplyr::select(source = "from",target = "to") %>% 
  dplyr::left_join(ica_space_df %>% 
                     dplyr::select(source = sample_name,
                                   source_prin_graph_dim_1 = prin_graph_dim_1,
                                   source_prin_graph_dim_2 = prin_graph_dim_2),
                   by = "source") %>% 
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(target = sample_name, 
                                   target_prin_graph_dim_1 = prin_graph_dim_1,
                                   target_prin_graph_dim_2 = prin_graph_dim_2),
                   by = "target")

pdf("Naive_to_Activated_Tcells_pseudotime_UMAP_root_nodes2.pdf", width = 1.75, height = 1.75)
  ggplot() +
  geom_point_rast(data = colData(cds) %>% 
                  as.data.frame(), 
                  aes( x = UMAP1, y =  UMAP2,  color = Pseudotime),
                  size=0.4, stroke = 0) +
  monocle3:::monocle_theme_opts() +
  #theme(legend.position = "right", text = element_text(size = 8), 
  #      legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
  geom_segment(data = edge_df,
               aes_string(x="source_prin_graph_dim_1", 
                          y="source_prin_graph_dim_2", 
                          xend="target_prin_graph_dim_1", 
                          yend="target_prin_graph_dim_2"), 
               size=.7, linetype="solid", na.rm=TRUE) +
  geom_point(data = colData(cds) %>% 
               as.data.frame() %>%
               filter(root_node) %>%
               dplyr::select(source = closest_vertex) %>%
               mutate(source = paste("Y_", source, sep = "")) %>%
               left_join(edge_df, by = "source"), 
             aes(x = source_prin_graph_dim_1, y = source_prin_graph_dim_2), 
             color = "red", size = 1, stroke = 0) +
  theme_void() +
  scale_color_viridis(option = "inferno") +
  theme(legend.position = "none") 
dev.off() 

###############################################
cat("Plotting dose to pseudo bin comparisons (assess heterogeneity post stimulation)...\n")

cdat = as.data.frame(colData(cds)) %>% 
  mutate(responder_PT_bin = paste0(Responder, "_", Pseudotime_bin), 
         Stimulator_factor = factor(Stimulator, levels = c("noStim", "stimA", "stimB", "stimC", "stimD", "Bead")))

# get number of cells per PT bin per responder
Pseudotime_bin_summary = 
  cdat %>% 
  group_by(responder_PT_bin) %>%
  summarise(num_in_bin = n())

Pseudotime_bin_summary =
  left_join(cdat,
            Pseudotime_bin_summary, 
            by = "responder_PT_bin")

# get fraction of responder cells from each PT bin from a specific Stimulator 
Pseudotime_bin_summary = 
  Pseudotime_bin_summary %>%
  group_by(responder_PT_bin, Stimulator) %>%
  add_tally() %>%
  mutate(fraction_Stimulator = n/num_in_bin) %>%
  dplyr::select(cell, Responder, Stimulator, Stimulator_factor,
                responder_PT_bin, Pseudotime_bin, num_in_bin, fraction_Stimulator) %>%
  ungroup() %>% 
  filter(Responder != "stimAlone") # remove stim alone samples from this analysis. 

# stacked bars showing proportion of cells in each pseudodose based on Stimulator
# faceted by responder 
Pseudotime_bin_summary %>%
  dplyr::select(Pseudotime_bin, Responder, Stimulator_factor, fraction_Stimulator ) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = Pseudotime_bin, y = fraction_Stimulator, fill = Stimulator_factor ), 
           color = "black", size = .25, stat = "identity") +
  facet_wrap(~Responder, ncol = 2) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Pseudotime Bin") +
  ylab("Proportion") +
  ggsave("stackedbarplot_pseudotime_proportions.pdf",
         width = 3, height = 2.25, unit = "in", dpi = 1200)

###############################################
cat( "plotting pseudotime trajectory ridges for each stimulus...\n")
combined.df = cdat %>% 
  dplyr::select(Pseudotime, Pseudotime_bin, cellType_broad, Responder, Stimulator_factor) %>% 
  filter(Responder != "stimAlone")

combined.df %>%
  filter(is.finite(Pseudotime)) %>%
  ggplot() +
  geom_density_ridges(aes( x = Pseudotime, y = Stimulator_factor, fill = Stimulator_factor), size = .15) + 
  monocle3:::monocle_theme_opts() +
  facet_wrap(~Responder,ncol = 2) +
  theme(legend.position = "none", 
        text = element_text(size = 6),
        #axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
        #strip.text.x = element_blank()) + 
  xlab("Pseudotime") +
  ylab("Frequency per stimulation") + 
  ggsave(paste0("PseudotimeRidges_Tcells_byResponder.pdf"), height = 2 , width = 2, unit = "in")


###############################################
# plot Frags/cell by pseudotime bins 

ggplot(cdat) +
  geom_boxplot(aes(x = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)),
                   y = log10(nFrags), 
                   fill = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)))) +
  theme_bw() +
  xlab("Pseudotime bin") +
  ylab("log10(nFrags)") +
  scale_fill_viridis_d(option ='inferno') +
  theme(legend.position = "none") +
ggsave(paste0("BoxPlot_nFrags_Tcells_byPT2.pdf"),
       height = 2 , width = 2, unit = "in")

######################################

if(assay == "ATAC"){
  cat("Create a aggregate cds (colums = number of bins)\n")
  colData(cds)$cell_subtype <- kmeans(colData(cds)$Pseudodose, centers=10)$cluster
  table(colData(cds)$cell_subtype)
  
  grouping <- plyr::ddply(as.data.frame(colData(cds)), plyr::.(cell_subtype), function(x) {
    cells <- x$cell
    if(length(cells) < 100) {
      return(data.frame(cell = x$cell, 
                        num_genes_expressed = x$num_genes_expressed,
                        Pseudodose = x$Pseudodose, 
                        cell_bin = paste(x$cell_subtype, 1, sep="_")))
    } else {
      num_cuts <- floor(length(cells)/50)
      return(data.frame(cell=x$cell, 
                        num_genes_expressed = x$num_genes_expressed,
                        Pseudodose = x$Pseudodose,
                        cell_bin = paste(x$cell_subtype, ntile(x$Pseudodose, n = num_cuts), sep="_")))
    }
  })
  
  row.names(grouping) <- grouping$cell
  grouping <- grouping[row.names(colData(cds)),]
  colData(cds)$cell_bin <- grouping$cell_bin
  
  table(colData(cds)$cell_bin)
  
  # Aggregates cells by pseudo dose bin
  binned_cds <-aggregate_by_cell_bin(cds, "cell_bin")
  
  group_info <- plyr::ddply(as.data.frame(colData(cds)), plyr::.(cell_bin), function(x) {
    data.frame(mean_pseudo = mean(x$Pseudodose), 
               mean_num_genes = mean(x$num_genes_expressed),
               ncells = nrow(x))
  })
  row.names(group_info) <- group_info$cell_bin
  group_info$group_name <- row.names(group_info)
  
  colData(binned_cds) <- merge(colData(binned_cds), group_info)
  
  row.names(colData(binned_cds)) <- colData(binned_cds)$cell_bin
  colData(binned_cds) <- colData(binned_cds)[row.names(group_info),]
  
  rowData(binned_cds) <- rowData(cds)[row.names(binned_cds),]
  
  binned_cds <- detect_genes(binned_cds, min_expr=0.1)
  binned_cds <- estimate_size_factors(binned_cds)
  
  rowData(binned_cds)$use_for_ordering <- FALSE
  
  cat("saving agrregate CDS with cells binned by pseudodose_bins")
  saveRDS(binned_cds, file = paste0(drug, "/", assay, "_", drug,"_agg_pseudodose_cds"))
}


# function to return models and coefficients for aggregated ATAC data by pseudodose (from binned_cds generated with "analyze_drug_psuedodose")
get_binned_DA_sites_spline = function(drug, outdir = paste0(out_dir, "results/NB7/")){
  cds=readRDS(file = paste0(outdir, drug, "/ATAC_", drug, "_agg_pseudodose_cds"))
  cat("filter CDS features (only those found in >20 cells)...\n")
  idx_r = rowData(cds)$num_cells_expressed > 20
  cds_f = cds[idx_r,]
  cat("modeling accessibility by mean_pseudodose (per bulked bin). This will take a few minutes ...\n")
  # here counts are estimated to be distributed according to negative binomial b/c counts are no longer binary (cells are binned)
  gene_fits = fit_models(cds_f, expression_family = "negbinomial", 
                         model_formula_str = "~ splines::ns(mean_pseudo, df = 3) + splines::ns(mean_num_genes, df = 3)")
  gene_fits$site_name = row.names(rowData(cds_f))
  fit_coefs = coefficient_table(gene_fits)
  dose_terms = fit_coefs %>% dplyr::filter(grepl("mean_pseudo", term))
  cat("writing dose terms to file: ", drug,"_binned_spline_terms_by_site", "...\n")
  saveRDS(gene_fits, file  = paste0(outdir, drug, "/ATAC_", drug, "_binned_spline_models_by_site"))
  saveRDS(dose_terms, file = paste0(outdir, drug, "/ATAC_", drug, "_binned_spline_terms_by_site"))
}


