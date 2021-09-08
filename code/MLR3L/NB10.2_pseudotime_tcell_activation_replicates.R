basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB10_2"))
setwd(paste0(out_dir, "results/NB10_2"))

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

###################
###################

# load cds
cds_in = readRDS(file = paste0(out_dir,"results/NB7/cds_NaiveActive_Tcells_PT"))

replicate_PT_analysis = function(CDS){
  rep_cds_list = list()
  
  # run PT analysis separately for each rep
  for (rep in 1:3){
    cds = CDS[,colData(CDS)$Replicate == rep]
    
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
    rep_cds_list[[rep]] = cds
  }
  return(rep_cds_list)
}


plot_umaps = function(cds, rep = "rep1"){
  # plot UMAPs with trajectories
  cat("plotting UMAP by cellType ...\n") 
  pdf(paste0(rep, "_Naive_to_Activated_Tcells_umap.pdf"), width = 2.25, height = 1.75)
  print(
    colData(cds) %>%
      as.data.frame() %>%
      ggplot(aes( x = UMAP1, y =  UMAP2, color = cellType_broad))+
      geom_point_rast(size=0.4, stroke = 0) +
      monocle3:::monocle_theme_opts() +
      theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) 
  )
  dev.off()   

  cat("plotting UMAP by pseudotime bin ...\n") 

  pdf(paste0(rep, "_Naive_to_Activated_Tcells_pseudotime_umap2.pdf"), width = 1.75, height = 1.75)
  print(
    colData(cds) %>%
      as.data.frame() %>%
      ggplot(aes( x = UMAP1, y =  UMAP2, color = Pseudotime))+
      geom_point_rast(size=0.4, stroke = 0) +
      monocle3:::monocle_theme_opts() +
      theme_void() +
      theme(legend.position = "none") +
      scale_color_viridis(option = "inferno") 
  )
  dev.off() 
}

plot_frags_over_PT = function(cds, rep = "rep1"){
  pdf(paste0(rep, "_BoxPlot_nFrags_Tcells_byPT2.pdf"), width = 2, height = 2)
  print(
    colData(cds) %>%
      as.data.frame() %>%
      ggplot() +
      geom_boxplot(aes(x = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)),
                       y = log10(nFrags), 
                       fill = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)))) +
      theme_bw() +
      xlab("Pseudotime bin") +
      ylab("log10(nFrags)") +
      scale_fill_viridis_d(option ='inferno') +
      theme(legend.position = "none") 
  )
  dev.off() 
}


plot_PT_ridges = function(cds, rep = "rep1"){
  cdat = as.data.frame(colData(cds)) %>% 
    mutate(responder_PT_bin = paste0(Responder, "_", Pseudotime_bin), 
           Stimulator_factor = factor(Stimulator, levels = c("noStim", "stimA", "stimB", "stimC", "stimD", "Bead")))
  
  combined.df = cdat %>% 
    dplyr::select(Pseudotime, Pseudotime_bin, cellType_broad, Responder, Stimulator_factor) %>% 
    filter(Responder != "stimAlone")
  
  pdf(paste0(rep, "_PseudotimeRidges_Tcells_byResponder.pdf"), width = 2, height = 2)
  print(
    combined.df %>%
      filter(is.finite(Pseudotime)) %>%
      ggplot(aes( x = Pseudotime, y = Stimulator_factor, fill = Stimulator_factor)) +
      geom_density_ridges( size = .15) + 
      monocle3:::monocle_theme_opts() +
      facet_wrap(~Responder,ncol = 2) +
      theme(legend.position = "none", 
            text = element_text(size = 6),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Pseudotime") +
      ylab("Frequency per stimulation") 
  )
  dev.off()
}


# get fraction of responder cells from each PT bin from a specific Stimulator 
plot_PT_stacked_bar = function(cds, rep = "rep1"){
  cdat = as.data.frame(colData(cds)) %>% 
    mutate(responder_PT_bin = paste0(Responder, "_", Pseudotime_bin), 
           Stimulator_factor = factor(Stimulator, levels = c("noStim", "stimA", "stimB", "stimC", "stimD", "Bead")))
  
  Pseudotime_bin_summary = 
    cdat %>% 
    group_by(responder_PT_bin) %>%
    summarise(num_in_bin = n())
  
  Pseudotime_bin_summary =
    left_join(cdat,
              Pseudotime_bin_summary, 
              by = "responder_PT_bin")
  
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
  pdf(paste0(rep, "_stackedbarplot_pseudotime_proportions.pdf.pdf"), width = 3, height = 2.25)
  print(
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
      ylab("Proportion")
  )
  dev.off()
}


#################
rep_cds_list_processed = replicate_PT_analysis(CDS=cds_in)

# plot separate umaps
plot_umaps(cds = rep_cds_list_processed[[1]], rep = "rep1")
plot_umaps(cds = rep_cds_list_processed[[2]], rep = "rep2")
plot_umaps(cds = rep_cds_list_processed[[3]], rep = "rep3")

# plot separate umis per cell across PT bins
plot_frags_over_PT(cds = rep_cds_list_processed[[1]], rep = "rep1")
plot_frags_over_PT(cds = rep_cds_list_processed[[2]], rep = "rep2")
plot_frags_over_PT(cds = rep_cds_list_processed[[3]], rep = "rep3")

# plot separate ridge plots across PT bins
plot_PT_ridges(cds = rep_cds_list_processed[[1]], rep = "rep1")
plot_PT_ridges(cds = rep_cds_list_processed[[2]], rep = "rep2")
plot_PT_ridges(cds = rep_cds_list_processed[[3]], rep = "rep3")

# plot separate stacked  barplots across PT bins
plot_PT_stacked_bar(cds = rep_cds_list_processed[[1]], rep = "rep1")
plot_PT_stacked_bar(cds = rep_cds_list_processed[[2]], rep = "rep2")
plot_PT_stacked_bar(cds = rep_cds_list_processed[[3]], rep = "rep3")



