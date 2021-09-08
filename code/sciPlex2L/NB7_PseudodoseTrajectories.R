basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190521_scichem2_AllPlates/"
bin_directory = paste0(basepath, "analysis/bin/")
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB7"))
setwd(paste0(out_dir, "results/NB7/"))

suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(ggridges)
  library(tidymodels)
  library(devtools)
  library(monocle3)
  library(cicero)
  library(furrr)
  library(viridis)
  library(piano)
  library(UpSetR)
  library(snowfall)
  plan(multicore)
})

set.seed(1) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

suppressPackageStartupMessages({
  source(paste0(bin_directory,
                "cell_cycle.R",
                sep = ""))
  source(paste0(bin_directory,
                "dose_response.R",
                sep = ""))
  source(paste0(bin_directory,
                "viability.R",
                sep = ""))
  cc.genes <- readRDS(paste0(bin_directory,
                             "cc.genes.RDS",
                             sep = ""))
  source(paste0(bin_directory,
                "viability.R",
                sep = ""))
  source(paste0(bin_directory,
                "dispersions_functions.R",
                sep = ""))
  source(paste0(bin_directory,
                "GSA_helper_functions.R",
                sep = ""))
  source(paste0(bin_directory,
                "loadGSCSafe.R",
                sep = ""))
  
})

# load ArchR project 
prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered"),
                       showLogo = FALSE)


###########################################
# check available matrices in archr project 
getAvailableMatrices(prj)

# load PeakMatrix into memory
# load PeakMatrix into memory
pMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"sc2_filtered/ArrowFiles/sc2.arrow"),
                          ArchRProj = prj, 
                          cellNames = prj$cellNames,
                          useMatrix = "PeakMatrix", 
                          binarize = TRUE)

# Create CDS from peak Matrix (SummarizedExperiment)
rd = data.frame(rowRanges(pMat)) %>% 
  rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end))

# peaks in this list are in a DIFFERENT ORDER 
# relative to pMat (rowRanges above is correct)
pd = data.frame(getPeakSet(prj)) %>% 
  rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end), 
         gene_short_name = ifelse(peakType == "Promoter" & distToGeneStart < 500, nearestGene, "NA")) %>% 
  dplyr::select(rowname, score, replicateScoreQuantile, groupScoreQuantile, 
                Reproducibility, GroupReplicate, distToGeneStart,
                nearestGene, peakType, distToTSS, nearestTSS,
                gene_short_name, GC, idx, peak)

a_rd = left_join(rd, pd, by = "peak") %>% 
  select(rowname = rowname.x, idx = idx.x, score, replicateScoreQuantile, groupScoreQuantile, 
         Reproducibility, GroupReplicate, distToGeneStart,
         nearestGene, peakType, distToTSS, nearestTSS,
         gene_short_name, GC, peak)

row.names(a_rd) <- a_rd$peak
row.names(pMat) <- a_rd$peak
cds_atac = monocle3::new_cell_data_set(assays(pMat)$PeakMatrix, 
                                       cell_metadata = colData(pMat), 
                                       gene_metadata = a_rd)

#################################################
# run analysis on ATAC cds
sciPlex_cds <- cds_atac

dim(sciPlex_cds)

colData(sciPlex_cds)$cell =
  row.names(colData(sciPlex_cds))

colData(sciPlex_cds)$new_treatment_label =
  sapply(colData(sciPlex_cds)$treatment, function(x) {
    if (grepl("BMS", x))
      return("BMS345541")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
  })

colData(sciPlex_cds)$solvent =
  sapply(pData(sciPlex_cds)$treatment, function(x) {
    if (grepl("Dex", x))
      return("Ethanol")
    else
      return("DMSO")
  })

cellDat = data.frame(colData(sciPlex_cds))

######################################################
######################################################
# function to iterate analysis over each drug. 
# outputs drug specific results to it's own directory. 

analyze_drug_psuedodose = function(cds_in = sciPlex_cds, assay = "ATAC", 
                                   drug, dose_levels, clr_values, Seed = 2017){ 
  set.seed(Seed)
  dir.create(paste0(drug, "/"))
  # isolate cells from target drug 
  cds = cds_in[,colData(cds_in)$new_treatment_label == drug]
  colData(cds)$dose_character = factor(colData(cds)$dose, levels = dose_levels)
  cds = detect_genes(cds)
  cds = estimate_size_factors(cds)
  cds <- estimate_cell_cycle(cds,
                             g1s_markers = cc.genes$s.genes,
                             g2m_markers = cc.genes$g2m.genes)
  
  cat("preprocessing ", assay, "  CDS of ", drug, "-treated cells...\n")
  if(assay == "RNA"){
    cds <- preprocess_cds(cds,
                          method = 'PCA',
                          num_dim = 25,
                          norm_method = 'log',
                          verbose = T)
  
    cds = reduce_dimension(cds,
                           max_components = 2,
                           preprocess_method = "PCA",
                           reduction_method = 'UMAP',
                           umap.metric = 'cosine',
                           umap.n_neighbors = 10,
                           umap.min_dist = 0.1,
                           verbose = TRUE)
  
    cds <- cluster_cells(cds,reduction_method = "PCA")
    colData(cds)$Cluster = clusters(cds, reduction_method="PCA")
  }
  else{
    cds <- preprocess_cds(cds,
                          method = 'LSI',
                          num_dim = 30,
                          #norm_method = 'log',
                          verbose = T)
    
    cds = reduce_dimension(cds,
                           max_components = 2,
                           preprocess_method = "LSI",
                           reduction_method = 'UMAP',
                           umap.metric = 'cosine',
                           umap.n_neighbors = 10,
                           umap.min_dist = 0.1,
                           verbose = TRUE)
    
    cds <- cluster_cells(cds,reduction_method = "LSI")
    colData(cds)$Cluster = clusters(cds, reduction_method="LSI")
  }
  
  cds <- cluster_cells(cds,reduction_method = "UMAP")
  colData(cds)$louvain_component = cds@clusters[["UMAP"]]$partitions
  
  graph_parameters = list()
  graph_parameters[["minimal_branch_len"]] = 10
  graph_parameters[["ncenter"]] = 750
  
  cds <- learn_graph(cds, 
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
    dplyr::count(dose) %>%
    dplyr::mutate(total_cells = sum(n), fraction_dose = n / total_cells)
  
  root_nodes =  
    ordering_summary %>% 
    filter(dose == 0 & fraction_dose > 0.5)
  
  root_nodes = root_nodes$closest_vertex
  
  colData(cds)$root_node = colData(cds)$closest_vertex %in% root_nodes
  
  root_cells  = 
    colData(cds) %>% 
    as.data.frame() %>% 
    filter(root_node) %>% 
    pull(cell) %>% 
    as.character()
  
  cds <- order_cells(cds, root_cells=root_cells)
  
  colData(cds)$Pseudodose = cds@principal_graph_aux[["UMAP"]]$pseudotime
  
  #colData(cds)$Pseudodose_bin = cut(colData(cds)$Pseudodose, breaks=quantile(colData(cds)$Pseudodose, seq(0, 1, 0.1)), labels=F)
  # edit: rather than just binning, partition cells into pseudodose groups using clustering. 
 
  colData(cds)$Pseudodose_bin <- kmeans(colData(cds)$Pseudodose, centers=10)$cluster
  # rename Pseudodose_bin to reflect ordering of average pseudodose
  x = as.data.frame(colData(cds)) %>% group_by(Pseudodose_bin) %>% 
    mutate(mean_pd = mean(Pseudodose)) %>% ungroup() %>% 
    mutate(Pseudodose_bin_ordered = dense_rank((mean_pd))) %>% 
    dplyr::select(treatment, Pseudodose, Pseudodose_bin, Pseudodose_bin_ordered)
  colData(cds)$Pseudodose_bin = x$Pseudodose_bin_ordered
  
  cat("number of cells in each pseudodose bin = \n")
  print(table(colData(cds)$Pseudodose_bin))
  
  #####################
  cat(" plotting proliferation index as a function of dose...\n")
  colData(cds)$drug_dose = paste0(colData(cds)$new_treatment_label, "_", colData(cds)$dose)
  
  colData(cds)$proliferation_cutoff  = 1.664 #median value for my data
  
  coldata_cds = 
    colData(cds) %>%
    as.data.frame() %>%
    group_by(cell_type) %>%
    mutate(low_proliferation = proliferation_index < proliferation_cutoff) %>%
    ungroup() %>%
    group_by(Pseudodose_bin,drug_dose) %>%
    add_tally() %>%
    mutate(fraction_in_low = sum(low_proliferation)/n)
  
  coldata_cds %>% 
    filter(!is.na(Pseudodose_bin)) %>%
    dplyr::select(treatment,Pseudodose_bin,fraction_in_low) %>%
    distinct() %>%
    ggplot(aes(x = factor(Pseudodose_bin), y = 100*(fraction_in_low))) +
    geom_boxplot(fill = "grey80", outlier.stroke = 0, outlier.size = 0.5, outlier.colour =  "grey80") +
    geom_jitter(size = 0.5, stroke = 0) +
    monocle3:::monocle_theme_opts() +
    theme(text = element_text(size = 6)) +
    xlab("Pseudodose Bin") +
    ylab("Percent Low") +
    ggsave(paste0(drug, "/", assay, "_", drug, "_percent_cells_in_low_proliferation_index.pdf"), height = 1.5, width = 2.5)
  
  ################
  cat("plotting UMAP by dose ...\n") 
  colData(cds) %>%
    as.data.frame() %>%
    ggplot()+
    geom_point(aes( x = UMAP1, y =  UMAP2, color = dose_character), stroke = 0, size = 0.5)  +
    monocle3:::monocle_theme_opts() +
    theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
    scale_color_manual("Dose (uM)", values = clr_values) +
    #theme_void() +
    #theme(legend.position = "right", text = element_text(size = 6)) +
    ggsave(paste0(drug, "/", assay, "_", drug,"_dose_umap.png"),
           width = 2.5, height = 1.5, unit = "in", dpi = 900)
  
  ################
  cat("plotting UMAP by pseudodose bin ...\n") 
  colData(cds) %>%
    as.data.frame() %>%
    ggplot()+
    geom_point(aes( x = UMAP1, y =  UMAP2, color = Pseudodose_bin), stroke = 0, size = 0.5)  +
    monocle3:::monocle_theme_opts() +
    theme(legend.position = "right", text = element_text(size = 8), legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
    #theme_void() +
    #theme(legend.position = "right", text = element_text(size = 6)) +
    ggsave(paste0(drug, "/", assay, "_",drug,"_pseudodose_umap.png"),
           width = 2.5, height = 1.5, unit = "in", dpi = 900)
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
  
  ggplot() +
    geom_point(data = 
                 colData(cds) %>%
                 as.data.frame(),
               aes(x = UMAP1, y = UMAP2, color = dose_character), 
               size = .5, stroke = 0)+
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
               color = "red", size = 1, stroke = 0)+ 
    monocle3:::monocle_theme_opts() +
    theme(legend.position = "right", text = element_text(size = 8), 
          legend.key.width = unit(0.4,"line"), legend.key.height = unit(0.5,"line")) + 
    scale_color_manual("Dose (uM)", values = clr_values) +
    #theme_void()+
    #theme(legend.position = "none", text = element_text(size = 6)) +
    ggsave(paste0(drug, "/", assay, "_",drug, "_dose_UMAP_root_nodes.png"),
           width = 2.5, height = 1.5, unit = "in", dpi = 900)  
  
  ###############################################
  cat("Plotting dose to pseudo bin comparisons (assess heterogeneity at each dose)...\n")
  
  Pseudodose_bin_summary = 
    pData(cds) %>% 
    as.data.frame() %>%
    group_by(Pseudodose_bin) %>%
    summarise(num_in_bin = n())
  
  Pseudodose_bin_summary =
    left_join(as.data.frame(pData(cds)),
              Pseudodose_bin_summary, 
              by = "Pseudodose_bin")
  
  Pseudodose_bin_summary = 
    Pseudodose_bin_summary %>%
    group_by(Pseudodose_bin, cell_type) %>%
    add_tally() %>%
    mutate(fraction_cell_type = n/num_in_bin) %>%
    dplyr::select(-n) %>%
    ungroup()
  
  # not critical for me since I only have one HDACi (SAHA)
  Pseudodose_bin_summary = 
    Pseudodose_bin_summary %>%
    group_by(Pseudodose_bin, treatment) %>% 
    add_tally() %>%
    mutate(fraction_product_name = n/num_in_bin) %>%
    dplyr::select(-n) %>%
    ungroup()
  
  Pseudodose_bin_summary = 
    Pseudodose_bin_summary %>%
    group_by(Pseudodose_bin, dose_character) %>%
    add_tally() %>%
    mutate(fraction_dose = n/num_in_bin) %>%
    dplyr::select(-n) %>%
    ungroup()
  
  # stacked bars showing proportion of cells in each pseudodose bin
  Pseudodose_bin_summary %>%
    dplyr::select(Pseudodose_bin,dose_character, fraction_dose ) %>%
    distinct() %>%
    ggplot() +
    geom_bar(aes(x = Pseudodose_bin, y =fraction_dose, fill = dose_character ), 
             color = "black", size = .25, stat = "identity") +
    scale_fill_manual("Dose (uM)", values = clr_values) +
    monocle3:::monocle_theme_opts() +
    theme(legend.position = "right", 
          text = element_text(size = 6),        
          legend.key.width = unit(0.15,"line"), 
          legend.key.height = unit(0.1,"line")) +    
    scale_x_continuous(breaks = seq(1,10,1))+
    xlab("Pseudodose Bin") +
    ylab("Proportion") +
    ggsave(paste0(drug, "/", assay, "_",drug, "_stackedbarplot_pseudodose_proportions.pdf"),
           width = 2, height = 1.5, unit = "in", dpi = 1200)
  
  ###############################################
  cat( "plotting psuedodose trajectory ridges for each real dose...\n")
  combined.df = data.frame(colData(cds)) %>% dplyr::select(cell_type, new_treatment_label,  Pseudodose, dose_character)
  
  combined.df %>%
    filter(is.finite(Pseudodose)) %>%
    filter(new_treatment_label == drug) %>%
    ggplot() +
    geom_density_ridges(aes( x = Pseudodose, y = dose_character, fill = dose_character), size = .15) + 
    monocle3:::monocle_theme_opts() +
    facet_wrap(~cell_type,ncol = 3) +
    scale_fill_manual("Dose (uM)", values = clr_values) + 
    theme(legend.position = "none", 
          text = element_text(size = 6),
          #axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text.x = element_blank()) + 
    ggtitle(drug) + 
    xlab("Pseudodose") +
    ylab("Frequency per Dose (uM)") + 
    ggsave(paste0(drug, "/", assay, "_",drug, "_Pseudodose_ridges.pdf"), height = 1.5 , width = 1.5, unit = "in")
  
  cat("saving processed CDS with pseudodose_bins\n")
  saveRDS(cds, file = paste0(drug, "/", assay, "_", drug, "_cds_PD"))
  
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

######################################################
######################################################
# use function to analyze each drug separately 
#SAHA
analyze_drug_psuedodose(cds_in = sciPlex_cds, assay = "ATAC", drug = "SAHA", Seed = 1,
                        dose_levels = c("0", "1", "5", "10", "50", "100", "500", "1000"),
                        clr_values = c("0"="gray", "1"="#1D1147FF", "5"="#51127CFF", "10"="#822681FF", "50"="#B63679FF",
                                       "100"="#E65164FF", "500" = "#FB8861FF", "1000"="#FEC287FF"))

get_binned_DA_sites_spline(drug = "SAHA")

#Dex
analyze_drug_psuedodose(cds_in = sciPlex_cds, assay = "ATAC", drug = "Dex", Seed = 1,
                        dose_levels = c("0", "0.05", "0.25", "0.5", "2.5", "5", "25", "50"),
                        clr_values = c("0"="gray", "0.05"="#1D1147FF", "0.25"="#51127CFF", "0.5"="#822681FF", "2.5"="#B63679FF",
                                       "5"="#E65164FF", "25" = "#FB8861FF", "50"="#FEC287FF"))

get_binned_DA_sites_spline(drug = "Dex")
##BMS
analyze_drug_psuedodose(cds_in = sciPlex_cds, assay = "ATAC", drug = "BMS345541", Seed = 1,
                        dose_levels = c("0", "1", "5", "10", "50", "100", "500", "1000"),
                        clr_values = c("0"="gray", "1"="#1D1147FF", "5"="#51127CFF", "10"="#822681FF", "50"="#B63679FF",
                                       "100"="#E65164FF", "500" = "#FB8861FF", "1000"="#FEC287FF"))

get_binned_DA_sites_spline(drug = "BMS345541")
#Nutlin
analyze_drug_psuedodose(cds_in = sciPlex_cds, assay = "ATAC", drug = "Nutlin3A", Seed = 1,
                        dose_levels = c("0", "0.25", "1.25", "2.5", "12.5", "25", "125", "250"),
                        clr_values = c("0"="gray", "0.25"="#1D1147FF", "1.25"="#51127CFF", "2.5"="#822681FF", 
                                       "12.5"="#B63679FF", "25"="#E65164FF", "125" = "#FB8861FF", "250"="#FEC287FF"))

get_binned_DA_sites_spline(drug = "Nutlin3A")
##########################
# run parallel analysis on sciRNA data 
# note: integrated matrix is not raw counts, and 
# multiple ATAC cells assigned to same RNA cell. 

# load integratedGeneMatrix into memory
#rnaMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"ArrowFiles/sc2.arrow"),
#                            ArchRProj = prj, 
#                            useMatrix = "GeneIntegrationMatrix", 
#                            binarize = FALSE)
# Create CDS from gene Integration Matrix (SummarizedExperiment)
# MATCHED sciRNA-seq cells
#r_rd = data.frame(rowData(rnaMat))
#r_rd$gene_short_name <- r_rd$name
#row.names(r_rd) <- r_rd$name
#row.names(rnaMat) <- r_rd$name
#cds_rna_integrated = monocle3::new_cell_data_set(assays(rnaMat)$GeneIntegrationMatrix, 
#                                                 cell_metadata = colData(pMat),
#                                                 gene_metadata = r_rd)

# load raw sciPlex-RNA cds. 
cds.rna <- readRDS(paste0(basepath, "analysis/archr_revised/scRNA/sciPlex2_cds_NB5processed.RDS"))
recode_key = c(BMS_0 = "0", BMS_0.1 = "1", BMS_0.5 = "5", BMS_1 = "10", 
               BMS_5 = "50", BMS_10 = "100", BMS_50 = "500", BMS_100 = "1000",
               Dex_0 = "0", Dex_0.1 = "0.5", Dex_0.5 = "2.5", Dex_1 = "5",
               Dex_5 ="25", Dex_10 = "50", Dex_50 = "250", Dex_100 = "500",
               Nutlin_0 = "0", Nutlin_0.1 = "0.25", Nutlin_0.5 = "1.25",Nutlin_1 = "2.5",
               Nutlin_5 = "12.5",Nutlin_10 = "25", Nutlin_50 = "125", Nutlin_100 = "250",
               SAHA_0 = "0", SAHA_0.1 = "1", SAHA_0.5 = "5", SAHA_1 = "10",
               SAHA_5 = "50",SAHA_10 = "100",SAHA_50 = "500", SAHA_100 = "1000")

colData(cds.rna)$dose = recode(colData(cds.rna)$treatment_RNAbased, !!!recode_key)


# use function to analyze sciRNA each drug separately 
#SAHA
analyze_drug_psuedodose(cds_in = cds.rna, assay = "RNA", drug = "SAHA", Seed = 1,
                        dose_levels = c("0", "1", "5", "10", "50", "100", "500", "1000"),
                        clr_values = c("0"="gray", "1"="#1D1147FF", "5"="#51127CFF", "10"="#822681FF", "50"="#B63679FF",
                                       "100"="#E65164FF", "500" = "#FB8861FF", "1000"="#FEC287FF")
)

#Dex
analyze_drug_psuedodose(cds_in = cds.rna, assay = "RNA", drug = "Dex", Seed = 1,
                        dose_levels = c("0", "0.5", "2.5", "5", "25", "50", "250", "500"),
                        clr_values = c("0"="gray", "0.5"="#1D1147FF", "2.5"="#51127CFF", "5"="#822681FF", "25"="#B63679FF",
                                       "50"="#E65164FF", "250" = "#FB8861FF", "500"="#FEC287FF")
)

##BMS
analyze_drug_psuedodose(cds_in = cds.rna, assay = "RNA", drug = "BMS345541", Seed = 1,
                        dose_levels = c("0", "1", "5", "10", "50", "100", "500", "1000"),
                        clr_values = c("0"="gray", "1"="#1D1147FF", "5"="#51127CFF", "10"="#822681FF", "50"="#B63679FF",
                                       "100"="#E65164FF", "500" = "#FB8861FF", "1000"="#FEC287FF")
)

#Nutlin
analyze_drug_psuedodose(cds_in = cds.rna, assay = "RNA", drug = "Nutlin3A", Seed = 1,
                        dose_levels = c("0", "0.25", "1.25", "2.5", "12.5", "25", "125", "250"),
                        clr_values = c("0"="gray", "0.25"="#1D1147FF", "1.25"="#51127CFF", "2.5"="#822681FF", 
                                       "12.5"="#B63679FF", "25"="#E65164FF", "125" = "#FB8861FF", "250"="#FEC287FF")
)

















