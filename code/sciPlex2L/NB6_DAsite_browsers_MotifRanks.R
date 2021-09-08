basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190521_scichem2_AllPlates/"
bin_directory = paste0(basepath, "analysis/bin/")
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB6"))
setwd(paste0(out_dir, "results/NB6"))


suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(cicero)
  library(ggplot2)
  library(ggridges)
  library(ggrastr)
  library(tidymodels)
  library(devtools)
  library(furrr)
  library(viridis)
  library(piano)
  library(UpSetR)
  library(snowfall)
  library(dplyr)
  plan(multicore)
})

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)

#########################################################

# load full ArchR project
prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered"),
                       showLogo = FALSE)

# create peak by motif matrix


# load DA site list: 
peak_model_res = read.table(file = paste0(out_dir, "results/NB5/sciPlex_ATAC_A549_small_screen_DApeaks.txt"), head = T)

################################################
# make drug specific plots 
################################################

dir.create("Plots")

for(d in c("Dex", "SAHA", "Nutlin", "BMS")){
  dir.create(d)

  prj_d = prj[prj@cellColData$treatment == d]
  prj_d$treatment_hashmod = paste0(prj_d$treatment, "_", prj_d$dose)

  dose_term = paste0("dose_", d)
  # isolate sig changed sites
  das = peak_model_res %>%
    dplyr::filter(grepl(dose_term, term)) %>%
    dplyr::filter(q_value <= 0.05) %>% 
    dplyr::arrange(-estimate) %>% 
    dplyr::mutate(seqnames = stringr::str_split_fixed(peak, "_", n = 3)[,1],
           start = stringr::str_split_fixed(peak, "_", n = 3)[,2], 
           end = stringr::str_split_fixed(peak, "_", n = 3)[,3]) %>% 
    dplyr::select(seqnames, start, end, peak, term, estimate, q_value) 

  #make GRanges object
  das_GR = makeGRangesFromDataFrame(das, keep.extra.columns=TRUE)

  # plot browser tracks of Dex treated cells by dose 
  TopUpPeaks = head(das_GR, 6)
  TopDownPeaks = tail(das_GR, 6)

  # extend windows for plotting
  TopUpPeaks_extended = TopUpPeaks + 25000
  TopDownPeaks_extended = TopDownPeaks + 25000

  p_up <- plotBrowserTrack(
    ArchRProj = prj_d,
    groupBy = "treatment_hashmod", 
    region = TopUpPeaks_extended,
    plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
    sizes = c(10, 1.5, 2),
    tileSize = 250
  )

  p_down <- plotBrowserTrack(
    ArchRProj = prj_d, 
    groupBy = "treatment_hashmod", 
    region = TopDownPeaks_extended,
    plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
    sizes = c(10, 1.5, 2),
    tileSize = 250
  )

  # plot up peaks
  plotPDF(plotList = p_up,
          name = paste0(d, "BrowserTrack_TopUpsites_new"), 
          width = 5, 
          height = 5, 
          ArchRProj = NULL, 
          addDOC = FALSE)

  plotPDF(plotList = p_down,
          name = paste0(d, "BrowserTrack_TopDownsites_new"), 
          width = 5, 
          height = 5, 
          ArchRProj = NULL, 
          addDOC = FALSE)
} 

###########################
###########################
# pull peak motif data for enrichment analysis. 
# function to retrieve and format peak by motif matrix
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

for(d in c("Dex", "SAHA", "Nutlin", "BMS")){
  dir.create(d)

  dose_term = paste0("dose_", d)
  
  # get all coefficients based on drug dose
  peak_coefs = peak_model_res %>%
    dplyr::filter(grepl(dose_term, term)) %>% 
    dplyr::mutate(direction = ifelse(q_value < 0.05 & estimate > 0, "Opening", 
                              ifelse(q_value < 0.05 & estimate < 0, "Closing", "Unchanged")))

  # get peak x motif matrix 
  motif_mat = pull_peak_by_motif_mat(prj)

  # filter for only used peaks
  motif_mat_f = motif_mat[peak_coefs$peak,]

  # list of motif names
  motifs = colnames(motif_mat_f)

  # reformat from logical to numeric
  mmf = matrix(as.numeric(motif_mat_f), nrow = nrow(x=motif_mat_f))
  row.names(mmf) = row.names(motif_mat_f)
  colnames(mmf) = colnames(motif_mat_f)
  motif_df = as.data.frame(mmf)
  motif_df$peak = row.names(motif_df)
  # add direction info from coef_table
  motif_df = dplyr::inner_join(peak_coefs, motif_df, by = "peak")
  # add binary columns describing directionality of each site 
  motif_df = cbind(motif_df, as.data.frame(model.matrix(~ 0 + direction, motif_df)))

  cat("fit linear regression predicting DA direction from motif presence. \nThis will take a few minutes...\n")

  # set container
  peak.DA.stats = list()

  if(dplyr::filter(motif_df, direction == "Opening") %>% 
     nrow > 0){
    peak.DA.stats$promoter.opening = get.peak.DA.logistic.reg.stats(
      "directionOpening", motifs, motif_df)
    #filter for only the significant motifs
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

  if(dplyr::filter(motif_df, direction == "Closing") %>% 
     nrow > 0){
    peak.DA.stats$promoter.closing = get.peak.DA.logistic.reg.stats(
      "directionClosing", motifs, motif_df)
    #filter for only the significant motifs
    peak.DA.stats$promoter.closing = data.frame(
      motif = row.names(peak.DA.stats$promoter.closing),
      beta = peak.DA.stats$promoter.closing[,1],
      p.val = peak.DA.stats$promoter.closing[,2]) %>%
      mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
      dplyr::filter(!is.na(Padj), beta > 0) %>% 
      dplyr::arrange(Padj) %>% 
      dplyr::mutate(rank = rank(Padj), 
             sig = ifelse(Padj <0.05, 1, 0), 
             mlog10Padj = -log10(Padj)) 
  }

  # plot motifs ranked by Pvals 
  ggUp <- ggplot(peak.DA.stats$promoter.opening, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point_rast(size = 1) +
    ggrepel::geom_label_repel(
      #data = peak.DA.stats$promoter.opening[peak.DA.stats$promoter.opening$sig ==1, ],
      data = peak.DA.stats$promoter.opening[1:5, ],
      aes(x = rank, y = mlog10Padj, label = motif), 
      size = 1.5,
      nudge_x = 2,
      color = "black"
    ) + 
    theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))


  ggDo <- ggplot(peak.DA.stats$promoter.closing, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point_rast(size = 1) +
    ggrepel::geom_label_repel(
      #data = peak.DA.stats$promoter.closing[peak.DA.stats$promoter.closing$sig ==1, ],
      data = peak.DA.stats$promoter.closing[1:5, ],
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
         filename = paste0(d, "/", d, "_motifs_openingSites_rast_Padj.pdf"),
         width = 3.75, 
         height = 3.75)

  ggsave(plot = ggDo, 
         filename = paste0(d, "/", d, "_motifs_closingSites_rast_Padj.pdf"),
         width = 3.75, 
         height = 3.75)
}

###############################################################
# Compare expression of TFs (from sciPlex-RNA) with sig. motifs 
###############################################################

# first process sciPlexRNA cds object
set.seed(2017)
cds.rna <- 
  readRDS(paste0(basepath, "analysis/archr_revised/scRNA/sciPlex2_cds_NB5processed.RDS"))

# process cds 
cds.rna <- preprocess_cds(cds.rna,
                          method = 'PCA',
                          num_dim = 25,
                          norm_method = 'log',
                          verbose = T)

cds.rna = reduce_dimension(cds.rna,
                           max_components = 2,
                           preprocess_method = "PCA",
                           reduction_method = 'UMAP',
                           umap.metric = 'cosine',
                           umap.n_neighbors = 10,
                           umap.min_dist = 0.1,
                           verbose = TRUE)

cds.rna <- cluster_cells(cds.rna,reduction_method = "PCA")
colData(cds.rna)$Cluster = clusters(cds.rna, reduction_method="PCA")

# manually selected list of TFs targetting top sig motifs for each drug
# SAHA
TFgenes_SAHA_up = c("LHX8", "CEBPZ", "NFYA", "NFYB", "NFYC", "PBX3")
TFgenes_SAHA_down = c("FOXA1", "FOXA2", "HNF1B", "FOXA3", "FOXO3")

#Dex
TFgenes_Dex_up = c("NR3C1", "NR2E3", "PROX1", "ZNF282", "TGIF1")
TFgenes_Dex_down = c("FOXA3", "FOXS1", "ATF7", "MESP2", "MESP1")

#Nutlin
TFgenes_Nutlin_up = c("TP63", "TP53", "MEF2C", "TP73", "HMGA1")
#TFgenes_Nutlin_down = c("ETV4", "ETV5", "ETV1", "OTX2", "ELK1")

#BMS
TFgenes_BMS_up = c("TFAP4", "MYF5", "PRDM4", "PAX9")
TFgenes_BMS_down = c("CEBPB", "CEBPE", "NKX25", "TEAD3", "AR")



d = "SAHA"
dose_term = paste0("dose_", d)

# Opening Motif TF expression
TF_genes_opening <-
  cds.rna[rowData(cds.rna)$gene_short_name %in% TFgenes_SAHA_up,
          colData(cds.rna)$treatment == d |
            colData(cds.rna)$vehicle == TRUE]

TF_genes_closing <-
  cds.rna[rowData(cds.rna)$gene_short_name %in% TFgenes_SAHA_down,
          colData(cds.rna)$treatment == d |
            colData(cds.rna)$vehicle == TRUE]


pdf(
  paste0(d, "/", d, "_OpeningMotif_TF_percent_RNApositive.pdf"),
  width = 2,
  height = 2
)
plot_percent_cells_positive(
  TF_genes_opening,
  group_cells_by = "dose_character",
  ncol = 2,
  plot_as_count = FALSE
) +
  scale_fill_manual(
    "Dose (µM)",
    values = c(
      "0" = "gray",
      "0.1" = "#1D1147FF",
      "0.5" = "#51127CFF",
      "1" = "#822681FF",
      "5" = "#B63679FF",
      "10" = "#E65164FF",
      "50" = "#FB8861FF",
      "100" = "#FEC287FF"
    )
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 4)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)")
dev.off()


# Closing Motif TF expression 

pdf(
  paste0(d, "/", d, "_ClosingMotif_TF_percent_RNApositive.pdf"),
  width = 2,
  height = 2
)
plot_percent_cells_positive(
  TF_genes_closing,
  group_cells_by = "dose_character",
  ncol = 2,
  plot_as_count = FALSE
) +
  scale_fill_manual(
    "Dose (µM)",
    values = c(
      "0" = "gray",
      "0.1" = "#1D1147FF",
      "0.5" = "#51127CFF",
      "1" = "#822681FF",
      "5" = "#B63679FF",
      "10" = "#E65164FF",
      "50" = "#FB8861FF",
      "100" = "#FEC287FF"
    )
  ) +
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 4)
  ) +
  xlab("Log(Dose [µM])") +
  ylab("Cells (percent)")
dev.off()






