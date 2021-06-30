basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB10"))
setwd(paste0(out_dir, "results/NB10"))

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(monocle3)
  library(cicero)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
  library(Gviz)
  library(biomaRt)
  library(RColorBrewer)
  library(viridis)
})

set.seed(2017)

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

# load processed cds (with pseudotime info)
cds = readRDS(file = paste0(out_dir,"results/NB7/cds_NaiveActive_Tcells_PT"))

#broadly define allo and bead stimulations
cdat = colData(cds) %>% 
  data.frame() %>% 
  mutate(stim_type = 
           ifelse(Responder == "stimAlone", "stimAlone", 
                  ifelse(Stimulator == "Bead", "bead", 
                         ifelse(Stimulator == "noStim", "noStim",
                                ifelse(stringr::str_sub(Stimulator, -1) == stringr::str_sub(Responder, -1), "auto", "allo")))))

colData(cds)$stim_type = cdat$stim_type


# count Activated Tcells from each broad stim_type
# including no stim etc. 

T_act_counts = filter(cdat, cellType_broad == "Activated_Tcell") %>% 
  group_by(stim_type) %>% 
  summarise(Tcells_act = n()) 

TcellCounts = group_by(cdat, stim_type) %>% 
  summarise(Tcells_total = n()) %>% 
  left_join(T_act_counts, by = "stim_type")

############
# Isolate only the bead or allo-Activated Tcells 
idx_cells = 
  (colData(cds)$cellType_broad == "Activated_Tcell" & 
  (colData(cds)$stim_type == "allo" | colData(cds)$stim_type == "bead"))

cds_a = cds[,idx_cells]

# Pseudotime ridges plot with ONLY Activated Tcells
combined.df = colData(cds_a) %>% 
  data.frame() %>% 
  mutate(Stimulator_factor = factor(Stimulator, levels = c("stimA", "stimB", "stimC", "stimD", "Bead"))) %>% 
  dplyr::select(Pseudotime, Pseudotime_bin, cellType_broad, Responder, Stimulator, Stimulator_factor , stim_type)

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
  ggsave(paste0("PseudotimeRidges_ONLY_ActivatedTcells_byResponder.pdf"), height = 2 , width = 2, unit = "in")

# Pseudotime ridges plot comparing allo (combined) vs beads 
combined.df %>%
  filter(is.finite(Pseudotime)) %>% 
  ggplot() +
  geom_density_ridges(aes( x = Pseudotime, y = stim_type, fill = stim_type), size = .15) + 
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
  ggsave(paste0("PseudotimeRidges_ONLY_ActivatedTcells_byResponder_alloVSbead.pdf"), height = 2 , width = 2, unit = "in")

# box plots comparing pseudotimes positions of allo and bead stimulated Tcells 
combined.df %>%
  filter(is.finite(Pseudotime)) %>% 
  ggplot(aes(x = stim_type, y = Pseudotime)) +
  geom_violin(aes(fill = stim_type)) +
  geom_boxplot(width=0.1) +
  facet_wrap(~Responder,ncol = 2) +
  scale_fill_brewer(palette='Set1') +
  coord_flip() +
  theme_bw()+ 
  xlab("Stimulation") +
  ylab("Pseudotime Position") +
  theme(legend.position = "none") +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0("BoxPlot_PT_byStimType2.pdf"),
       height = 3 , width = 3, unit = "in")


####################
####################

# remake UMAP only using activated Tcells
set.seed(2017)
cds_pl <- detect_genes(cds_a)
cds_pl <- estimate_size_factors(cds_pl)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
reducedDim(cds_pl) <- reducedDim(cds_pl)[,2:50]
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# prepare relevant custom UMAPs 
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))

# plot UMAP (by stim type)
pdf("ActivatedTcells_UMAP_StimType2.pdf", width = 2.5, height = 1.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = stim_type)) +
  geom_point_rast(size= 1.25, stroke = 0) + 
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
       axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  scale_color_brewer(palette ="Set1")
dev.off()

# plot UMAP (by pseudotime)
pdf("ActivatedTcells_UMAP_Pseudotime2.pdf", width = 1.75, height = 1.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
  geom_point_rast(size= 1.25, stroke = 0) + 
  theme_void() +
  theme(legend.position = "none") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  scale_color_viridis(option ="inferno")
dev.off()

#######################
#######################

# perform DA analysis on activated T-cells 2 ways
# 1) Sites influenced by Bead vs Allo stim (activated Tcells only).
# 2) Sites influenced by position in pseudotime (activated Tcells only). 

# Only need to run fit_models once: load the saved output in next section 

# subset to peaks open in > 0.5% cells
cds_a <- detect_genes(cds_a)
cutoff = ncol(cds_a)*.005
cds_f = cds_a[rowData(cds_a)$num_cells_expressed > cutoff, ]

# 1)
# fit regression models for each peak across with stimulation type as predictor 
peak_fits_stimType <- fit_models(
  cds = cds_f, 
  model_formula_str = "~stim_type + nFrags", 
  expression_family = "binomial")

# get coeficients from fit_models
peak_fit_coefs_ST = coefficient_table(peak_fits_stimType)

peak_fit_coefs_ST_f = dplyr::select(peak_fit_coefs_ST, 
                                 peak, num_cells_expressed, status, term,
                                 estimate, std_err, test_val, p_value, 
                                 normalized_effect, model_component, q_value)
write.table(peak_fit_coefs_ST_f, 
            file =  "Tcell_StimType_DApeak_coefs.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")

# 2)
# fit regression models for each peak across with Pseudotime as predictor 
peak_fits_PT <- fit_models(
  cds = cds_f, 
  model_formula_str = "~Pseudotime + nFrags", 
  expression_family = "binomial")

# get coeficients from fit_models
peak_fit_coefs_PT = coefficient_table(peak_fits_PT)

peak_fit_coefs_PT_f = dplyr::select(peak_fit_coefs_PT, 
                                 peak, num_cells_expressed, status, term,
                                 estimate, std_err, test_val, p_value, 
                                 normalized_effect, model_component, q_value)
write.table(peak_fit_coefs_PT_f, 
            file =  "Tcell_Pseudotime_DApeak_coefs.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")

##########################
##########################

# functions for looking for element enrichment within DA peaks (taken from NB6.5)
pull_peak_by_feature_mat = function(archr_project, file = "Motif-Matches-In-Peaks.rds"){
  Ann = readRDS(paste0(getOutputDirectory(archr_project), "/Annotations/", file))
  pf_mat = assay(Ann)
  peaks = data.frame(rowRanges(Ann)) %>% # GB Corrected 11/17/2020
    dplyr::mutate(peak = paste(seqnames, start, end, sep = "_"))
  row.names(pf_mat) <- peaks$peak
  return(pf_mat)
}

# Function:
# fit a linear regression to determine if features predict binary peak behavior.  
# each peak is an observation, features are predictors and the response variable is binary classification 
# for each peak (i.e. sig opening vs not). This gives a coef estimate for each feature (i.e. motif)
get.feature.logistic.reg.stats = function(response, features, data) {
  t(sapply(features, function(feature) {
    res = summary(glm(paste(response, feature, sep=" ~ "), family = "binomial", data = data))
    if (nrow(res$coefficients) > 1) {
      return(c(res$coefficients[2,1], res$coefficients[2,4]))
    } else {
      return(c(NA, NA))
    }
  }))
}

# function tests weather DA peaks are enriched for features 
feature_enrichment = function(
  Peak_coefs, 
  feature_mat,
  pVal_cutoff = 0.05)
{
  # Deafine peaks as DA or unchanged for cell type 
  peak_coefs = 
    Peak_coefs %>%
    dplyr::mutate(direction = ifelse(p_value <= pVal_cutoff & estimate > 0, "Opening",
                                     ifelse(p_value <= pVal_cutoff & estimate < 0, "Closing",
                                            "Unchanged")))
  
  # filter peak x feature matrix for only used peaks (in DA analysis)
  feature_mat_f = feature_mat[peak_coefs$peak,]
  
  # list of feature names
  features = colnames(feature_mat_f)
  
  # replace any characters that will screw up formulas
  features = stringr::str_replace_all(features, "[+-.*]", "_")
  
  colnames(feature_mat_f) <- features
  
  # reformat from logical to numeric
  fmf = matrix(as.numeric(feature_mat_f), nrow = nrow(x=feature_mat_f))
  row.names(fmf) = row.names(feature_mat_f)
  colnames(fmf) = colnames(feature_mat_f)
  feature_df = as.data.frame(fmf)
  feature_df$peak = row.names(feature_df)
  # add direction info from coef_table
  feature_df = dplyr::inner_join(peak_coefs, feature_df, by = "peak")
  # add binary columns (dummy variables) describing directionality of each site 
  feature_df = cbind(feature_df, as.data.frame(model.matrix(~ 0 + direction, feature_df)))
  
  cat("fit linear regression predicting DA direction from motif presence. \nThis will take a few minutes...\n")
  
  # set container
  peak.DA.stats = list()
  
  # test for feature association with OPENING peaks
  if(dplyr::filter(feature_df, direction == "Opening") %>% 
     nrow > 0){
    peak.DA.stats$opening = get.feature.logistic.reg.stats(
      "directionOpening", features, data = feature_df)
    
    peak.DA.stats$opening = data.frame(
      feature = row.names(peak.DA.stats$opening),
      beta = peak.DA.stats$opening[,1],
      p.val = peak.DA.stats$opening[,2]) %>%
      mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
      dplyr::filter(!is.na(Padj), beta > 0) %>% 
      dplyr::arrange(Padj) %>% 
      dplyr::mutate(rank = rank(Padj), 
                    sig = ifelse(Padj <0.05, 1, 0), 
                    mlog10Padj = -log10(Padj)) 
  }
  
  # test for feature association with CLOSING peaks
  if(dplyr::filter(feature_df, direction == "Closing") %>% 
     nrow > 0){
    peak.DA.stats$closing = get.feature.logistic.reg.stats(
      "directionClosing", features, data = feature_df)
    
    peak.DA.stats$closing = data.frame(
      feature = row.names(peak.DA.stats$closing),
      beta = peak.DA.stats$closing[,1],
      p.val = peak.DA.stats$closing[,2]) %>%
      mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
      dplyr::filter(!is.na(Padj), beta > 0) %>% 
      dplyr::arrange(Padj) %>% 
      dplyr::mutate(rank = rank(Padj), 
                    sig = ifelse(Padj <0.05, 1, 0), 
                    mlog10Padj = -log10(Padj)) 
  }
  
  return(peak.DA.stats) 
}

# funciton for plotting  motifs ranked by Pvals
plot_motif_enrichment = function(motif_test_results, fname = "ranked_sig_motifs.pdf")
{
  ggUp <- ggplot(motif_test_results, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point_rast(size = 2) +
    ggrepel::geom_label_repel(
      data = motif_test_results[1:5, ],
      aes(x = rank, y = mlog10Padj, label = feature), 
      size = 3,
      nudge_x = 2,
      color = "black"
    ) + 
    theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet"))
  
  ggsave(plot = ggUp, 
         filename = fname,
         width = 3.75, 
         height = 3.75)
}

######################
######################

######################
# test for feature enrichments in DA sites
# using above functions.
######################

# load ArchR project
prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 

# load peak x motif matrix
motif_mat = pull_peak_by_feature_mat(
  archr_project = prj, 
  file = "Motif-Matches-In-Peaks.rds")

# pull peak x encode TFBS matrix 
Encode_mat = pull_peak_by_feature_mat(
  archr_project = prj, 
  file = "EncodeTFBS-Matches-In-Peaks.rds")
# correct_formatting of encode colnames
cn = paste0(stringr::str_split_fixed(colnames(Encode_mat), pattern = "\\.", n = 3)[,2],
            "_",
            stringr::str_split_fixed(colnames(Encode_mat), pattern = "\\.", n = 3)[,1])
colnames(Encode_mat)<- cn

# load DA coefficients comparing Allo vs. Bead Stim. T-cells. 
# DA based on Stimulation type (ST)
# NOTE:  Here positive coefs mean open in Bead relative to Allo
coefs_ST = read.table("Tcell_StimType_DApeak_coefs.txt", head = TRUE) %>% 
  filter(term == "stim_typebead") 

# run Motif enrichment 
motif_enrichment_ST = feature_enrichment(
  Peak_coefs = coefs_ST, 
  feature_mat = motif_mat,
  pVal_cutoff = 0.05)

# run Encode enrichment 
Encode_enrichment_ST = feature_enrichment(
  Peak_coefs = coefs_ST, 
  feature_mat = Encode_mat,
  pVal_cutoff = 0.05)

#plot Opening sig.Features by rank  
plot_motif_enrichment(motif_enrichment_ST$opening, fname = "ActivatedTcell_BeadOpening_DAsites_ranked_motifs2.pdf")
plot_motif_enrichment(Encode_enrichment_ST$opening, fname = "ActivatedTcell_BeadOpening_DAsites_ranked_EncodeTFBS2.pdf")

#plot Closing sig. TF motifs by rank  
plot_motif_enrichment(motif_enrichment_ST$closing, fname = "ActivatedTcell_BeadClosing_DAsites_ranked_motifs2.pdf")
plot_motif_enrichment(Encode_enrichment_ST$closing, fname = "ActivatedTcell_BeadClosing_DAsites_ranked_EncodeTFBS2.pdf")


############
# repeat analysis with DA sites based on pseudotime
# DA based on Pseudotime (PT) %>% 
coefs_PT = read.table("Tcell_Pseudotime_DApeak_coefs.txt", head = TRUE) %>% 
  filter(term == "Pseudotime")

motif_enrichment_PT = feature_enrichment(
  Peak_coefs = coefs_PT, 
  feature_mat = motif_mat,
  pVal_cutoff = 0.05)

# run Encode enrichment 
Encode_enrichment_PT = feature_enrichment(
  Peak_coefs = coefs_PT, 
  feature_mat = Encode_mat,
  pVal_cutoff = 0.05)

#plot Opening sig. TF motifs by rank  
plot_motif_enrichment(motif_enrichment_PT$opening, fname = "ActivatedTcell_PseudotimeOpening_DAsites_ranked_motifs2.pdf")
plot_motif_enrichment(Encode_enrichment_PT$opening, fname = "ActivatedTcell_PseudotimeOpening_DAsites_ranked_EncodeTFBS2.pdf")

#plot Closing sig. TF motifs by rank  
plot_motif_enrichment(motif_enrichment_PT$closing, fname = "ActivatedTcell_PseudotimeClosing_DAsites_ranked_motifs2.pdf")
plot_motif_enrichment(Encode_enrichment_PT$closing, fname = "ActivatedTcell_PseudotimeClosing_DAsites_ranked_EncodeTFBS2.pdf")

##### 
# make a barplot of "DA sites" between bead and allo stim. 
DAsites_ST = filter(coefs_ST, p_value <= 0.05) %>%
  mutate(direction = ifelse(estimate < 0, "allo_open", "allo_closed")) %>% # Note: that coefs relative to allo condition (i.e. opposite signs)
  group_by(direction) %>% 
  summarise(nSites = n())

ggplot(DAsites_ST) +
  geom_col(aes(x = direction, y = nSites, fill = direction)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "barplot_DAsites_nominal_AlloVsBead.pdf",
       width = 2.5, 
       height = 2.5)
  
# make a barplot of "DA sites" across Pseudotime (Activated Tcells only). 
DAsites_PT = filter(coefs_PT, p_value <= 0.05) %>%
  mutate(direction = ifelse(estimate > 0, "PT_open", "PT_closed")) %>% # Note: that coefs relative to allo condition (i.e. opposite signs)
  group_by(direction) %>% 
  summarise(nSites = n())

ggplot(DAsites_PT) +
  geom_col(aes(x = direction, y = nSites, fill = direction))+ 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "barplot_DAsites_nominal_Pseudotime.pdf",
       width = 2.5, 
       height = 2.5)  

##############
##############
# What genes are closest / connected to Allo specific peaks?

# pull in peak annotations from ArchR project. 
peak_dat = getPeakSet(prj) %>% 
  data.frame() %>% 
  mutate(peak = paste(seqnames, start, end, sep = "_")) %>% 
  dplyr::select(peak, distToGeneStart, nearestGene, peakType, distToTSS, nearestTSS, GC)

# add info to coef tables 
coefs_ST = left_join(coefs_ST, peak_dat, by = "peak")
coefs_PT = left_join(coefs_PT, peak_dat, by = "peak")

# create output for inputting into online tools (GREAT etc.)
# sites by Bead-Allo comparison
nearbyGenes_AlloUp_ST = filter(coefs_ST, p_value <= 0.05, estimate < 0) %>% 
  dplyr::select(nearestGene)
backgroundGenes_ST = dplyr::select(coefs_ST, nearestGene)
write.table(nearbyGenes_AlloUp_ST, file = "nearbyGenes_AlloUp_ST.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(backgroundGenes_ST, file = "backgroundGenes_ST.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# sites by Pseudotime Term
nearbyGenes_AlloUp_PT = filter(coefs_PT, p_value <= 0.05, estimate > 0) %>% 
  dplyr::select(nearestGene)
backgroundGenes_PT = dplyr::select(coefs_PT, nearestGene)

write.table(nearbyGenes_AlloUp_PT, file = "nearbyGenes_AlloUp_PT.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
write.table(backgroundGenes_PT, file = "backgroundGenes_PT.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# create bed files from the sig and background peaks (for input into GREAT) 
# sites by Bead-Allo comparison
AlloUp_sites_ST = filter(coefs_ST, p_value <= 0.05, estimate < 0) %>% 
  mutate(chr = stringr::str_split_fixed(peak, "_", 3)[,1], 
         start = stringr::str_split_fixed(peak, "_", 3)[,2], 
         end = stringr::str_split_fixed(peak, "_", 3)[,3]) %>% 
  dplyr::select(chr, start, end, peak)

background_sites_ST = 
  mutate(coefs_ST, 
         chr = stringr::str_split_fixed(peak, "_", 3)[,1], 
         start = stringr::str_split_fixed(peak, "_", 3)[,2], 
         end = stringr::str_split_fixed(peak, "_", 3)[,3]) %>% 
  dplyr::select(chr, start, end, peak)

write.table(AlloUp_sites_ST, file = "AlloUp_sites_ST.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(background_sites_ST, file = "background_sites_ST.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# sites by Pseudotime Term
AlloUp_sites_PT = filter(coefs_PT, p_value <= 0.05, estimate > 0) %>% 
  mutate(chr = stringr::str_split_fixed(peak, "_", 3)[,1], 
         start = stringr::str_split_fixed(peak, "_", 3)[,2], 
         end = stringr::str_split_fixed(peak, "_", 3)[,3]) %>% 
  dplyr::select(chr, start, end, peak)

background_sites_PT = 
  mutate(coefs_PT, 
         chr = stringr::str_split_fixed(peak, "_", 3)[,1], 
         start = stringr::str_split_fixed(peak, "_", 3)[,2], 
         end = stringr::str_split_fixed(peak, "_", 3)[,3]) %>% 
  dplyr::select(chr, start, end, peak)

write.table(AlloUp_sites_PT, file = "AlloUp_sites_PT.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
write.table(background_sites_PT, file = "background_sites_PT.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)



