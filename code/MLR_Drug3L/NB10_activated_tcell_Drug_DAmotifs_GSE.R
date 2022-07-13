basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/211103_3Level_scichem_MLR_Drugs/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB10"))
setwd(paste0(out_dir, "results/NB10"))

suppressPackageStartupMessages({
  library(ArchR)
  #library(Seurat)
  library(monocle3)
  library(dplyr)
  library(tidymodels)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
  library(RColorBrewer)
  #library(ggpubr)
  library(viridis)
  library(snowfall)
  library(furrr)
  library(gt)
  library(piano)
  library(pheatmap)
  plan(multicore)
})

suppressPackageStartupMessages({
  source(paste0(out_dir, "external_resources/bin/GSA_helper_functions.R"))
  source(paste0(out_dir, "external_resources/bin/loadGSCSafe.R"))
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

#########################
#########################

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
      dplyr::mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
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
      dplyr::mutate(Padj = p.adjust(p.val, method = "BH")) %>% 
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
# test for feature enrichments in DA sites
# using above functions.
######################

# load ArchR project
prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered_annotated"), 
                       showLogo = FALSE) 

# Add annotation matrices to project 
# Motifs
#prj <- addMotifAnnotations(ArchRProj = prj, motifSet = "cisbp", name = "Motif")
# Encode TFBS 
#prj <- addArchRAnnotations(ArchRProj = prj, collection = "EncodeTFBS")


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

# load Dose-dependent DA coefficients from Drug-Treated Allo-Activated T-cells. 
# NOTE:  Here positive coefs mean openning in allo-activated T-cells with increasing drug
# load model output from NB9
# DA sites
Peak_coefs_Drug = read.table(
  file =  paste0(out_dir, "results/NB9/sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DApeaks.txt"), 
  head = TRUE) %>% 
  dplyr::filter(grepl("dose", term))

# 
drug_terms =  c(dplyr::distinct(Peak_coefs_Drug, term)) 
for(drug in drug_terms$term){
  cat("looking for motif and TFBS envrichment within DA sites affected by", drug, "\n")

  coefs_drug_f = dplyr::filter(Peak_coefs_Drug, grepl(drug, term))
  # run Motif enrichment 
  motif_enrichment_ST = feature_enrichment(
    Peak_coefs = coefs_drug_f, 
    feature_mat = motif_mat,
    pVal_cutoff = 0.05)
  
  write.table(motif_enrichment_ST$opening, file = paste0("AlloActivatedTcell_",drug, "_DAsites_opening_motifEnrichment.txt"),
                                         sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  write.table(motif_enrichment_ST$closing, file = paste0("AlloActivatedTcell_",drug, "_DAsites_closing_motifEnrichment.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  # run Encode enrichment 
  Encode_enrichment_ST = feature_enrichment(
    Peak_coefs = coefs_drug_f, 
    feature_mat = Encode_mat,
    pVal_cutoff = 0.05)
  
  write.table(Encode_enrichment_ST$opening, file = paste0("AlloActivatedTcell_",drug, "_DAsites_opening_ENCODEEnrichment.txt"),
                                         sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  write.table(Encode_enrichment_ST$closing, file = paste0("AlloActivatedTcell_",drug, "_DAsites_closing_ENCODEEnrichment.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  #plot Opening sig.Features by rank  
  if(!is.null(motif_enrichment_ST$opening)){
    print(
    plot_motif_enrichment(motif_enrichment_ST$opening, fname = paste0("AlloActivatedTcell_",drug, "_opening_DAsites_ranked_motifs.pdf"))
    )
    }
  if(!is.null(Encode_enrichment_ST$opening)){
    print(
    plot_motif_enrichment(Encode_enrichment_ST$opening, fname =  paste0("AlloActivatedTcell_",drug, "_opening_DAsites_ranked_EncodeTFBS.pdf"))
    )
    }
  #plot Closing sig. TF motifs by rank  
  if(!is.null(motif_enrichment_ST$closing)){
    print(
    plot_motif_enrichment(motif_enrichment_ST$closing, fname =  paste0("AlloActivatedTcell_",drug, "_closing_DAsites_ranked_motifs.pdf"))
    )  
    }
  if(!is.null(Encode_enrichment_ST$closing)){
    print(
    plot_motif_enrichment(Encode_enrichment_ST$closing, fname =  paste0("AlloActivatedTcell_",drug, "_closing_DAsites_ranked_EncodeTFBS.pdf"))
    )
    }
}

#############################################
# DA genes 
DAGs = read.table( 
  file =  paste0(out_dir, "results/NB9/sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DAgenes.txt"), 
  head = TRUE)



##########################
########## GSEA ##########
##########################

# load DA coefficients comparing Allo vs. Bead Stim. T-cells. 
# DA based on Stimulation type (ST)
# NOTE:  Here positive coefs mean open in Bead relative to Allo
  
coefs_ST = read.table( 
  file =  paste0(out_dir, "results/NB9/sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DAgenes.txt"), 
  head = TRUE)

## Load Gene Set Collection
hallmarks <-loadGSCSafe(file= paste0(out_dir,"external_resources/GeneSets/h.all.v6.0.symbols.gmt"))

# This function is designed to collect GSA stats on DAgenes 
get_gsa_res <- function(DAgene_test_res, gsc, min_genes_in_set=5, ncores=1){
  DAgene_test_res = DAgene_test_res %>% 
    dplyr::group_by(gene_short_name) %>%
    dplyr::filter(abs(test_val) == max(abs(test_val))) 
  directions=data.frame(row.names=as.character(DAgene_test_res$gene_short_name), stats=DAgene_test_res$test_val)
  geneLevelStats=data.frame(row.names=as.character(DAgene_test_res$gene_short_name), stats=DAgene_test_res$test_val)
  tryCatch({
    gsares <- runGSA(geneLevelStats=geneLevelStats, direction=geneLevelStats, geneSetStat="mean", gsc=gsc, nPerm=10000, ncpu=ncores, gsSizeLim=c(min_genes_in_set, Inf), verbose=T)
    return (gsares)
  }, error = function(e) { print(e); return (NA) } )
}
Sys.setenv(OMP_NUM_THREADS = 1)


# run GSEA using DA genes (comparing Allo vs Bead directly)
drug_terms =  c(distinct(coefs_ST, term))
for(drug in drug_terms$term){
  d = stringr::str_split_fixed(drug, "_", 2)[2]
  cat("running GSEA analysis for DAgenes within allo-activated T-cells affected by", d)
  
  DAgene_gsa_res = filter(coefs_ST, term == drug) %>% 
   nest(-term) %>% 
    dplyr::mutate(
      gsares = purrr::map(.f = get_gsa_res, .x = data, hallmarks, ncores=4)
    ) 

  # Extract p values and GSA test statistics into a table
  DAgene_gsa_summary = 
    DAgene_gsa_res %>%
    dplyr::mutate(gsa_summary = purrr::map(.f = GSAsummaryTable, .x = gsares)) %>% 
    unnest(gsa_summary)


  ## Make a heatmap showing the effects of all the drugs on one cell line:
  simple_model_gsa_summary = DAgene_gsa_summary %>%
    dplyr::mutate(up_stat = `Stat (mix.dir.up)`, dn_stat = `Stat (mix.dir.dn)`) %>%
    dplyr::select(term,
                  Name,
                  `Genes (tot)`,
                  `p (mix.dir.up)`,
                  up_stat,
                  `p (mix.dir.dn)`,
                  dn_stat)
  simple_model_gsa_summary = simple_model_gsa_summary %>% 
    dplyr::mutate(gsa_stat = ifelse((is.na(up_stat) == FALSE & `p (mix.dir.up)` < `p (mix.dir.dn)`),
                             up_stat, ifelse(is.na(dn_stat) == FALSE, -dn_stat, 0)
    ))
  # save summary output from GSEA 
  write.table(simple_model_gsa_summary, file = paste0("AlloActivatedTcells_", d, "_Response_GSEAsummary.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 
  
  simple_model_gsa_stat = simple_model_gsa_summary %>% 
    dplyr::select(Name, term, gsa_stat) %>% 
    spread(term, gsa_stat)

  gsa_effects = simple_model_gsa_stat

  # Get the top gene sets for both
  num_top_gene_sets = 10

  top_simple_gene_sets =
    simple_model_gsa_summary %>%
    group_by(term) %>%
    top_n(num_top_gene_sets, abs(gsa_stat))

  top_simple_gene_sets =
    top_simple_gene_sets %>%
    arrange(term, gsa_stat) %>%
    dplyr::select(Name, term, gsa_stat) %>%
    dplyr::mutate(Name = stringr::str_split_fixed(Name, "%", 3)[, 1], 
           direction = ifelse(gsa_stat > 0, "up", "down"), 
           Name_short = stringr::str_split_fixed(Name, "HALLMARK_", 2)[, 2]) %>% 
    as.data.frame
  sig_pathway_names = unique(as.character(top_simple_gene_sets$Name))

  # prepare barplot of Gene set enrichments
  top_simple_gene_sets$Name_fact <- 
    factor(top_simple_gene_sets$Name_short, levels = top_simple_gene_sets$Name_short[order(top_simple_gene_sets$gsa_stat)])
  
  ggplot(top_simple_gene_sets, aes(x = Name_fact, y = gsa_stat, fill = direction)) +
    geom_col() +
    scale_fill_brewer(palette = "Set1") +
    xlab("pathway") + 
    ylab("directionality") +
    coord_flip() +
    theme_bw()

  ggsave(file = paste0("AlloActivatedTcells_", d, "_Response_GSEA_barplot.pdf"), width = 5, height = 3)
}


##### GSEA heatmap #############
# Make into a matrix
gsa_effect_matrix = gsa_effects %>% dplyr::select(-Name) %>% as.matrix
gsa_effect_matrix[is.na(gsa_effect_matrix)] = 0
row.names(gsa_effect_matrix) = stringr::str_split_fixed(gsa_effects$Name, "%", 3)[, 1]

# Make into a data frame
gsa_effect_df =
  gsa_effects %>%
  gather(key = term, value = gsa_stat, contains("stim"))

gsa_effect_df[is.na(gsa_effect_df)] = 0
gsa_effect_df$Name = stringr::str_split_fixed(gsa_effects$Name, "HALLMARK_", 2)[, 2]
gsa_effect_df$Name = stringr::str_replace(gsa_effect_df$Name, "_", " ")


gsa_effect_matrix[gsa_effect_matrix < -10] = -10
gsa_effect_matrix[gsa_effect_matrix > 10] = 10

matrix_to_plot = as.matrix(gsa_effect_matrix[sig_pathway_names, ])

colnames(matrix_to_plot)  = (colnames(matrix_to_plot) %>% stringr::str_split_fixed("_", 2))[, 2]

rn =   stringr::str_split_fixed(rownames(matrix_to_plot), pattern = "HALLMARK_", n =  2)[, 2]
rownames(matrix_to_plot) = stringr::str_replace_all(rn, "_", " ")
pheatmap(
  matrix_to_plot,
  cex = 0.5,
  filename = "Pseudotime_GOheatmap.pdf",
  width = 3,
  height = 2.25,
  cellwidth = 8,
  cellheight = 8,
  fontsize = 10,
  fontsize_row = 25,
  fontsize_col = 25,
  legend = F,
  cluster_cols = FALSE,
  treeheight_row = 15,
  treeheight_col = 15,
  border_color = "black"
)
##########################

# write bed files for DA peaks for each drug. 
Peak_coefs_Drug = read.table(
  file =  paste0(out_dir, "results/NB9/sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DApeaks.txt"), 
  head = TRUE) %>% 
  dplyr::filter(grepl("dose", term))

drug_terms =  c(distinct(Peak_coefs_Drug, term))

for (drug in drug_terms$term){
  DAsites = dplyr::filter(Peak_coefs_Drug, term == drug, p_value < 0.01) %>% 
    dplyr::mutate("chr" = stringr::str_split_fixed(peak,"_", 3)[,1],
           "start" = stringr::str_split_fixed(peak, "_", 3)[,2],
           "end" = stringr::str_split_fixed(peak, "_", 3)[,3]) %>% 
    dplyr::select(chr, start, end, term, estimate, p_value, q_value) 
  write.table(DAsites, file = paste0(drug, "_Affected_DAsites_inAlloActivatedTcells.bed"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE) 
}

head(DAsites)


# write bed files for DA peaks for each drug. 
coefs_ST = read.table( 
  file =  paste0(out_dir, "results/NB9/sciPlex_ATAC3_MLR_Drug_AlloActivatedTcell_DAgenes.txt"), 
  head = TRUE) %>% 
  dplyr::filter(grepl("dose", term))


drug_terms =  c(distinct(coefs_ST, term))

for (drug in drug_terms$term){
  DAgenes = dplyr::filter(coefs_ST, term == drug, p_value < 0.01) %>% 
    dplyr::mutate("chr" = stringr::str_split_fixed(gene_coords,"_", 3)[,1],
                  "start" = stringr::str_split_fixed(gene_coords, "_", 3)[,2],
                  "end" = stringr::str_split_fixed(gene_coords, "_", 3)[,3]) %>% 
    dplyr::select(chr, start, end, gene_id, term, estimate, p_value, q_value) 
  write.table(DAgenes, file = paste0(drug, "_Affected_DAgenes_inAlloActivatedTcells.bed"),
              sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE) 
}

head(DAgenes)


