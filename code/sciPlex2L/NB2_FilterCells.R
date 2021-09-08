# This script is intended to isolate high-quality cells 
# with good hash labels for identifying treatment conditions. 
# It also aims to identify and filter doublets. 
# Since doublets should be identifiable by hashing, I also compare 
# several doublet calling approaches with hashing data. 
# 
# Various ATAC and Hashing metrics are evaluated and plotted before/after filtering.
# 
# Ultimately, the filtered cells are written to a new, processed 
# ArchR project which will be used for downstream analysis. 

# load scripts containing custom ATAC scrublet code
source("/net/trapnell/vol1/home/gtb7/sciatac_pipeline/src/sylvia_scripts/scripts/atac_helper_functions.R")
basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190521_scichem2_AllPlates/"
out_dir = paste0(basepath, "analysis/archr_revised/")
setwd(out_dir)
dir.create("results")
dir.create("results/NB2")

# load hash info
hash_assmnts = read.table(file = paste0(out_dir, "hash/hashCellAssignments.txt"), 
                          head = TRUE, 
                          stringsAsFactors = FALSE)

suppressPackageStartupMessages(c(
  library(ArchR),
  library(monocle3),
  library(dplyr), 
  library(VennDiagram),
  library(ggrastr)
))

addArchRGenome("hg19")
set.seed(1)

ArrowFiles =  "sc2.arrow"

####################################################
# Create Raw ArchR project with hashing and doublet info.
####################################################

prj <- ArchRProject(ArrowFiles = ArrowFiles, 
                    copyArrows = TRUE,
                    outputDirectory = "sc2_raw")

getAvailableMatrices(prj)

# calculate ArchR doublet likelihood metrics
prj <- addDoubletScores(input = prj,
                         k = 10, 
                         knnMethod = "UMAP", 
                         LSIMethod = 3, 
                         nTrials = 5)

# add Hash-based cell data 
# conform cell names to ArchR convention
hash_assmnts = hash_assmnts %>% 
  dplyr::mutate(cellNames = paste0("sc2#", cell),
                top_to_second_best_ratio_fix = ifelse(is.infinite(top_to_second_best_ratio), 
                                                      hash_umis, top_to_second_best_ratio))

#filter hash assignment table for cells in ArchR project. 
idx = BiocGenerics::which(hash_assmnts$cellNames %in% row.names(prj@cellColData))
HA = hash_assmnts[idx,]
HA = HA[HA$cell_type == "A549",] # remove lingering barnyard cells.

# add hash info columns to project
prj <- addCellColData(ArchRProj = prj, data = HA$hash_umis,
                      cells = HA$cellNames, name = "hash_umis", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$pval,
                      cells = HA$cellNames, name = "hash_pval", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$qval,
                      cells = HA$cellNames, name = "hash_qval", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$top_to_second_best_ratio_fix,
                      cells = HA$cellNames, name = "top_to_second_best_ratio", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$cell_type,
                      cells = HA$cellNames, name = "cell_type", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$sample,
                      cells = HA$cellNames, name = "sample", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$replicate,
                      cells = HA$cellNames, name = "replicate", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$culture_plate,
                      cells = HA$cellNames, name = "culture_plate", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$well_oligo,
                      cells = HA$cellNames, name = "well_oligo", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$treatment,
                      cells = HA$cellNames, name = "treatment", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$dose,
                      cells = HA$cellNames, name = "dose", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Relative_dose,
                      cells = HA$cellNames, name = "Relative_dose", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$vehicle,
                      cells = HA$cellNames, name = "vehicle", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$vehicle_type,
                      cells = HA$cellNames, name = "vehicle_type", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$doublet,
                      cells = HA$cellNames, name = "hash_doublet", 
                      force = TRUE)

#########################################
# run custom scrublet
#########################################
# retrieve tile x cell matrix 
BMAT = getMatrixFromArrow(ArrowFile = paste0(out_dir,"sc2_raw/ArrowFiles/sc2.arrow"),
                          ArchRProj = prj, 
                          useMatrix = "TileMatrix", 
                          binarize = TRUE)
# remove lowly used and unused genomic bins
BMAT = BMAT[rowSums(BMAT@assays@data$TileMatrix) > (0.005 * ncol(BMAT)),] 
bMat = assay(BMAT)

scrub = atac_scrublet(bmat = bMat, k=NULL, fraction_sim_doublets=1, 
                      estimated_doublet_rate=0.1, dims=2:50)
scrub_res = filter(scrub, simulated_doublet == FALSE)
threshold = quantile(scrub_res$doublet_score, .9)
scrub_res$doublet = sapply(scrub_res$doublet_score, function(x){
  ifelse(x < threshold, "singlet", "doublet")})
# add scrublet metrics to ArchR project
prj <- addCellColData(ArchRProj = prj, data = scrub_res$doublet_score,
                      cells = scrub_res$cell, name = "scrublet_score", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = scrub_res$doublet_likelihood,
                      cells = scrub_res$cell, name = "scrublet_likelihood", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = scrub_res$doublet,
                      cells = scrub_res$cell, name = "scrub_doublet", 
                      force = TRUE)

# save pre-filtered archr project object
prj = saveArchRProject(ArchRProj= prj, 
                       load = TRUE)

################################################### 
# compare doublet metrics. 
################################################### 
prj = loadArchRProject(path = paste0(out_dir, "sc2_raw"), 
                       showLogo = FALSE) 

doublet_info = data.frame(prj@cellColData) %>% 
  filter(!is.na(hash_doublet)) %>% 
  select(DoubletScore, DoubletEnrichment, hash_umis, hash_qval, 
         top_to_second_best_ratio, hash_doublet, scrublet_likelihood,
         scrublet_score, scrub_doublet)


pdf(paste0(out_dir,"results/NB2/Boxplot_hash_vs_archr_doubletEnrich.pdf"), width = 3, height = 3)
ggplot(doublet_info, aes(x = hash_doublet, y = DoubletEnrichment, fill = hash_doublet)) + 
  geom_boxplot()
dev.off()

pdf(paste0(out_dir,"results/NB2/Boxplot_hash_vs_scrubletEnrich.pdf"), width = 3, height = 3)
ggplot(doublet_info, aes(x = hash_doublet, y = scrublet_score, fill = hash_doublet)) + 
  geom_boxplot()
dev.off()

# venn diagram of doublet calls by all methods 
h_d = filter(doublet_info, hash_doublet == "doublet") %>% 
  rownames()
s_d = filter(doublet_info, scrub_doublet == "doublet") %>% 
  rownames()
a_d = filter(doublet_info, DoubletEnrichment > 1) %>% 
  rownames()

venn.diagram(
  x = list(h_d, s_d, a_d),
  category.names = c("hash" , "Scrub" , "archr"),
  filename = paste0(out_dir,"results/NB2/VennDiagram_doublet_calls.png"),
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1)

########################################
# assess hashing metrics

df = data.frame(prj@cellColData)

# Hash Enrichment plots 
# Hash UMIs
filter(df, !is.na(hash_umis)) %>% 
  ggplot() +
  geom_histogram(aes(x = log10(hash_umis))) +
  xlab("log10(Hash UMIs)") +
  ylab("Frequency") +
  theme_bw() +
  geom_vline(xintercept=log10(10), color = "red")
ggsave(filename = paste0(out_dir,"results/NB2/Histogram_hashUMIs.pdf"), 
       width = 2.25, height = 2.25)

# Hash Enrichment
filter(df, !is.na(hash_umis)) %>% 
  ggplot() +
  geom_histogram(aes(x = log10(top_to_second_best_ratio))) +
  xlab("log10(Enrichment Ratio)") +
  ylab("Frequency") +
  theme_bw() +
  geom_vline(xintercept=log10(2), color = "red")
ggsave(filename = paste0(out_dir,"results/NB2/Histogram_hashEnrichment.pdf"), 
       width = 2.25, height = 2.25)


# filter cells based on hashing metrics
fcells = df[df$hash_doublet == "singlet" & 
              df$hash_qval <= 0.05 &
              df$hash_umis >= 10 & 
              df$top_to_second_best_ratio >= 2 &
              df$treatment %in% c("BMS", "Dex", "Nutlin", "SAHA"),]
idx = BiocGenerics::which(prj$cellNames %in% row.names(fcells))

prj_f <- prj[idx]

#######################################
# Before removing any doublets, prepare UMAP 
# and assess distribution of Scrublet doublets. 
#######################################

# load tile matrix 
bMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"sc2_raw/ArrowFiles/sc2.arrow"),
                          ArchRProj = prj_f, 
                          cellNames= prj_f$cellNames,
                          useMatrix = "TileMatrix", 
                          binarize = TRUE)

# format rowData 
rd = data.frame(rowData(bMat)) %>% 
  dplyr::mutate(stop = start + 499,
                bin = paste(seqnames, start, stop, sep = "_")) %>% 
  dplyr::select(bin)
row.names(rd) <- rd$bin
row.names(bMat) <- rd$bin

# Create CDS from tile Matrix (SummarizedExperiment)
cds_b = monocle3::new_cell_data_set(assays(bMat)$TileMatrix, 
                                    cell_metadata = colData(bMat),
                                    gene_metadata = rd)

# preprocess bin cds
# reduce dimensions
set.seed(2017)
cds_pl <- detect_genes(cds_b)
# restrict to bins accessible in >0.1% of cells
ncells = ncol(cds_b)*0.005
cds_pl = cds_pl[rowData(cds_pl)$num_cells_expressed > ncells,]
cds_pl <- estimate_size_factors(cds_pl)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
#cds_pl = align_cds(cds_pl, preprocess_method = "LSI", residual_model_formula_str = "~log(nFrags)")
#reducedDim(cds_pl) <- reducedDim(cds_pl)[,2:50]
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# custom plot
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))
TCdat_pr$dose_character <- as.character(TCdat_pr$Relative_dose)
TCdat_pr$dose_character <- factor(TCdat_pr$dose_character,
                                  levels = c("0", "0.1", "0.5", "1",
                                             "5","10", "50", "100")) 

#pdf by clusters
pdf(paste0(out_dir,"results/NB2/Mon3_UMAP_scrublet.pdf"), width = 2, height = 1.25)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = scrub_doublet)) +
  geom_point_rast(size=0.2, stroke = 0) + 
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts()
dev.off()

###########################################################
###########################################################

# lastly, filter out scrublet-called doublets
df = data.frame(prj_f@cellColData)
fcells = df[df$scrub_doublet == "singlet",]
idx = BiocGenerics::which(prj_f$cellNames %in% row.names(fcells))

prj_ff <- prj_f[idx]

saveArchRProject(ArchRProj= prj_ff,
                 outputDirectory = "sc2_filtered",
                 dropCells = TRUE,
                 overwrite = TRUE,
                 load = TRUE)

####################
# Run ArchR iterative LSI on tile matrix (only filtered cells)

# run iterative LSI 
prj_ff <- addIterativeLSI(ArchRProj = prj_ff, 
                       useMatrix = "TileMatrix", 
                       name = "IterativeLSI", 
                       LSIMethod = 3, # 3 is the method most similar to monocle3
                       varFeatures = 100000,
                       dimsToUse = 1:50,
                       iterations = 2,
                       force = TRUE)

# Add clusters
prj_ff <- addClusters(input = prj_ff, 
                   name = "Clusters",
                   reducedDims = "IterativeLSI", 
                   method = "Seurat",
                   corCutOff = 0.75,
                   knnAssign = 10,
                   force = TRUE)

# Add UMAP embedding
prj_ff <- addUMAP(ArchRProj = prj_ff, 
               reducedDims = "IterativeLSI", 
               name = "archrUMAP",
               nNeighbors = 40,
               minDist = 0.4,
               metric = "cosine",
               force = TRUE)

# plot archrUMAPs to pdf 
# by cluster
p1 <- plotEmbedding(ArchRProj = prj_ff, 
                    colorBy = "cellColData", 
                    name = "Clusters",
                    embedding = "archrUMAP")
# by drug
p2 <-  plotEmbedding(ArchRProj = prj_ff, 
                     colorBy = "cellColData", 
                     name = "treatment",
                     embedding = "archrUMAP")

# print to pdf
plotPDF(p1,p2,
        name = "archrPlot-TileUMAP-LSICusters-Drug.pdf",
        ArchRProj = prj_ff,
        addDOC = FALSE, 
        width = 5, 
        height = 5)

# save final filtered ArchR project to new directory
saveArchRProject(ArchRProj= prj_ff, 
                 overwrite = TRUE,
                 load = FALSE)

##################################################
# plot QC metrics by condition. 
##################################################
prj_ff = loadArchRProject(path = paste0(out_dir, "sc2_filtered"), 
                       showLogo = FALSE) 
cd = as.data.frame(prj_ff@cellColData)

cd$new_treatment_label =
  sapply(cd$treatment, function(x) {
    if (grepl("BMS", x))
      return("BMS345541")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
  })

# number of cells by dose
ncells = group_by(cd, new_treatment_label, Relative_dose) %>% 
  dplyr::summarise(ncells = n())
ncells$Relative_dose = factor(ncells$Relative_dose,
                              levels = c("0", "0.1", "0.5", "1", 
                                         "5", "10", "50", "100"))
ggplot(ncells) +
  geom_col(aes(x = Relative_dose, y = ncells)) +
  facet_wrap("~new_treatment_label") + 
  monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Relative Dose") +
  ylab("cells") + 
  theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/cellsbydose_drugFacet.png"),
       width = 4, height = 3)

# TSS Enrichment by dose
cd$Relative_dose_factor =  
  factor(cd$Relative_dose, levels = c("0", "0.1", "0.5", "1",
                                      "5", "10", "50", "100"))
ggplot(cd) +
  geom_boxplot(aes(x =  Relative_dose_factor, y = TSSEnrichment)) +
  facet_wrap("~new_treatment_label") +
  monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Relative Dose") +
  ylab("TSS Enrichment") +
  theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/TSSEnrich_byDose_drugFacet.png", sep = ""), 
       width = 4, height = 3)

# UMI by dose
ggplot(cd) +
  geom_boxplot(aes(x =  Relative_dose_factor, y = log10(nFrags))) +
  facet_wrap("~new_treatment_label") +
  monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Relative Dose") +
  ylab("log10(Frags per cell)") + 
  theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/UMIbyDose_drugFacet.png", sep = ""),
       width = 4, height = 3)

#Hash UMI by dose
ggplot(cd) +
  geom_boxplot(aes(x =  Relative_dose_factor, y = log10(hash_umis))) +
  facet_wrap("~new_treatment_label") +
  monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Relative Dose") +
  ylab("log10(Hash counts per cell)") + 
  theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/HashUMIbyDose_drugFacet.png", sep = ""),
       width = 4, height = 3)

#Hash Enrichment score by dose
ggplot(cd) +
  geom_boxplot(aes(x =  Relative_dose_factor, y = log10(top_to_second_best_ratio))) +
  facet_wrap("~new_treatment_label") +
  monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Relative Dose") +
  ylab("log10(Hash Enrichment score)") + 
  theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/HashEnrichmentbyDose_drugFacet.png", sep = ""),
       width = 4, height = 3)

#run regression to test significance of log(dose) on various metrics. 
dose_TSSenr_relations = function(df, treat ){
  df_filt = filter(df, treatment == treat) %>% 
    mutate(log_dose = log(dose + 0.01)) #add pseudodose and log
  
  doseCoef = summary(glm("TSSEnrichment~log_dose + nFrags",family = "gaussian", data = df_filt))$coefficients["log_dose", "Estimate"]
  pVal = summary(glm("TSSEnrichment~log_dose + nFrags",family = "gaussian", data = df_filt))$coefficients["log_dose", "Pr(>|t|)"]
  res = data.frame(treatment = treat, dose_coef = doseCoef, p_val = pVal)
  return(res)
}

coef_table = data.frame(treatment = character(), dose_coef = numeric(), p_val = numeric())
for(d in c("SAHA", "Dex", "BMS", "Nutlin")){
    coefs = dose_TSSenr_relations(df = cd, treat = d)
    coef_table = rbind(coef_table, coefs)}


####################################################
# compare cells with hill_filtered set 
load(paste0(basepath, "analysis/archr/monocle3_hill/cds_p"))
hillCells = data.frame(colData(cds_p_h)) 
hillCells = hillCells %>% dplyr::mutate(cellNames = paste0("sc2#", cell))
hillCells_f = hillCells[hillCells$umi >= 1000, ]

idx_hill = BiocGenerics::which(row.names(prj@cellColData) %in% hillCells$cellNames)
prj_t = prj[idx_hill]


cd = data.frame(prj_t@cellColData)
cd$cellNames = row.names(cd)

cd_j = left_join(hillCells, cd, by = "cellNames") %>% 
  select(cellNames, total, umi, nMultiFrags, nMonoFrags, nFrags, nDiFrags) %>% 
  mutate(frac = nFrags/umi)








