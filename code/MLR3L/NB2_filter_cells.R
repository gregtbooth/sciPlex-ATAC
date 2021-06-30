# This script is intended to isolate high-quality cells 
# with good hash labels for identifying treatment conditions. 
# It also aims to identify and filter doublets. 
# Since doublets should be identifiable by hashing, I also compare 
# several doublet calling approaches with hashing data. 
# 
# Various ATAC and Hashing metrics are evaluated and plotted before/after filtering.
# 
# Ulimately, the filtered cells are written to a new, processed 
# ArchR project which will be used for downstream analysis. 

#set working environment in R
source("/src/atac_helper_functions.R")
basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
setwd(out_dir)
dir.create("results")
dir.create("results/NB2")

suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(VennDiagram)
  library(ggrastr)
})

hash_assmnts = read.table(file = paste0(out_dir, "hash/hashCellAssignments.txt"), 
                          head = TRUE, stringsAsFactors = FALSE)

addArchRGenome("hg19")
set.seed(1)

ArrowFiles = c("MLR" = "MLR.arrow")

####################################################
# section 1
# Create Raw ArchR project with hashing and doublet info for all cells.
# This only needs to be run once.  Can skip to next section for filtering.
####################################################

prj <- ArchRProject(ArrowFiles = ArrowFiles, 
                    copyArrows = TRUE,
                    outputDirectory = "mlr_raw")

getAvailableMatrices(prj)

# calculate ArchR doublet likelihood metrics
prj <- addDoubletScores(input = prj,
                        k = 10, 
                        knnMethod = "UMAP", 
                        LSIMethod = 3, 
                        nTrials = 5)

######################################
# add hash assignment data project 
# adjust hash cell naming convention to archr
hash_assmnts = hash_assmnts %>% 
  dplyr::mutate(cellNames = paste0("MLR#", Cell),
                top_to_second_best_ratio_fix = ifelse(is.infinite(top_to_second_best_ratio), 
                                                      hash_umis, top_to_second_best_ratio))
#filter hash assigment table for cells in ArchR project. 
idx = BiocGenerics::which(hash_assmnts$cellNames %in% row.names(prj@cellColData))
HA = hash_assmnts[idx,]

# add hash info columns to the project based by cellNames
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
prj <- addCellColData(ArchRProj = prj, data = HA$Responder,
                      cells = HA$cellNames, name = "Responder", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Stimulator,
                      cells = HA$cellNames, name = "Stimulator", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Replicate,
                      cells = HA$cellNames, name = "Replicate", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$doublet,
                      cells = HA$cellNames, name = "hash_doublet", 
                      force = TRUE)


#########################################
# run custom scrublet

# retrieve tile x cell matrix 
BMAT = getMatrixFromArrow(ArrowFile = paste0(out_dir,"mlr_raw/ArrowFiles/MLR.arrow"),
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
                       overwrite = TRUE,
                       dropCells = FALSE,
                       load = TRUE)

########################################
# Section 2 
# Filter cells based on hashing and doublet scores
########################################

# load processed ArchR project from NB2
addArchRGenome("hg19")
prj = loadArchRProject(path = paste0(out_dir, "mlr_raw"), showLogo = FALSE)

cdat = data.frame(prj@cellColData)

#compare doublet scores 
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

#########################################
# assess hashing metrics

# nHashed cells 
f1 = cdat %>% 
  dplyr::mutate(hashed = ifelse(!is.na(hash_umis), "hashed", "not_hashed")) %>% 
  dplyr::group_by(hashed) %>% 
  dplyr::summarise(cells = n()) 

ggplot(f1) +
  geom_col(aes(x = hashed, y = cells)) +
  ylab("cells") + 
  theme_bw() 
ggsave(filename = paste0(out_dir,"results/NB2/BarPlot_nHashedCells.pdf"), 
       width = 2.25, height = 2.25)

# hash umis/cell 
cdat %>% 
  filter(!is.na(hash_umis)) %>% 
  ggplot() +
  geom_histogram(aes(x = log10(hash_umis)), bins = 50) +
  xlab("log10(hash_umis)") +
  ylab("cells") +
  geom_vline(xintercept=log10(10), color = "red") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"results/NB2/histogram_hash_umis.pdf"), width = 2.5, height = 2.5)

# hash enrichment
cdat %>% 
  dplyr::filter(!is.na(hash_umis)) %>% 
  dplyr::mutate(hash_enrichment = ifelse(is.infinite(top_to_second_best_ratio), 
                                         hash_umis, top_to_second_best_ratio)) %>% 
  ggplot() +
  geom_histogram(aes(x = log10(top_to_second_best_ratio)), bins = 50) +
  xlab("log10(enrichment ratio)") +
  ylab("cells") +
  geom_vline(xintercept=log10(2), color = "red") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"results/NB2/histogram_hash_enrichment.pdf"), width = 2.5, height = 2.5)

# TSS enrichment
cdat %>% 
  dplyr::filter(!is.na(hash_umis)) %>% 
  ggplot() +
  geom_histogram(aes(x = log10(TSSEnrichment)), bins = 50) +
  xlab("log10(TSS Enrichment)") +
  ylab("cells") +
  geom_vline(xintercept=log10(4), color = "red") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"results/NB2/histogram_TSS_enrichment.pdf"), width = 2.5, height = 2.5)

# Fragments per cell
cdat %>% 
  dplyr::filter(!is.na(hash_umis)) %>% 
  ggplot() +
  geom_histogram(aes(x = log10(nFrags)), bins = 50) +
  xlab("log10(Fragments)") +
  ylab("cells") +
  geom_vline(xintercept=log10(500), color = "red") +
  theme_bw() 
ggsave(filename = paste0(out_dir,"results/NB2/histogram_Frags.pdf"), width = 2.5, height = 2.5)

####################################

# filter cells based on hashing metrics
fcells = cdat[cdat$hash_doublet == "singlet" & 
                cdat$hash_qval <= 0.05 &
                cdat$hash_umis >= 10 & 
                cdat$top_to_second_best_ratio >= 2,] 
idx = BiocGenerics::which(prj$cellNames %in% row.names(fcells))

prj_f <- prj[idx]

#######################################
# section 3: Filter doublets
#######################################
# Before removing any doublets, prepare UMAP 
# and assess distribution of Scrublet doublets. 

# load tile matrix 
bMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"mlr_raw/ArrowFiles/MLR.arrow"),
                          ArchRProj = prj_f, 
                          cellNames= prj_f$cellNames,
                          useMatrix = "TileMatrix", 
                          binarize = TRUE)

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
ncells = ncol(cds_b)*0.002
cds_pl = cds_pl[rowData(cds_pl)$num_cells_expressed > ncells,]
cds_pl <- estimate_size_factors(cds_pl)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
#cds_pl = align_cds(cds_pl, preprocess_method = "LSI", residual_model_formula_str = "~log(nFrags)")
reducedDim(cds_pl) <- reducedDim(cds_pl)[,2:50]
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# custom plot
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))

#pdf by clusters
pdf(paste0(out_dir,"results/NB2/Mon3_UMAP_scrublet_noDim1.pdf"), width = 2, height = 1.25)
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

# lastly, filter out scrublet-called doublets
df = data.frame(prj_f@cellColData)
fcells = df[df$scrub_doublet == "singlet",]
idx = BiocGenerics::which(prj_f$cellNames %in% row.names(fcells))

prj_ff <- prj_f[idx]

saveArchRProject(ArchRProj= prj_ff,
                 outputDirectory = "mlr_filtered",
                 dropCells = TRUE,
                 overwrite = TRUE,
                 load = TRUE)

# add ArchR embeddings to project (usefull for some downstreeam functions)
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
dir.create("results/NB2/Plots")
# by cluster
p1 <- plotEmbedding(ArchRProj = prj_ff, 
                    colorBy = "cellColData", 
                    name = "Clusters",
                    embedding = "archrUMAP")
# by Stimulator
p2 <-  plotEmbedding(ArchRProj = prj_ff, 
                     colorBy = "cellColData", 
                     name = "Stimulator",
                     embedding = "archrUMAP")

# by Responder
p3 <-  plotEmbedding(ArchRProj = prj_ff, 
                     colorBy = "cellColData", 
                     name = "Responder",
                     embedding = "archrUMAP")

# print to pdf
plotPDF(p1,p2,p3,
        name = "archrPlot-TileUMAP-LSICusters-conditions.pdf",
        ArchRProj = prj_ff,
        addDOC = FALSE, 
        width = 5, 
        height = 5)

# save final filtered ArchR project to new directory
saveArchRProject(ArchRProj= prj_ff,
                 outputDirectory = "mlr_filtered",
                 dropCells = TRUE,
                 overwrite = TRUE,
                 load = TRUE)

##################################################
# plot QC metrics by condition. 
##################################################

prj_fff = loadArchRProject(path = paste0(out_dir, "mlr_filtered"), 
                          showLogo = FALSE) 

cd = as.data.frame(prj_fff@cellColData)

# number of cells by Stimulator
ncells = group_by(cd, Responder, Stimulator) %>% 
  dplyr::summarise(ncells = n())

ggplot(ncells) +
  geom_col(aes(x = Stimulator, y = ncells)) +
  facet_wrap("~Responder") + 
  #monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Stimulator") +
  ylab("cells") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  #theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/cells_byStim_ResponderFacet.png"),
       width = 4, height = 3)

# number of cells by Stimulator (broken down by replicates)
ncells_rep = group_by(cd, Responder, Stimulator, Replicate) %>% 
  dplyr::summarise(ncells = n())

ggplot(ncells_rep) +
  geom_boxplot(aes(x =  Stimulator, y = ncells)) +
  facet_wrap("~Responder") +
  #monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Stimulator") +
  ylab("cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
#theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/cells_byStim_ResponderFacet_BoxPlot.png", sep = ""), 
       width = 4, height = 3)

# TSS Enrichment by Stimulator
ggplot(cd) +
  geom_boxplot(aes(x =  Stimulator, y = TSSEnrichment)) +
  facet_wrap("~Responder") +
  #monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Stimulator") +
  ylab("TSS Enrichment") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  #theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/TSSEnrich_byStim_ResponderFacet.png", sep = ""), 
       width = 4, height = 3)

# UMI by dose
ggplot(cd) +
  geom_boxplot(aes(x = Stimulator, y = log10(nFrags))) +
  facet_wrap("~Responder") +
  #monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Stimulator") +
  ylab("log10(Frags per cell)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  #theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/UMIbyStim_ResponderFacet.png", sep = ""),
       width = 4, height = 3)

#Hash UMI by dose
ggplot(cd) +
  geom_boxplot(aes(x = Stimulator, y = log10(hash_umis))) +
  facet_wrap("~Responder") +
  #monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Stimulator") +
  ylab("log10(Hash counts per cell)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  #theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/HashUMIbyStim_ResponderFacet.png", sep = ""),
       width = 4, height = 3)

#Hash Enrichment score by dose
ggplot(cd) +
  geom_boxplot(aes(x =  Stimulator, y = log10(top_to_second_best_ratio))) +
  facet_wrap("~Responder") +
  #monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Stimulator") +
  ylab("log10(Hash Enrichment score)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  #theme_bw()
ggsave(filename = paste0(out_dir,"results/NB2/HashEnrichmentbyStim_ResponderFacet.png", sep = ""),
       width = 4, height = 3)




