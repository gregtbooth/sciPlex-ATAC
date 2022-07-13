basepath = "/home/"
out_dir = paste0(basepath, "analysis/archr/")
setwd(out_dir)

suppressPackageStartupMessages(c(
  library(ArchR),
  library(dplyr),
  library(monocle3), 
  library(ggrastr)
))

source("/src/atac_helper_functions.R")


prefix = "MLR_drug"

addArchRGenome("hg19")
set.seed(1)

ArrowFiles = c("MLR_IlluminaTn5" = paste0(out_dir, "MLR_IlluminaTn5.arrow"),
               "MLR_homeTn5" =  paste0(out_dir, "MLR_homeTn5.arrow"),
               "MLR_DiagenodeTn5" =  paste0(out_dir,"MLR_DiagenodeTn5.arrow"),
               "MLR_5um_homeTn5" =  paste0(out_dir,"MLR_5um_homeTn5.arrow"))

# initiate ArchR project object
prj <- ArchRProject(ArrowFiles = ArrowFiles, 
                    copyArrows = FALSE,
                    outputDirectory = paste0(out_dir, prefix, "_raw"))

getAvailableMatrices(prj)

# calculate doublet likelihood metrics
prj <- addDoubletScores(input = prj,
                         k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                         knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
                         LSIMethod = 3, 
                         nTrials = 5)

# add Hash-based cell data 
# conform cell names to ArchR convention

hash_assmnts = read.table(file = paste0(basepath, "analysis/hash/hashCellAssignments.txt"), 
                          head = TRUE, 
                          stringsAsFactors = FALSE)

hash_assmnts = hash_assmnts %>% 
  mutate(top_to_second_best_ratio_fix = ifelse(is.infinite(top_to_second_best_ratio), 
                                               hash_umis, top_to_second_best_ratio),
         Responder = as.character(stringr::str_split_fixed(top_oligo, "_", 7)[,1]),
         Stimulator = as.character(stringr::str_split_fixed(top_oligo, "_", 7)[,2]), 
         Drug = as.character(stringr::str_split_fixed(top_oligo, "_", 7)[,3]),
         Dose = as.character(stringr::str_split_fixed(top_oligo, "_", 7)[,4]),
         Unit = as.character(stringr::str_split_fixed(top_oligo, "_", 7)[,5]),
         Replicate = as.character(stringr::str_split_fixed(top_oligo, "_", 7)[,6]),
         Relative_Dose = as.character(stringr::str_split_fixed(top_oligo, "_", 7)[,7]), 
         Dose_Unit = paste(Dose, Unit, sep = "_"))

hash_assmnts_cd = data.frame(prj@cellColData) %>% 
  tibble::rownames_to_column(var = "cellName") %>% 
  mutate(barcode = stringr::str_split_fixed(cellName, "#", 2)[,2]) %>% 
  left_join(hash_assmnts, by = c("barcode" = "Cell")) %>% 
  tidyr::replace_na(list(hash_umis = 1, pval = 1, qval = 1, top_to_second_best_ratio_fix = 1))


#filter hash assigment table for cells in ArchR project. 
idx = BiocGenerics::which(hash_assmnts_cd$cellName %in% row.names(prj@cellColData))
HA = hash_assmnts_cd[idx,]

# add hash info columns to project
prj <- addCellColData(ArchRProj = prj, data = HA$hash_umis,
                      cells = HA$cellName, name = "hash_umis", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$pval,
                      cells = HA$cellName, name = "hash_pval", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$qval,
                      cells = HA$cellName, name = "hash_qval", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$top_to_second_best_ratio_fix,
                      cells = HA$cellName, name = "top_to_second_best_ratio", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$top_oligo,
                      cells = HA$cellName, name = "top_oligo", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$doublet,
                      cells = HA$cellName, name = "hash_doublet", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Responder,
                      cells = HA$cellName, name = "Responder", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Stimulator,
                      cells = HA$cellName, name = "Stimulator", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Drug,
                      cells = HA$cellName, name = "Drug", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Dose,
                      cells = HA$cellName, name = "Dose", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Dose_Unit,
                      cells = HA$cellName, name = "Dose_Unit", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Replicate,
                      cells = HA$cellName, name = "Replicate", 
                      force = TRUE)
prj <- addCellColData(ArchRProj = prj, data = HA$Relative_Dose,
                      cells = HA$cellName, name = "Relative_Dose", 
                      force = TRUE)


# save pre-filtered archr project object
prj = saveArchRProject(ArchRProj= prj, 
                       load = TRUE)

################# 
# Run Custom scrublet 
# retrieve tile x cell matrix 
BMAT = getMatrixFromProject(ArchRProj = prj, 
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



# compare doublet metrics. 
#df = data.frame(prj@cellColData) %>% 
#  filter(!is.na(hash_doublet))

#pdf(paste0(out_dir,"sc2_raw/Plots/DoubletEnrich_vs_hashEnrich.pdf"), 
#    width = 2, 
#    height = 1.25)
#ggplot(data = df, aes(x = log10(top_to_second_best_ratio), y = DoubletEnrichment)) +
#  geom_point(size = 0.2, stroke = 0) +
#  theme(legend.position = "right", text = element_text(size = 6),  
#        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
#  xlab("log10(Hash Enrich. score)") +
#  ylab("LSI Doublet Enrich.") 
#dev.off()

#pdf(paste0(out_dir,"sc2_raw/Plots/Boxplot_hash_vs_DoubletEnrich.pdf"), width = 3, height = 3)
#ggplot(df, aes(x = hash_doublet, y = DoubletEnrichment, fill = hash_doublet)) + 
#  geom_boxplot()
#dev.off()

###################
# filter based on LSI doublets
#prj_f = filterDoublets(ArchRProj = prj,
#                       cutEnrich = 1,
#                       cutScore = -Inf,
#                       filterRatio = 1)
#prj_f

# final cell filters. 
cdat = prj@cellColData
fcells = cdat[!is.na(cdat$top_to_second_best_ratio) &
                cdat$hash_doublet == "singlet" & 
                cdat$hash_umis >= 10 &
                cdat$top_to_second_best_ratio >= 2 &
                cdat$TSSEnrichment >= 4 &
                cdat$nFrags >= 500,] 
idx = BiocGenerics::which(prj$cellNames %in% row.names(fcells))
prj_ff <- prj[idx]

saveArchRProject(ArchRProj= prj_ff, 
                 outputDirectory = paste0(prefix, "_filtered"),
                 overwrite = TRUE,
                 load = FALSE)

####################
# run iterative LSI on tile matrix (only filtered cells)
prj_ff <- addIterativeLSI(ArchRProj = prj_ff, 
                       useMatrix = "TileMatrix", 
                       name = "IterativeLSI", 
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
               name = "UMAP",
               nNeighbors = 40,
               minDist = 0.4,
               metric = "cosine",
               force = TRUE)

# plot UMAPs to pdf 
# by cluster
p1 <- plotEmbedding(ArchRProj = prj_ff, 
                    colorBy = "cellColData", 
                    name = "Clusters",
                    embedding = "UMAP")
# by target gene (hash)
p2 <-  plotEmbedding(ArchRProj = prj_ff, 
                     colorBy = "cellColData", 
                     name = "target_gene",
                     embedding = "UMAP")
# by Dex time point (hash)
p3 <-  plotEmbedding(ArchRProj = prj_ff, 
                     colorBy = "cellColData", 
                     name = "dex_time",
                     embedding = "UMAP")

# print to pdf
plotPDF(p1,p2,p3,
        name = "Plot-TileUMAP-LSICusters-treatment.pdf",
        ArchRProj = prj_ff,
        addDOC = FALSE, 
        width = 5, 
        height = 5)

# save final filtered ArchR project to new directory
saveArchRProject(ArchRProj= prj_ff, 
                 overwrite = TRUE,
                 load = FALSE)


####################################################
####################################################

# Port bin x cell matrix to Monocle and store monocle3 embedding in ArchR project 
prj_ff = loadArchRProject(paste0("archr/", prefix, "_filtered"))

bMat = getMatrixFromProject(ArchRProj = prj_ff, 
                          useMatrix = "TileMatrix", 
                          binarize = TRUE)

# format rowData from peakset
bins = data.frame(rowData(bMat)) %>% 
  dplyr::mutate(end = start + 499, 
                bin = paste(seqnames, start, end, sep = "_")) %>% 
  dplyr::select(bin)
row.names(bins) <- bins$bin
row.names(bMat) <- bins$bin

# Create CDS from peak Matrix (SummarizedExperiment)
cds_archr = monocle3::new_cell_data_set(assays(bMat)$TileMatrix, 
                                        cell_metadata = colData(bMat),
                                        gene_metadata = bins)
#########################

set.seed(2017) 

# select bins accessible in > 0.1% of cells
cds_pl <- detect_genes(cds_archr)
ncells = ncol(cds_pl)*0.002
cds_pl = cds_pl[rowData(cds_pl)$num_cells_expressed > ncells,]


##########################
# preproces with monocle3 (default)
cds_pl <- estimate_size_factors(cds_pl)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
reducedDim(cds_pl) <- reducedDim(cds_pl)[,2:50]
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl = cluster_cells(cds_pl)

# custom plot
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))

pdf(paste0(out_dir,"archr/dexTC_filtered/Plots/Monocle3_bMat_UMAP_byTn5Condition.pdf"), width = 3, height = 2.5)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point_rast(size = 0.4, stroke = 0) +
  theme(legend.position = "right", text = element_text(size = 12),  
        legend.key.width = unit(0.5,"line"), legend.key.height = unit(0.5,"line")) + 
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts()
dev.off()


#########################################
# assess hashing metrics

# nHashed cells 
cdat = data.frame(prj@cellColData)

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

###############
# Facet by Tn5 conditions: 
# TSS enrichment
cdat %>% 
  dplyr::filter(!is.na(hash_umis)) %>% 
  ggplot() +
  geom_boxplot(aes(x = Sample, y = log10(TSSEnrichment))) +
  ylab("log10(TSS Enrichment)") +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=0.5))
ggsave(filename = paste0(out_dir,"results/NB2/boxplot_TSS_enrichment_byTn5.pdf"), width = 2.5, height = 3)

# Fragments per cell
cdat %>% 
  dplyr::filter(!is.na(hash_umis)) %>% 
  ggplot() +
  geom_boxplot(aes(x = Sample, y = log10(nFrags))) +
  ylab("log10(nFrags)") +
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=0.5))
ggsave(filename = paste0(out_dir,"results/NB2/boxplot_Frags_byTn5.pdf"), width = 2.5, height = 3)

# Recovered Cells
cdat %>% 
  dplyr::group_by(Sample) %>% 
  dplyr::summarise(cells = n()) %>% 
  ggplot() +
  geom_col(aes(x = Sample, y = log10(cells))) +
  ylab("log10(cells)") + 
  xlab(NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=0.5))
ggsave(filename = paste0(out_dir,"results/NB2/BarPlot_Cells_byTn5.pdf"), 
       width = 2.25, height = 3)

