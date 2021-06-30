basepath = "github/"
bin_directory = paste0(basepath, "analysis/bin/")
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB000"))
setwd(paste0(out_dir, "results/NB000/"))

suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(dplyr)
  library(cicero)
  library(Matrix)
  library(ggrastr)
})

prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered"),
                       showLogo = FALSE)


getAvailableMatrices(prj)


#############
# create a peak x cell CDS and raw matrix files

# load PeakMatrix into memory
pMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"sc2_filtered/ArrowFiles/sc2.arrow"),
                          ArchRProj = prj, 
                          cellNames = prj$cellNames,
                          useMatrix = "PeakMatrix", 
                          binarize = TRUE)

# Create CDS from peak Matrix (SummarizedExperiment)
rd = data.frame(rowRanges(pMat)) %>% 
  tibble::rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end))

# peaks in this list are in a DIFFERENT ORDER 
# relative to pMat (rowRanges above is correct)
pd = data.frame(getPeakSet(prj)) %>% 
  tibble::rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end), 
         gene_short_name = ifelse(peakType == "Promoter" & distToGeneStart < 500, nearestGene, "NA")) %>% 
  dplyr::select(rowname, score, replicateScoreQuantile, groupScoreQuantile, 
                Reproducibility, GroupReplicate, distToGeneStart,
                nearestGene, peakType, distToTSS, nearestTSS,
                gene_short_name, GC, idx, peak)

a_rd = left_join(rd, pd, by = "peak") %>% 
  select(peak, score, replicateScoreQuantile, groupScoreQuantile, 
         Reproducibility, GroupReplicate, distToGeneStart,
         nearestGene, peakType, distToTSS, nearestTSS,
         gene_short_name, GC)

row.names(a_rd) <- a_rd$peak
row.names(pMat) <- a_rd$peak

cds_p = monocle3::new_cell_data_set(assays(pMat)$PeakMatrix, 
                                       cell_metadata = colData(pMat), 
                                       gene_metadata = a_rd)

# write processed files to folder
#cds object (peaks)
saveRDS(cds_p, file = "sciPlex2_peak_cds.RDS")

# peak set
peakset = a_rd
write.table(peakset, file = "sciPlex2_peakset.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# cell Meta data 
cellMetaData = data.frame(colData(cds_p)) %>% 
  tibble::rownames_to_column(var = "cell")
write.table(cellMetaData, file = "sciPlex2_cellMetaData.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)


# create raw matrix (sparse) files
peak_by_cell_matrix_sp = Matrix(assays(pMat)$PeakMatrix, sparse = TRUE)
peak_by_cell_matrix_rows = row.names(assays(pMat)$PeakMatrix)
peak_by_cell_matrix_cols = colnames(assays(pMat)$PeakMatrix)

writeMM(peak_by_cell_matrix_sp,file='sciPlex2_peak_matrix.mtx.txt')
write.table(peak_by_cell_matrix_rows, file='sciPlex2_peak_matrix.rows.txt', sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(peak_by_cell_matrix_cols, file='sciPlex2_peak_matrix.cols.txt', sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)



#################################
# Create CDS from gene score Matrix (SummarizedExperiment)
# load GeneScoreMatrix into memory
gMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"sc2_filtered/ArrowFiles/sc2.arrow"),
                          ArchRProj = prj, 
                          cellNames= prj$cellNames, # restricts to filtered cells
                          useMatrix = "GeneScoreMatrix", 
                          binarize = FALSE)

# Create CDS from peak Matrix (SummarizedExperiment)
g_rd = data.frame(rowData(gMat)) %>% 
  dplyr::mutate(gene_coords = paste0(seqnames, "_", start, "_", end)) %>% 
  dplyr::rename(gene_short_name = name) %>% 
  dplyr::select(gene_coords, gene_short_name, idx)
row.names(g_rd) <- g_rd$gene_short_name
row.names(gMat) <- g_rd$gene_short_name

cds_g = monocle3::new_cell_data_set(assays(gMat)$GeneScoreMatrix, 
                                    cell_metadata = colData(gMat), 
                                    gene_metadata = g_rd)

#cds object (GeneScores)
saveRDS(cds_g, file = "sciPlex2_GeneScores_cds.RDS")

# create raw matrix (sparse) files
gene_by_cell_matrix_sp = Matrix(assays(gMat)$GeneScoreMatrix, sparse = TRUE)
gene_by_cell_matrix_rows = row.names(assays(gMat)$GeneScoreMatrix)
gene_by_cell_matrix_cols = colnames(assays(gMat)$GeneScoreMatrix)

writeMM(peak_by_cell_matrix_sp,file='sciPlex2_GeneScore_matrix.mtx.txt')
write.table(peak_by_cell_matrix_rows, file='sciPlex2_GeneScore_matrix.rows.txt', sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(peak_by_cell_matrix_cols, file='sciPlex2_GeneScore_matrix.cols.txt', sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)


#################################
# get peak x motif matrix 

pull_peak_by_motif_mat = function(archr_project){
  motif_dat = readRDS(paste0(getOutputDirectory(archr_project), "/Annotations/Motif-In-Peaks-Summary.rds"))
  motif_mat = assay(motif_dat$motifMatches)
  peaks = data.frame(rowRanges(motif_dat$motifMatches)) %>% # GB Corrected 11/17/2020
    dplyr::mutate(peak = paste(seqnames, start, end, sep = "_"))
  row.names(motif_mat) <- peaks$peak
  return(motif_mat)
}

motif_mat = pull_peak_by_motif_mat(prj)

# make numeric
mmf = matrix(as.numeric(motif_mat), nrow = nrow(x=motif_mat))

peak_motif_matrix_sp = Matrix(mmf, sparse = TRUE)
peak_motif_matrix_rows = row.names(motif_mat)
peak_motif_matrix_cols = colnames(motif_mat)

writeMM(peak_motif_matrix_sp,file='peak_motif_matrix.mtx.txt')
write.table(peak_motif_matrix_cols, file='peak_motif_matrix.cols.txt', sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(peak_motif_matrix_rows, file='peak_motif_matrix.rows.txt', sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)


##############################
##############################
# example) loading sparse matrices to cds and preprocessing to UMAP

pm = readMM("sciPlex2_peak_matrix.mtx.txt")
rn = read.table("sciPlex2_peak_matrix.rows.txt", header = FALSE, comment.char = "")
cn = read.table("sciPlex2_peak_matrix.cols.txt", header = FALSE, comment.char = "")
peaks = read.table("sciPlex2_peakset.txt", header = TRUE)
rownames(peaks) <- peaks$peak
meta_data = read.table("sciPlex2_cellMetaData.txt", head = TRUE,comment.char = "")
rownames(meta_data) <- meta_data$cell

pm_num = pm*1
rownames(pm_num) <- rn$V1
colnames(pm_num) <- cn$V1

new_cds = monocle3::new_cell_data_set(
  pm_num,
  cell_metadata = meta_data, 
  gene_metadata = peaks)

# preprocess
set.seed(2017)
cds_pl <- detect_genes(new_cds)
cds_pl <- estimate_size_factors(cds_pl)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# custom plot
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))
TCdat_pr$dose_character <- as.character(TCdat_pr$Relative_dose)
TCdat_pr$dose_character <- factor(TCdat_pr$dose_character,
                                  levels = c("0", "0.1", "0.5", "1", "5","10", "50", "100")) 

#pdf by drug
pdf("UMAP_peak_drug.pdf", width = 2, height = 1.25)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = treatment)) +
  #geom_point(size = 0.2, stroke = 0) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  scale_color_manual("Treatment",
                     labels = c("BMS" = "BMS345541", "Dex" = "Dex",
                                "Nutlin" = "Nutlin3A", "SAHA" = "SAHA"),
                     values = c("BMS" = "dimgrey", "Dex" = "deepskyblue3",
                                "Nutlin" = "firebrick3", "SAHA" = "springgreen4")) +
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts()
dev.off()







