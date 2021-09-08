#set working environment in R

basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190521_scichem2_AllPlates/"
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB9"))
setwd(paste0(out_dir, "results/NB9/"))

suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

addArchRGenome("hg19")

# load ArchR project 
prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered"),
                       showLogo = FALSE)



# create new project to save integrated info to
#prj_int = saveArchRProject(ArchRProj = prj,
#                           outputDirectory = paste0(basepath, "analysis/archr/sc2_integratedRNA_v2"),
#                           load = TRUE)

prj_int = prj  
# add column to metadata with full treatment info. 
cdat = data.frame(prj_int@cellColData)
cdat$cellNames = row.names(prj_int@cellColData)
cdat = mutate(cdat, treatment_hashbased = paste0(treatment, "_", Relative_dose))
# to boost cell numbers in several treatments, collapse higher doses into one group (BMS_5p and Nutlin_50p)
cdat$treatment_hashmod = ifelse(cdat$treatment_hashbased %in% c("BMS_100", "BMS_50", "BMS_10", "BMS_5"),
                                "BMS_5p", ifelse(cdat$treatment_hashbased %in% c("Nutlin_100", "Nutlin_50"), 
                                                 "Nutlin_50p", ifelse(cdat$treatment_hashbased %in% c("SAHA_0", "Dex_0", "Nutlin_0", "BMS_0"), 
                                                                      "Vehicle", cdat$treatment_hashbased)))

prj_int <- addCellColData(ArchRProj = prj_int, data = cdat$treatment_hashbased,
                      cells = cdat$cellNames, name = "treatment_hashbased", 
                      force = TRUE)
prj_int <- addCellColData(ArchRProj = prj_int, data = cdat$treatment_hashmod,
                      cells = cdat$cellNames, name = "treatment_hashmod", 
                      force = TRUE)

########################
# load cds.rna and convert to seurat object 
#load SciPlex_RNA dataset (Srivatsan et. al. 2020)
# need to change feature names to gene_short_names for matching to ATAC activity matrix
cds.rna <- readRDS(paste0(basepath, "analysis/archr_revised/scRNA/sciPlex2_cds_NB5processed.RDS"))
rd = rowData(cds.rna)
cm = counts(cds.rna)
cd = colData(cds.rna)
row.names(rd) <- rd$gene_short_name
row.names(cm) <- rd$gene_short_name
cds.rna = new_cell_data_set(expression_data = cm, 
                            cell_metadata = cd, 
                            gene_metadata = rd)


# Create Seurat object for RNA dataset
rna.matrix <- counts(cds.rna)
seRNA <- CreateSeuratObject(counts = rna.matrix, assay = "RNA", project = "SC2_RNA")
seRNA <- subset(seRNA, subset = nFeature_RNA > 200)
# add metaData to Seurat object
meta.rna = as.data.frame(colData(cds.rna))
seRNA <- AddMetaData(seRNA, metadata = meta.rna)
seRNA$tech <- "rna"
seRNA$treatment = stringr::str_split_fixed(string = seRNA$top_oligo, pattern = "_", n=3)[,1]
seRNA$Relative_dose = stringr::str_split_fixed(string = seRNA$top_oligo, pattern = "_", n=3)[,2]
seRNA$treatment_RNAbased = paste(seRNA$treatment, seRNA$Relative_dose, sep = "_")
seRNA$treatmentRNA_broad = ifelse(seRNA$Relative_dose == 0, "vehicle", seRNA$treatment)
# to boost cell numbers in several treatments, collapse higher doses into one group (BMS_5p and Nutlin_50p)
seRNA$treatment_RNAmod = ifelse(seRNA$treatment_RNAbased %in% c("BMS_100", "BMS_50", "BMS_10", "BMS_5"),
                                "BMS_5p", ifelse(seRNA$treatment_RNAbased %in% c("Nutlin_100", "Nutlin_50"), 
                                                 "Nutlin_50p", ifelse(seRNA$treatment_RNAbased %in% c("SAHA_0", "Dex_0", "Nutlin_0", "BMS_0"), 
                                                                        "Vehicle", seRNA$treatment_RNAbased)))

##################################
# use hash-based treatment labels to constrain 
# integration with corresponding RNA data

# format a nested Simplelist with groups of cells for each treatment 
# modified to collapse high BMS and nutlin doses 
atac_cDat = data.frame(prj_int@cellColData)
atac_cDat$cellNames = row.names(prj_int@cellColData)

rna_cDat = data.frame(seRNA@meta.data)
rna_cDat$cellNames = row.names(seRNA@meta.data)

groupList <- SimpleList()
for(treat in dplyr::distinct(atac_cDat, treatment_hashmod)$treatment_hashmod){
  print(treat)
  groupList[[treat]] = SimpleList(
    ATAC = atac_cDat$cellNames[atac_cDat$treatment_hashmod == treat],
    RNA = rna_cDat$cellNames[rna_cDat$treatment_RNAmod == treat])
}

prj_int <- addIterativeLSI(ArchRProj = prj_int, 
                          useMatrix = "TileMatrix", 
                          name = "IterativeLSI", 
                          iterations = 2,
                          force = TRUE)

prj_int <- addGeneIntegrationMatrix(
  ArchRProj = prj_int, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE, # keep FALSE until verified that it worked well
  groupList = groupList,
  groupRNA = "treatment_RNAmod",
  nameCell = "predictedCell_Co",
  nameGroup = "predictedGroup_Co",
  nameScore = "predictedScore_Co", 
  k.filter = 50, # this modifies Seurat's CCA method (because I don't have enough cells for default 200)
  force = TRUE
)

# confusion matrix of the constrained predictions vs actual
cM_co <- as.matrix(confusionMatrix(prj_int$treatment_hashmod, prj_int$predictedGroup_Co))
preClust_co <- colnames(cM_co)[apply(cM_co, 1 , which.max)]
cbind(preClust_co, rownames(cM_co)) #Assignments

# transform to fraction of rowSums (total cells/condition)
cM_co_frac <- apply(cM_co, 2, function(i) i/sum(i))
cM_co_melt = melt(cM_co_frac)
colnames(cM_co_melt) = c("hashBased", "predictionBased", "fracPredicted")
# reorder cols and rows
cM_co_melt$hashBased = factor(cM_co_melt$hashBased, levels = c("Vehicle", "BMS_0.1", "BMS_0.5", "BMS_1", "BMS_5p", 
                                                               "Dex_0.1", "Dex_0.5", "Dex_1", "Dex_5", "Dex_10", "Dex_50", "Dex_100", 
                                                               "Nutlin_0.1", "Nutlin_0.5", "Nutlin_1", "Nutlin_5", "Nutlin_10", "Nutlin_50p", 
                                                               "SAHA_0.1", "SAHA_0.5", "SAHA_1", "SAHA_5", "SAHA_10", "SAHA_50", "SAHA_100"))
cM_co_melt$predictionBased = factor(cM_co_melt$predictionBased, levels = c("Vehicle", "BMS_0.1", "BMS_0.5", "BMS_1", "BMS_5p", 
                                                                           "Dex_0.1", "Dex_0.5", "Dex_1", "Dex_5", "Dex_10", "Dex_50", "Dex_100", 
                                                                           "Nutlin_0.1", "Nutlin_0.5", "Nutlin_1", "Nutlin_5", "Nutlin_10", "Nutlin_50p", 
                                                                           "SAHA_0.1", "SAHA_0.5", "SAHA_1", "SAHA_5", "SAHA_10", "SAHA_50", "SAHA_100"))
# plot confusion matrix (all should be on diagonal)
ggplot(cM_co_melt, aes(x = hashBased, y = predictionBased, fill = fracPredicted)) + 
  geom_tile(color = "gray") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  ylab("prediction based") +
  xlab("hash based") + 
  labs(fill = "frac. predicted") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "Confusion_PredVsHash_constrained.pdf",
       width = 7, height = 5)

# save integrated ArchR project
#prj_int = saveArchRProject(ArchRProj = prj_int,
#                           outputDirectory = paste0(basepath, "analysis/archr/sc2_integratedRNA_v2"),
#                           load = TRUE)

##############################
# unconstrained integration of scRNA with scATAC
# not adding to Arrow, but just for assessing assignments
prj_int <- addGeneIntegrationMatrix(
  ArchRProj = prj_int, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "treatment_RNAmod",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un", 
  force = TRUE
)

# confusion matrix of the predictions vs actual
cM <- as.matrix(confusionMatrix(prj_int$treatment_hashmod, prj_int$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments

# transform to fraction of rowSums (total cells/condition)
cM_frac <- apply(cM, 2, function(i) i/sum(i))
cM_melt = melt(cM_frac)
colnames(cM_melt) = c("hashBased", "predictionBased", "fracPredicted")
# reorder cols and rows
cM_melt$hashBased = factor(cM_melt$hashBased, levels = c("Vehicle", "BMS_0.1", "BMS_0.5", "BMS_1", "BMS_5p", 
                                                         "Dex_0.1", "Dex_0.5", "Dex_1", "Dex_5", "Dex_10", "Dex_50", "Dex_100", 
                                                         "Nutlin_0.1", "Nutlin_0.5", "Nutlin_1", "Nutlin_5", "Nutlin_10", "Nutlin_50p", 
                                                          "SAHA_0.1", "SAHA_0.5", "SAHA_1", "SAHA_5", "SAHA_10", "SAHA_50", "SAHA_100"))
cM_melt$predictionBased = factor(cM_melt$predictionBased, levels = c("Vehicle", "BMS_0.1", "BMS_0.5", "BMS_1", "BMS_5p", 
                                                                     "Dex_0.1", "Dex_0.5", "Dex_1", "Dex_5", "Dex_10", "Dex_50", "Dex_100", 
                                                                     "Nutlin_0.1", "Nutlin_0.5", "Nutlin_1", "Nutlin_5", "Nutlin_10", "Nutlin_50p", 
                                                                     "SAHA_0.1", "SAHA_0.5", "SAHA_1", "SAHA_5", "SAHA_10", "SAHA_50", "SAHA_100"))
                                      
# plot confusion matrix
ggplot(cM_melt, aes(x = hashBased, y = predictionBased, fill = fracPredicted)) + 
  geom_tile(color = "gray") +
  scale_fill_gradient(low = "white", high = "steelblue") + 
  ylab("prediction based") +
  xlab("hash based") + 
  labs(fill = "frac. predicted") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "Confusion_PredVsHash_unconstrained.pdf",
       width = 7, height = 5)

saveRDS(cM, file = "treatment_by_prediction_CellCount_mtx")
saveRDS(cM_melt, file = "treatment_by_prediction_fracCells_melt")
#########################

pal <- paletteDiscrete(values = seRNA@meta.data$treatment_RNAmod)

pal

p1 <- plotEmbedding(
  prj_int, 
  embedding = "UMAP_mon3",
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal
)

p2 <- plotEmbedding(
  prj_int, 
  embedding = "UMAP_mon3",
  colorBy = "cellColData", 
  name = "predictedGroup_Co", 
  pal = pal
)

plotPDF(p1,p2, 
        name = "Plot-UMAP-RNA-Integration.pdf", 
        ArchRProj = prj_int, 
        addDOC = FALSE, 
        width = 5, 
        height = 5)


