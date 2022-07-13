## Note: you need to have MACS in your path to performe peak calling: 
## module load MACS/2.2.7.1

basepath = "/home/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB7"))
setwd(paste0(out_dir, "results/NB7"))

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
  library(RColorBrewer)
  library(ggpubr)
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered_annotated"), 
                       showLogo = FALSE) 


##############################################
## Call denovo peaks on full dataset (using broad Cell type annotation) 
#prj <- addGroupCoverages(ArchRProj = prj, groupBy = "cellType_broad")

#prj <- addReproduciblePeakSet(
#  ArchRProj = prj, 
#  groupBy = "cellType_broad", 
#  pathToMacs2 = pathToMacs2
#)

#peaks_denovo = getPeakSet(prj)

#prj <- addPeakSet(ArchRProj = prj, peakSet = peaks_denovo)

#prj <- addPeakMatrix(prj)

#saveArchRProject(ArchRProj = prj, 
#                 load = FALSE)


# print peaks to bed file: 
#peaks = getPeakSet(prj) %>% data.frame() %>% 
#  select(seqnames, start, end, width, score, peakType, GroupReplicate)

#write.table(peaks, file = "MLR_Drug_peaks.bed", col.names = FALSE, quote = FALSE, sep = "\t", row.names =  FALSE)
##################################################

#Note: More peaks were called on original dataset, therefore I'll just use those going forward 
## this should also facilitate comparison with the original dataset.

# load original mlr archr project: 
orig_data_path = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/analysis/archr/"
prj_orig = loadArchRProject(path = paste0(orig_data_path, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 

peaks_orig = getPeakSet(prj_orig)

# add peaks to this (drugged) dataset: 
prj <- addPeakSet(ArchRProj = prj, peakSet = peaks_orig, force = TRUE)

prj <- addPeakMatrix(prj)

saveArchRProject(ArchRProj = prj, 
                 load = FALSE)

########################################
# create pseudobulk bigwig files for 
# all celltypes from all treatment conditions

prj$celltype_condition = paste0(prj$Stimulator, "_", prj$Drug, "_", prj$Relative_Dose, "_", prj$cellType_broad)

# separate tracks for all cell types within each condition (139 total)
getGroupBW(
  ArchRProj = prj,
  groupBy = "celltype_condition",
  normMethod = "ReadsInTSS", # creates a group scale factor = 10k/sum(reads in TSS) 
  #normMethod = "nFrags", # creates a group scale factor = 10k/sum(nFrags)
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

# separate tracks for each cell type (combined conditions)
getGroupBW(
  ArchRProj = prj,
  groupBy = "cellType_broad",
  normMethod = "ReadsInTSS", # creates a group scale factor = 10k/sum(reads in TSS) 
  #normMethod = "nFrags", # creates a group scale factor = 10k/sum(nFrags)
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)

# print peaks to bed file: 
peaks = getPeakSet(prj) %>% data.frame() %>% 
  select(seqnames, start, end, width, score, peakType, GroupReplicate)

write.table(peaks, file = "MLR_peaks.bed", col.names = FALSE, quote = FALSE, sep = "\t", row.names =  FALSE)
