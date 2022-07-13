basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/211103_3Level_scichem_MLR_Drugs/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB14"))
setwd(paste0(out_dir, "results/NB14"))

suppressPackageStartupMessages({
  library(ArchR)
  #library(Seurat)
  #library(monocle3)
  library(devtools)
  load_all('/net/trapnell/vol1/home/gtb7/git/github/monocle3_dev') # latest developer branch of monocle3
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
  plan(multicore)
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered_annotated"), 
                       showLogo = FALSE) 

# Load CDS with T-cell trajectory pseudotime information (then pull only allo-activated cells)
cds = readRDS(file = paste0(out_dir, "results/NB8/cds_NaiveActive_Tcells_PT"))
cds_f = cds[,colData(cds)$cellType_broad == "Activated_Tcell" & colData(cds)$Stimulator == "StimA"]
