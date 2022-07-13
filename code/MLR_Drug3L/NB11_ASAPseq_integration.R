basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/211103_3Level_scichem_MLR_Drugs/"
basepath_ASAP = "/net/trapnell/vol1/home/gtb7/shared/nobackup/downloaded_data/GSE156478_ASAPseq/"

out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB11"))
setwd(paste0(out_dir, "results/NB11"))

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
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)



#################### Only perform once #########################
# build archr project from Hg38 aligned MLR data 
set.seed(1)
inputFiles_mlr = c(mlr_drug = paste0(basepath, "pipeline_output/analysis_MLR_hg38/get_unique_fragments/MLR_5um_homeTn5.fragments2.tsv"))

inputFiles_mlr

addArchRGenome("hg38")

# create arrow file 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles_mlr,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 500,
  bcTag = "qname",
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

# initiate ArchR project object
prj <- ArchRProject(ArrowFiles = ArrowFiles, 
                    copyArrows = FALSE,
                    outputDirectory =  "mlr_drug_raw")

getAvailableMatrices(prj)








# build archr project from ASAP data 
set.seed(1)
inputFiles = c(ASAPseq_pbmc = paste0(basepath_ASAP, "GSM4732109_CD28_CD3_control_ASAP_fragments.tsv.gz"))

inputFiles

addArchRGenome("hg38")

# create arrow file 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 500,
  bcTag = "qname",
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)

# initiate ArchR project object
prj_ASAP <- ArchRProject(ArrowFiles = ArrowFiles, 
                    copyArrows = FALSE,
                    outputDirectory =  "ASAP_raw")

getAvailableMatrices(prj)


################################################################
ADT_data = read.table(file = paste0(basepath_ASAP, "GSM4732109_CD28_CD3_control_ASAP_ADT.tsv.gz"))



