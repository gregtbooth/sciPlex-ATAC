# on cluster, initiate qlogin session with 16 cores. 
#  qlogin -q trapnell-login.q -l mfree=4G -pe serial 16
# start R (3.5.2)

basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/211103_3Level_scichem_MLR_Drugs/"
Fragpath = paste0(basepath, "pipeline_output/analysis_MLR/get_unique_fragments/")
out_dir = paste0(basepath, "analysis/archr/")
dir.create(out_dir)
setwd(out_dir)

library(ArchR)
set.seed(1)
inputFiles = c(MLR_5um_homeTn5 = paste0(Fragpath, "MLR_5um_homeTn5.fragments.txt.gz"),
               MLR_DiagenodeTn5 = paste0(Fragpath, "MLR_DiagenodeTn5.fragments.txt.gz"),
               MLR_homeTn5 = paste0(Fragpath, "MLR_homeTn5.fragments.txt.gz"),
               MLR_IlluminaTn5 = paste0(Fragpath, "MLR_IlluminaTn5.fragments.txt.gz"))

inputFiles

addArchRGenome("hg19")

# create arrow file 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 500,
  bcTag = "qname",
  #gsubExpression = ":.*", #If BamFile : Cell Barcode is everything after ":"
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)




