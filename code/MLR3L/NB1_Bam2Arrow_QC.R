# on cluster, initiate qlogin session with 16 cores. 
#  qlogin -q trapnell-login.q -l mfree=4G -pe serial 16
# start R (3.5.2)

basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
setwd(out_dir)

library(ArchR)
set.seed(1)
inputFiles = c(MLR = paste0(basepath, "pipeline_hill/analysis_MLR/merge_bams/MLR.q10.sorted.merged.bam"))
#inputFiles <- getTutorialData("Hematopoiesis")
inputFiles

addArchRGenome("hg19")
# create arrow file 

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 500, 
  bcTag = "qname",
  gsubExpression = ":.*", #If BamFile : Cell Barcode is everything after ":"
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
# note that the output gives qc_metrix for each input sample 




