# on cluster, initiate qlogin session with 16 cores. 
#  qlogin -q trapnell-login.q -l mfree=4G -pe serial 16
# start R (4.0)

basepath = "github/"
out_dir = paste0(basepath, "analysis/archr_test/")
dir.create(out_dir)
setwd(out_dir)

library(ArchR)
set.seed(1)
inputFiles = c(sc2 = paste0(basepath, "pipeline_hill/analysis_SC2/merge_bams/SC2.q10.sorted.merged.bam"))
#inputFiles = c(sc2 = paste0(basepath, "pipeline_hill/analysis_SC2/get_unique_fragments/SC2.fragments.txt.gz"))
inputFiles

addArchRGenome("hg19")

# create arrow file 
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 3, #Dont set this too high because you can always increase later
  filterFrags = 500, 
  minFrags = 300,
  bcTag = "qname",
  gsubExpression = ":.*", #If BamFile : Cell Barcode is everything after ":"
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)




