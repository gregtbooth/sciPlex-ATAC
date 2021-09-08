suppressPackageStartupMessages({
  library(devtools)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  })

binpath ="/net/trapnell/vol1/home/gtb7/scripts/GB_src/"
basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
source(paste0(binpath,"chiSq_2lvl_GB.R"))
out_dir = paste0(basepath,"analysis/archr/hash/")

samples = c("MLR")
#########
# rbind all samples count tables and define good cells and background cells based on UMI threshold
cell_counts <- data.frame(cell=character(),
                          total=numeric(), 
                          total_deduplicated=numeric(),
                          total_deduplicated_peaks=numeric(), 
                          total_deduplicated_tss=numeric(), 
                          stringsAsFactors=FALSE) 
for (s in samples){
  df = read.table(file = paste0(basepath, "pipeline_hill/analysis_MLR/count_report/", s, ".count_report.txt"), head = T)
  cell_counts = rbind(cell_counts, df)
}

# UMI cutoff is set in pipeline based on knee
cutoff = 182
# get background cells (below cutoff)
background_cells = filter(cell_counts, total_deduplicated < cutoff)
background_cells = as.character(background_cells[,1]) 
# get good cells (above cutoff)
real_test_cell_df = filter(cell_counts, total_deduplicated >= cutoff)
test_cells = as.character(real_test_cell_df[,1])

#########
# Read in hash Oligo count table and reformat Cell IDs 
hashTable = read.csv(paste0(basepath,"pipeline_hill/analysis_hash/hashRDS/hashTable.out"), sep = "\t", header = FALSE)
colnames(hashTable) = c("SampleName", "Cell", "Oligo", "Axis", "Count")

# filter out any hashes from row H (wasn't used in this experiment)
hashTable = filter(hashTable, !grepl("H", Oligo))

###############
# Run chi-squared test and asign hash labels to Cells

# Set up containers for result 
hash_df_all <- data.frame(hash_umis=as.numeric(),
                          pval=numeric(), 
                          qval=numeric(),
                          top_to_second_best_ratio = numeric(),
                          top_oligo=character()) 

downsample_rate =1
print(paste0("Starting hash assignment "))
background_cells_curr = background_cells[background_cells %in% hashTable$Cell]
background_sample_hashes = hashTable %>% filter(Cell %in% background_cells_curr, Axis == 1)
test_cell_hashes = subset(hashTable, Cell %in% test_cells & Axis == 1) 
hash_df_tmp = assign_hash_labels(test_cell_hashes, background_sample_hashes, downsample_rate=downsample_rate)
hash_df_all = rbind(hash_df_all, hash_df_tmp)

##################
# map hash labels to samples
hash_df_all_f = filter(hash_df_all, !is.na(top_oligo))
hash_df_all_f$Plate = stringr::str_split_fixed(hash_df_all_f$top_oligo, pattern = "_" ,n = 2)[,1]
hash_df_all_f$hashWell = stringr::str_split_fixed(hash_df_all_f$top_oligo, pattern = "_" ,n = 2)[,2]
hash_df_all_f$hashRow = substr(hash_df_all_f$hashWell,1,1)
hash_df_all_f$hashCol = substr(hash_df_all_f$hashWell,2,3)

Responder_code = function(hashplate, column, row){
  if (hashplate == "plate10"){
    if(row %in% c("A", "B", "C", "D", "E", "F")){
      if(column %in% c(1,2,3)) return("respA")
      if(column %in% c(4,5,6)) return("respB") 
      if(column %in% c(7,8,9)) return("respC")
      if(column %in% c(10,11,12)) return("respD")}
    else return("stimAlone")
  }}

Stimulator_code = function(hashplate, column, row){
  if (hashplate == "plate10"){
    if(column %in% c(1,2,3)){
      if(row == "A") return("noStim")
      if(row == "B") return("stimA")
      if(row == "C") return("stimB")
      if(row == "D") return("stimC")
      if(row == "E") return("stimD")
      if(row == "F") return("Bead")
      if(row == "G") return("stimA")}
    if(column %in% c(4,5,6)){
      if(row == "A") return("noStim")
      if(row == "B") return("stimB")
      if(row == "C") return("stimA")
      if(row == "D") return("stimC")
      if(row == "E") return("stimD")
      if(row == "F") return("Bead")
      if(row == "G") return("stimB")}
    if(column %in% c(7,8,9)){
      if(row == "A") return("noStim")
      if(row == "B") return("stimC")
      if(row == "C") return("stimA")
      if(row == "D") return("stimB")
      if(row == "E") return("stimD")
      if(row == "F") return("Bead")
      if(row == "G") return("stimC")}
    if(column %in% c(10,11,12)){
      if(row == "A") return("noStim")
      if(row == "B") return("stimD")
      if(row == "C") return("stimA")
      if(row == "D") return("stimB")
      if(row == "E") return("stimC")
      if(row == "F") return("Bead")
      if(row == "G") return("stimD")}
  }}

Replicate_code = function(hashplate, column, row){
  if (hashplate == "plate10"){
    if(column %in% c(1,4,7,10)) return(1)
    if(column %in% c(2,5,8,11)) return(2) 
    if(column %in% c(3,6,9,12)) return(3)
  }}


isDoublet = function(ratio, qval, FDR = 0.05, cutoff = 2){
  if (ratio > cutoff & qval < FDR) return("singlet")
  else(return("doublet"))
}

hash_df_all_f$Responder <- mapply(Responder_code, hash_df_all_f$Plate, hash_df_all_f$hashCol, hash_df_all_f$hashRow)
hash_df_all_f$Stimulator <- mapply(Stimulator_code, hash_df_all_f$Plate, hash_df_all_f$hashCol, hash_df_all_f$hashRow)
hash_df_all_f$Replicate <- mapply(Replicate_code, hash_df_all_f$Plate, hash_df_all_f$hashCol, hash_df_all_f$hashRow)
hash_df_all_f$doublet = mapply(isDoublet, hash_df_all_f$top_to_second_best_ratio, hash_df_all_f$qval,  FDR = 0.05, cutoff = 2)
# make characters
hash_df_all_f$Responder = as.character(hash_df_all_f$Responder)
hash_df_all_f$Stimulator = as.character(hash_df_all_f$Stimulator)
hash_df_all_f$doublet = as.character(hash_df_all_f$doublet)

hashOut <- select(hash_df_all_f, Cell, hash_umis, pval, qval, 
                  top_to_second_best_ratio, Responder, Stimulator, Replicate, doublet)
write.table(hashOut, file= paste0(out_dir, "hashCellAssignments.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)




