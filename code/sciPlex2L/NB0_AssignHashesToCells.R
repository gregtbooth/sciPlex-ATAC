
suppressPackageStartupMessages({
  library(devtools)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
})

source("~/scripts/GB_src/ConvertBCs_ToWells_V3.R")
source("~/scripts/GB_src/chiSq_2lvl_GB.R")

basepath = "github/"
out_dir = paste0(basepath,"analysis/archr_revised/hash/")

samples = c(SC2 = paste0(basepath, "pipeline_hill/analysis_SC2/count_report/SC2.count_report.txt"),
            BY = paste0(basepath, "pipeline_hill/analysis_barnyard/count_report/barnyard.count_report.txt"))
#########
# rbind all samples count tables and define good cells and background cells based on UMI threshold
cell_counts <- data.frame(cell=character(),
                          total=numeric(), 
                          total_deduplicated=numeric(),
                          total_deduplicated_peaks=numeric(), 
                          total_deduplicated_tss=numeric(), 
                          stringsAsFactors=FALSE) 
for (s in samples){
  message("loading", s)
  df = read.table(s, head = T)
  cell_counts = rbind(cell_counts, df)
}

# convert barcodes to wells for merging with HashTable
cell_counts$wellID = sapply(cell_counts$cell, BCtoWellIDs)
BCtoID = dplyr::select(cell_counts, cell, wellID) 
colnames(BCtoID)<- c("barcode", "Cell")

# headUMI cutoff is set in pipeline based on knee
cutoff = 336
# get background cells (below cutoff)
background_cells = filter(cell_counts, total_deduplicated < cutoff)
background_cells = as.character(background_cells$wellID)
# get good cells (above cutoff)
real_test_cell_df = filter(cell_counts, total_deduplicated >= cutoff)
test_cells = as.character(real_test_cell_df$wellID)

#########
# Read in hash Oligo count table and reformat Cell IDs 
hashTable = read.csv(paste0(basepath,"pipeline_hill/analysis_hash/hashRDS/hashTable.out"), sep = "\t", header = FALSE)
colnames(hashTable) = c("SampleName", "Cell", "Oligo", "Axis", "Count")
## Need to remap cell IDs based on different RT well layout
hashTable$PCRwell=substr(hashTable$Cell,1,3)
hashTable$PCRplate=substr(hashTable$Cell,5,7)
hashTable$RTwell =stringr::str_split_fixed(string =  hashTable$Cell, pattern = "_", n = 8)[,4]
remaptoNext = function(x){
  w = c(0,12,24,36,48,60,72,84)
  TruseqRTwell = 0:95
  nextRTwell = c(w, w+1, w+2, w+3, w+4, w+5, w+6, w+7, w+8, w+9, w+10, w+11)
  TruWell = match(x, TruseqRTwell)
  NextWell = nextRTwell[TruWell]
  return(NextWell)
}
hashTable_f <- mutate(hashTable, RTwell_N = remaptoNext(RTwell)) %>%
  mutate(Cell_N = paste(PCRwell,PCRplate,"RT",RTwell_N, sep = "_")) %>%
    select(SampleName, Cell_N, Count, Oligo, Axis)
colnames(hashTable_f) = c("SampleName", "Cell", "UMI", "Oligo", "Axis")
hashTable_f$Count = hashTable_f$UMI

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
background_cells_curr = background_cells[background_cells %in% hashTable_f$Cell]
background_sample_hashes = hashTable_f %>% 
  filter(Cell %in% background_cells_curr, Axis == 1)
test_cell_hashes = subset(hashTable_f, Cell %in% test_cells & Axis == 1)
hash_df_tmp = assign_hash_labels(test_cell_hashes, background_sample_hashes, downsample_rate=downsample_rate)
hash_df_all = rbind(hash_df_all, hash_df_tmp)

##################
# map hash labels to samples
hash_df_all$RTwell = stringr::str_split_fixed(hash_df_all$Cell, pattern = "_" ,n = 4)[,4]
hash_df_all$Plate = stringr::str_split_fixed(hash_df_all$top_oligo, pattern = "_" ,n = 2)[,1]
hash_df_all$hashWell = stringr::str_split_fixed(hash_df_all$top_oligo, pattern = "_" ,n = 2)[,2]
hash_df_all$hashRow = substr(hash_df_all$hashWell,1,1)
hash_df_all$hashCol = substr(hash_df_all$hashWell,2,3)

RTwell_coding <- sapply(hash_df_all$RTwell, function(x){
  if(x %in% as.character(0:87)) return("SciChem")
  if(x %in% as.character(88:91)) return("PreHashBarnyard")
  if(x %in% as.character(92:95)) return("PostHashBarnyard")
  return("NA")
})
platecode = function(hashplate, column, row){
  if (hashplate == "plate19"){
    if(row %in% c("A", "B", "C", "D", "E", "F", "G", "H")){
      if(column %in% as.character(1:3)) return("Nutlin")
      if(column %in% as.character(4:6)) return("Dex")
      if(column %in% as.character(7:9)) return("SAHA")
      if(column %in% as.character(10:12)) return("BMS")
    }
  }
  if(hashplate == "Plate1"){
    if (column == as.character(1)){
      if (row == "A") return("PREmixed")
      if (row == "B") return("mouse")
      if (row == "C") return("human")
    }
    #else return("NA")
  }
}
getCellType = function(hashplate, column, row){
  if (hashplate == "plate19") return("A549")
  if(hashplate == "Plate1"){
    if (column == as.character(1)){
      if (row == "A") return("Barnyard")
      if (row == "B") return("NIH-3T3")
      if (row == "C") return("HEK293T")
    }
    #else return("NA")
  }
}
dose_code =  function(hashplate, column, row){
  if (hashplate == "plate19"){
    if(row %in% c("A", "B", "C", "D", "E", "F", "G", "H")){
      if(column %in% as.character(1:3)){# return("Nutlin")
        if (row == "A") return(250)
        if (row == "B") return(125)
        if (row == "C") return(25)
        if (row == "D") return(12.5)
        if (row == "E") return(2.5)
        if (row == "F") return(1.25)
        if (row == "G") return(0.25)
        if (row == "H") return(0)
      }
      if(column %in% as.character(4:6)){# return("Dex")
        if (row == "A") return(50)
        if (row == "B") return(25)
        if (row == "C") return(5)
        if (row == "D") return(2.5)
        if (row == "E") return(0.5)
        if (row == "F") return(0.25)
        if (row == "G") return(0.05)
        if (row == "H") return(0)
      }
      if(column %in% as.character(7:12)){ # both BMS and SAHA
        if (row == "A") return(1000)
        if (row == "B") return(500)
        if (row == "C") return(100)
        if (row == "D") return(50)
        if (row == "E") return(10)
        if (row == "F") return(5)
        if (row == "G") return(1)
        if (row == "H") return(0)
      }
    }
    else return(NA)
  }
}
Relative_dose_code =  function(hashplate, row){
  if (hashplate == "plate19"){
    if (row == "A") return(100)
    if (row == "B") return(50)
    if (row == "C") return(10)
    if (row == "D") return(5)
    if (row == "E") return(1)
    if (row == "F") return(0.5)
    if (row == "G") return(0.1)
    if (row == "H") return(0)
  }
  else return(NA)
}
isDoublet = function(ratio, qval, FDR = 0.05, cutoff = 2){
  if (ratio > cutoff & qval < FDR) return("singlet")
  else(return("doublet"))
}
getReplicate = function(hashplate, column){
  if (hashplate == "plate19"){
    if(column %in% as.character(c(1,4,7,10))) return(1)
    if(column %in% as.character(c(2,5,8,11))) return(2)
    if(column %in% as.character(c(3,6,9,12))) return(3)
  }
  else return("NA")
}
getCulturePlate = function(hashplate){
  if (hashplate == "plate19") return(1)
  else return(2)
}
getWellOligo = function(hashplate){ ## Note:  This info is for when wells recieve more than one oligo (i.e. like axis)
  if (hashplate == "plate19") return(1)
  else return(2)
}
getVehicle = function(hashplate, column){
  if (hashplate == "plate19"){
    if(column %in% as.character(c(1:3,7:12))) return("DMSO")
    if(column %in% as.character(4:6)) return("ethanol")
  }
  else return("NA")
}
isVehicle =  function(hashplate, row){
  if (hashplate == "plate19"){
    if(row == "H") return(TRUE)
    else return(FALSE)
  }
  else return(FALSE)
}

hash_df_all$hash_coding <- mapply(platecode, hash_df_all$Plate, hash_df_all$hashCol, hash_df_all$hashRow)
hash_df_all$RTwell_coding = RTwell_coding
hash_df_all$dose_code = mapply(dose_code, hash_df_all$Plate, hash_df_all$hashCol, hash_df_all$hashRow)
hash_df_all$Relative_dose = mapply(Relative_dose_code, hash_df_all$Plate, hash_df_all$hashRow)
hash_df_all$doublet = mapply(isDoublet, hash_df_all$top_to_second_best_ratio, hash_df_all$qval,  FDR = 0.05, cutoff = 2)
hash_df_all$replicate = mapply(getReplicate, hash_df_all$Plate, hash_df_all$hashCol)
hash_df_all$cell_type <- mapply(getCellType, hash_df_all$Plate, hash_df_all$hashCol, hash_df_all$hashRow)
hash_df_all$culture_plate <- mapply(getCulturePlate, hash_df_all$Plate)
hash_df_all$well_oligo <- mapply(getWellOligo, hash_df_all$Plate)
hash_df_all$vehicle_type = mapply(getVehicle, hash_df_all$Plate, hash_df_all$hashCol)
hash_df_all$vehicle = mapply(isVehicle, hash_df_all$Plate, hash_df_all$hashRow)
hash_df_all$Cell = as.character(hash_df_all$Cell)
hash_df_all = inner_join(hash_df_all, BCtoID, by = "Cell")

# make characters
hash_df_all$sample = hash_df_all$RTwell_coding
hash_df_all$treatment = as.character(hash_df_all$hash_coding)
hash_df_all$dose = as.numeric(as.character(hash_df_all$dose_code))
hash_df_all$Relative_dose = as.numeric(as.character(hash_df_all$Relative_dose))
hash_df_all$doublet = as.character(hash_df_all$doublet)
hash_df_all$replicate = as.character(hash_df_all$replicate)
hash_df_all$cell_type = as.character(hash_df_all$cell_type)
hash_df_all$culture_plate = as.numeric(hash_df_all$culture_plate)
hash_df_all$well_oligo = as.numeric(hash_df_all$well_oligo)
hash_df_all$vehicle_type = as.character(hash_df_all$vehicle_type)

hashOut <- select(hash_df_all, barcode, Cell, hash_umis, pval,
                  qval, top_to_second_best_ratio, cell_type, 
                  sample, replicate, culture_plate, well_oligo,
                  treatment, dose, Relative_dose, vehicle, 
                  vehicle_type, doublet)

colnames(hashOut)<- c("cell", "wellID", "hash_umis", "pval", 
                      "qval", "top_to_second_best_ratio", "cell_type", 
                      "sample", "replicate", "culture_plate",
                      "well_oligo", "treatment", "dose", "Relative_dose",
                      "vehicle", "vehicle_type", "doublet")

write.table(hashOut, file= paste0(out_dir, "hashCellAssignments.txt"), 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)




