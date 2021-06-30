library(devtools)
library(ggplot2)
library(ggrastr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggridges)
library(UpSetR)
source("~/scripts/GB_src/io_functions.R")
source("~/scripts/GB_src/chiSq_2lvl_GB.R")
source("~/scripts/GB_src/ParseBCtoWells.R")
basepath = "github/"
matpath = paste0(basepath, "pipeline_hill/analysis_barnyard/make_matrices/")
out_dir = paste0(basepath,"analysis/results/NB1/")
dir.create(out_dir)
setwd(out_dir)


################################
# merge sample window matrices (only need to do once)

samples = c("polyA", "enhanced", "polyA_LNA", "enhanced_LNA")

matrices = list()
for(s in samples){
  matrices[[s]] = load_mtx_file(paste0(matpath, s, ".", "window_matrix.mtx.gz"))
}

# merge sample bin matrices and write as new matrix
if (length(matrices) < 2) {
  merged_data = matrices[[1]]
} else if (length(matrices) == 2) {
  merged_data = Seurat:::RowMergeSparseMatrices(matrices[[1]], matrices[[2]])
} else {
  merged_data = matrices[[1]]
  for (i in 2:length(matrices)) {
    merged_data <<- Seurat:::RowMergeSparseMatrices(merged_data, matrices[[i]])
  }
}

write_mtx_file(merged_data, file=paste0(matpath, "merged_data.window_matrix.mtx.gz"))

################################

DOUBLET_PERCENTAGE_THRESHOLD = 0.9

message('-> loading window matrix to calculate doublet stats...')
sample_window_matrix = load_mtx_file(paste0(matpath, "merged_data.window_matrix.mtx.gz"))

hg19_counts = Matrix::colSums(sample_window_matrix[grepl('hg19', rownames(sample_window_matrix)),])
mm9_counts = Matrix::colSums(sample_window_matrix[grepl('mm9', rownames(sample_window_matrix)),])

if (sum(hg19_counts) == 0) {
  stop('No hg19 counts found. Chromosomes must start with hg19 in name (e.g. hg19chr19) for human and mm9 for mouse.')
}
if (sum(mm9_counts) == 0) {
  stop('No mm9 counts found. Chromosomes must start with mm9 in name (e.g. mm9chr19) for mouse and hg19 for human.')
}

by_df = data.frame(cell=names(mm9_counts), hg19=hg19_counts, mm9=mm9_counts) %>% 
  mutate(count = hg19 + mm9)


#########
#define good cells and background cells based on UMI threshold

# UMI cutoff is set in pipeline based on knee
cutoff = 500

# get background cells (below cutoff)
background_cells = dplyr::filter(by_df, count < cutoff)
background_cells = as.character(background_cells$cell)
# get good cells (above cutoff)
real_test_cell_df = dplyr::filter(by_df, count >= cutoff)
test_cells = as.character(real_test_cell_df$cell)

#########
# Read in hash Oligo count table and reformat Cell IDs 
hashTable = read.csv(paste0(basepath,"pipeline_hill/hashRDS/hashTable.out"), sep = "\t", header = FALSE)
colnames(hashTable) = c("SampleName", "Cell", "Oligo", "Axis", "Count")

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
#hash_df_all$RTwell = stringr::str_split_fixed(hash_df_all$Cell, pattern = "_" ,n = 4)[,4]
hash_df_all_f$Plate = stringr::str_split_fixed(hash_df_all_f$top_oligo, pattern = "_" ,n = 2)[,1]
hash_df_all_f$hashWell = stringr::str_split_fixed(hash_df_all_f$top_oligo, pattern = "_" ,n = 2)[,2]
hash_df_all_f$hashRow = substr(hash_df_all_f$hashWell,1,1)
hash_df_all_f$hashCol = substr(hash_df_all_f$hashWell,2,3)
hash_df_all_f = cbind(hash_df_all_f, N7well = sapply(hash_df_all_f$Cell, function(z) get_N7_well(z)))

Cellcode = function(hashplate, column, row){
  if (hashplate == "plate19"){
    ifelse(row %in% c("A", "B", "C"), return("A549"), return("NIH-3T3"))
    }
  else {
    ifelse(row %in% c("A", "B", "C"), return("A549"), return("NIH-3T3"))
  }
}

capture_code = function(hashplate, column, row){
  ifelse (hashplate == "plate19", return("polyA_hash"), return("enhanced_hash"))
}

isDoublet = function(ratio, qval, FDR = 0.05, cutoff = 2){
  if (ratio > cutoff & qval < FDR) return("singlet")
  else(return("doublet"))
}

captureCode <- sapply(hash_df_all_f$N7well, function(x){
  colNum = as.numeric(stringr::str_sub(x, 9,10))
  ifelse(colNum <= 6, return("polydT"), return("Enhanced"))
})

OligoCode <- sapply(hash_df_all_f$N7well, function(x){
  rowLet = stringr::str_sub(x, 8,8)
  ifelse(rowLet %in% c("A", "B", "C", "D"), return("standard"), return("LNA"))
})


hash_df_all_f$Cell_coding <- mapply(Cellcode, hash_df_all_f$Plate, hash_df_all_f$hashCol, hash_df_all_f$hashRow)
hash_df_all_f$hash_capture_seq <- mapply(capture_code, hash_df_all_f$Plate, hash_df_all_f$hashCol, hash_df_all_f$hashRow)
hash_df_all_f$captureCode = captureCode
hash_df_all_f$OligoCode = OligoCode
hash_df_all_f$doublet = mapply(isDoublet, hash_df_all_f$top_to_second_best_ratio, hash_df_all_f$qval,  FDR = 0.05, cutoff = 5)

# make characters
hash_df_all_f$celltype = as.character(hash_df_all_f$Cell_coding)
hash_df_all_f$lysis = as.character(hash_df_all_f$hash_capture_seq)
hash_df_all_f$doublet = as.character(hash_df_all_f$doublet)

# add mouse and human frag counts to table
hash_df_all_f = inner_join(hash_df_all_f, by_df, by = c("Cell" = "cell"))

hashOut <- select(hash_df_all_f, Cell, hash_umis, pval, qval,
                  top_to_second_best_ratio, celltype, lysis, 
                  OligoCode, captureCode, doublet, 
                  hg19_frags = hg19, mm9_frags = mm9)

write.table(hashOut, file= paste0(out_dir, "hashCellAssignments.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


