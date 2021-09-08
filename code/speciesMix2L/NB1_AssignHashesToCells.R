library(Matrix)
library(argparse)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrastr)
source("~/scripts/GB_src/io_functions.R")
source("~/scripts/GB_src/chiSq_2lvl_GB.R")
source("~/scripts/GB_src/ConvertBCs_ToWells_V3.R")
basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190421_barnyard7/"
matpath = paste0(basepath, "pipeline_hill/analysis_barnyard/make_matrices/")
out_dir = paste0(basepath,"analysis/results/NB1/")
dir.create(out_dir)
setwd(out_dir)

################################
# merge sample window matrices (only need to do once)

samples = c("CLB_NEB", "CLB_SSRT", "OMNI_NEB", "OMNI_SSRT", "PS2_NEB", "PS2_SSRT")

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

# tabulate Frag counts/cell from each sepecies
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

# convert barcodes to wells for merging with HashTable
by_df$wellID = sapply(by_df$cell, BCtoWellIDs)

BCtoID = dplyr::select(by_df, cell, wellID) 
colnames(BCtoID)<- c("barcode", "Cell")


# generate quick knee plot to identify background cells 

kp = by_df %>% 
  mutate(count.rank = min_rank(-count)) %>% 
  arrange(-count) %>% 
  select(count, count.rank) %>% 
  distinct() 

knee = ggplot(kp, aes(x = count.rank, y = count)) +
  geom_line(size = 0.8) +
  scale_x_log10(limits = c(10, NA),
                breaks = c(10, 20, 40, 100, 200, 400, 1000, 2000, 4000, 8000, 16000,32000,64000,100000,250000,500000,750000, 1000000)) +
  scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000)) +
  xlab("# of barcodes") +
  ylab("UMI count threshold") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste(out_dir, "/allcells_knee.pdf", sep = ""),
       plot = knee, units = "in", width = 3, height = 3)


# headUMI cutoff is set in based on knee
cutoff = 500
# get background cells (below cutoff)
background_cells = dplyr::filter(by_df, count < cutoff)
background_cells = as.character(background_cells$wellID)
# get good cells (above cutoff)
real_test_cell_df = dplyr::filter(by_df, count >= cutoff)
test_cells = as.character(real_test_cell_df$wellID)

#########
# Read in hash Oligo count table and reformat Cell IDs 
hashTable = read.csv(paste0(basepath,"hash/hashRDS/hashTable.out"), sep = "\t", header = FALSE)
colnames(hashTable) = c("SampleName", "Cell", "Oligo", "Axis", "Count")

## Need to remap cell IDs based on different RT well layout (Nextera vs Truseq)
hashTable$PCRwell =stringr::str_split_fixed(string =  hashTable$Cell, pattern = "_", n = 4)[,1]
hashTable$PCRplate =stringr::str_split_fixed(string =  hashTable$Cell, pattern = "_", n = 4)[,2]
hashTable$RTwell =stringr::str_split_fixed(string =  hashTable$Cell, pattern = "_", n = 4)[,4]
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
background_sample_hashes = hashTable_f %>% filter(Cell %in% background_cells_curr, Axis == 1)
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
hash_df_all$Tn5well = stringr::str_split_fixed(string =  hash_df_all$Cell, pattern = "_", n = 4)[,4]

hash_coding <- sapply(hash_df_all$hashCol, function(x){
  if(x == 1) return("A549_CLB")
  if(x == 2) return("A549_CLB")
  if(x == 3) return("3T3_CLB")
  if(x == 4) return("3T3_CLB")
  if(x == 5) return("A549_OMNI")
  if(x == 6) return("A549_OMNI")
  if(x == 7) return("3T3_OMNI")
  if(x == 8) return("3T3_OMNI")
  if(x == 9) return("A549_PS2")
  if(x == 10) return("A549_PS2")
  if(x == 11) return("3T3_PS2")
  if(x == 12) return("3T3_PS2")
})

isDoublet = function(ratio, qval, FDR = 0.05, cutoff = 2){
  if (ratio > cutoff & qval < FDR) return("singlet")
  else(return("doublet"))
}


Tn5_WellTreat_coding <- sapply(hash_df_all$Tn5well, function(x){
  w= c(0:7) # well numbers for first 3 cols
  if(x %in% as.character(w)) return("CLB_SSRT")
  if(x %in% as.character(w+8)) return("CLB_SSRT")
  if(x %in% as.character(w+16)) return("CLB_NEB")
  if(x %in% as.character(w+24)) return("CLB_NEB")
  if(x %in% as.character(w+32)) return("OMNI_SSRT")
  if(x %in% as.character(w+40)) return("OMNI_SSRT")
  if(x %in% as.character(w+48)) return("OMNI_NEB")
  if(x %in% as.character(w+56)) return("OMNI_NEB")
  if(x %in% as.character(w+64)) return("PS2_SSRT")
  if(x %in% as.character(w+72)) return("PS2_SSRT")
  if(x %in% as.character(w+80)) return("PS2_NEB")
  if(x %in% as.character(w+88)) return("PS2_NEB")
  return("NA")
})

hash_df_all$hash_coding <- hash_coding
hash_df_all$cell_type <-  stringr::str_split_fixed(hash_df_all$hash_coding, pattern = "_" ,n = 2)[,1]
hash_df_all$doublet = mapply(isDoublet, hash_df_all$top_to_second_best_ratio, hash_df_all$qval,  FDR = 0.05, cutoff = 2)
hash_df_all$doublet = as.character(hash_df_all$doublet)
hash_df_all$Tn5_WellTreat_coding = Tn5_WellTreat_coding
hash_df_all = inner_join(hash_df_all, by_df, by = c("Cell" = "wellID"))

# write hash assignment table to file
hashOut <- select(hash_df_all, cell, wellID = Cell, hash_umis, well_oligo = top_oligo, pval, qval,
                  top_to_second_best_ratio, hash_coding, cell_type, doublet, Tn5_WellTreat_coding,  hg19_frags = hg19, mm9_frags = mm9)

write.table(hashOut, file= paste0(out_dir, "hashCellAssignments.txt"),
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)


# create raw count matrix (sparse) files
bin_by_cell_matrix_sp = Matrix(sample_window_matrix, sparse = TRUE)
bin_by_cell_matrix_rows = row.names(sample_window_matrix)
bin_by_cell_matrix_cols = colnames(sample_window_matrix)

writeMM(bin_by_cell_matrix_sp,file='BY2_bin_matrix.mtx.txt')
write.table(bin_by_cell_matrix_rows, file='BY2_bin_matrix.rows.txt', sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(bin_by_cell_matrix_cols, file='BY2_bin_matrix.cols.txt', sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)


system('gzip BY2_bin_matrix.mtx.txt')
system('gzip BY2_bin_matrix.rows.txt')
system('gzip BY2_bin_matrix.cols.txt')


# Test loading matrix
#pm = readMM("BY2_bin_matrix.mtx.txt.gz")
#rn = read.table("BY2_bin_matrix.rows.txt.gz", header = FALSE, comment.char = "")
#cn = read.table("BY2_bin_matrix.cols.txt.gz", header = FALSE, comment.char = "")
#rownames(pm) <- rn$V1
#colnames(pm) <- cn$V1
#pm[1:10,1:10]


