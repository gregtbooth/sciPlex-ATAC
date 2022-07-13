suppressPackageStartupMessages(c(
library(devtools),
library(ggplot2),
library(dplyr),
library(tidyr),
library(reshape2),
library(ggridges),
library(UpSetR))
)

binpath ="/bin/"
basepath = "/home/"
source(paste0(binpath,"chiSq_2lvl_GB.R"))
out_dir = paste0(basepath,"analysis/hash/")
clean_hashList = read.table(file = paste0(basepath, "pipeline/hash_pipeline/HashSampleSheet.txt"), head = FALSE) 
colnames(clean_hashList) = c("Oligo", "ID", "Axis")

samples = c("MLR_5um_homeTn5", "MLR_DiagenodeTn5", "MLR_homeTn5", "MLR_IlluminaTn5")
#########
# rbind all samples count tables and define good cells and background cells based on UMI threshold
cell_counts <- data.frame(cell=character(),
                          total=numeric(), 
                          total_deduplicated=numeric(),
                          total_deduplicated_peaks=numeric(), 
                          total_deduplicated_tss=numeric(), 
                          stringsAsFactors=FALSE) 
for (s in samples){
  df = read.table(file = paste0(basepath, "pipeline_output/analysis_MLR_hg19/count_report/", s, ".count_report.txt"), head = T)
  cell_counts = rbind(cell_counts, df)
}

# UMI cutoff is set in pipeline based on knee
cutoff = 142
# get background cells (below cutoff)
background_cells = filter(cell_counts, total_deduplicated < cutoff)
background_cells = as.character(background_cells[,1])
# get good cells (above cutoff)
real_test_cell_df = filter(cell_counts, total_deduplicated >= cutoff)
test_cells = as.character(real_test_cell_df[,1])

#########
# Read in hash Oligo count table and reformat Cell IDs 
hashTable = read.csv(paste0(basepath,"pipeline_output/hashRDS/hashTable.out"), sep = "\t", header = FALSE)
colnames(hashTable) = c("SampleName", "Cell", "Oligo", "Axis", "Count")

# filter out hash not used in experiment
hashTable = filter(hashTable, Oligo %in% clean_hashList$Oligo)

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

isDoublet = function(ratio, qval, FDR = 0.05, cutoff = 2){
  if (ratio > cutoff & qval < FDR) return("singlet")
  else(return("doublet"))
}

hash_df_all_f$doublet = mapply(isDoublet, hash_df_all_f$top_to_second_best_ratio, hash_df_all_f$qval,  FDR = 0.05, cutoff = 2)
# make characters
hash_df_all_f$doublet = as.character(hash_df_all_f$doublet)

hashOut <- select(hash_df_all_f, Cell, hash_umis, pval, qval, top_to_second_best_ratio, top_oligo, doublet)
write.table(hashOut, file= paste0(out_dir, "hashCellAssignments.txt"), 
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# quick plot of hash enrichments
p = qplot(log10(hashOut$top_to_second_best_ratio))
ggsave(p, file = paste0(out_dir,"hashEnrichment.png"), height = 4, width = 4)



