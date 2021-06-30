basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB15"))
setwd(paste0(out_dir, "results/NB15"))

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gt)
})

set.seed(2017) # ensures reproducibility of previous random number generation

prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 

# combined replicates meta data
meta_data = prj@cellColData %>% 
  data.frame() %>% 
  mutate(MLR_sample = paste0(Responder, "_", Stimulator)) %>% 
  group_by(MLR_sample) %>% 
  summarise(
    n_cells = n(),
    med_frags = median(nFrags), 
    med_TSSenr = median(TSSEnrichment), 
    med_FRIP = median(FRIP), 
    med_hash = median(hash_umis),
    med_hashEnr = median(top_to_second_best_ratio)
  )

write.table(format(data.frame(meta_data), digits=2), file = "MLRmetaData.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# save formatted table as image:
gtable = meta_data %>% gt() %>% 
  fmt_number(columns = vars(med_TSSenr, med_FRIP, med_hashEnr), decimals = 2) 
gtsave(gtable, filename = "MLRmetaData.pdf")

# separate replicates meta data
meta_data_reps = prj@cellColData %>% 
  data.frame() %>% 
  mutate(MLR_sample_rep = paste0(Responder, "_", Stimulator, "_rep", Replicate)) %>% 
  group_by(MLR_sample_rep) %>% 
  summarise(
    n_cells = n(),
    med_frags = median(nFrags), 
    med_TSSenr = median(TSSEnrichment), 
    med_FRIP = median(FRIP), 
    med_hash = median(hash_umis),
    med_hashEnr = median(top_to_second_best_ratio)
  )

write.table(format(data.frame(meta_data_reps), digits=2), file = "MLRmetaData_replicates.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# save formatted table as image:
gtable_reps = meta_data_reps %>% gt() %>% 
  fmt_number(columns = vars(med_TSSenr, med_FRIP, med_hashEnr), decimals = 2) 
gtsave(gtable_reps, filename = "MLRmetaData_reps.pdf")
