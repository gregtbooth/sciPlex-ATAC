basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/211103_3Level_scichem_MLR_Drugs/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB16"))
setwd(paste0(out_dir, "results/NB16"))


suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gt)
})

set.seed(2017) # ensures reproducibility of previous random number generation

prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered_annotated"), 
                                              showLogo = FALSE) 

# combined replicates meta data
meta_data = prj@cellColData %>% 
  data.frame() %>% 
  mutate(MLR_drug_sample = paste0(Responder, "_", Stimulator, "_", Drug, "_", Relative_Dose, "_", Replicate)) %>% 
  group_by(MLR_drug_sample) %>% 
  summarise(
    n_cells = n(),
    med_frags = median(nFrags), 
    med_TSSenr = median(TSSEnrichment), 
    med_FRIP = median(FRIP), 
    med_hash = median(hash_umis),
    med_hashEnr = median(top_to_second_best_ratio)
  )

write.table(format(data.frame(meta_data), digits=2), file = "MLRDrugMetaData.csv", sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)

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


##############################
# make table of drugs and doses 
cols = data.frame("drug", "vehicle", "dose_1", "dose_10", "dose_100")
BMS = c("BMS", "DMSO", "0.5uM", "5uM", "50uM")
CycA = c("CycA", "DMSO", "0.25uM", "2.5uM", "25uM")
CycA_Rapa = c("CycA, Rapa", "DMSO", "0.125uM, 0.25nM", "1.25uM, 2.5nM", "12.5uM, 25nM")
Dex = c("Dex", "Ethanol",  "0.25uM", "2.5uM", "25uM")
Meth = c("Meth", "DMSO", "0.1uM", "1uM", "10uM")
Meth_CycA = c(" Meth, CycA", "DMSO", "0.05uM, 0.125uM", "0.5uM, 1.25uM", "5uM, 12.5uM")
Rapa = c("Rapa", "DMSO", "0.5uM", "5uM", "50uM")
SAHA = c("SAHA", "DMSO", "5uM", "50uM", "500uM")


drug_df = data.frame(rbind(BMS, CycA, CycA_Rapa, Dex, Meth, Meth_CycA, Rapa, SAHA))
colnames(drug_df) = cols

ddf = drug_df %>% tibble::rownames_to_column(var = "Drug") %>% 
  select(-drug)

gtable = tibble(ddf) %>%  gt() 
  #fmt_number(columns = vars(med_TSSenr, med_FRIP, med_hashEnr), decimals = 2) 
gtsave(gtable, filename = "DoseData.pdf")






