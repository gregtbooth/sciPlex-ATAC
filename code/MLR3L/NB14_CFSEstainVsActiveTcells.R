basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB14"))
setwd(paste0(out_dir, "results/NB14"))

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggrastr)
  library(ggpubr)
})

set.seed(2017) # ensures reproducibility of previous random number generation

prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 

# load and reformat the FACS data
sample_key = c(A1 = "respA_noStim", A2 = "respA_stimA", A3 = "respA_stimB", A4 = "respA_stimC", A5 = "respA_stimD", A6 = "respA_Bead", A7 = "stimAlone_stimA", 
               B1 = "respB_noStim", B2 = "respB_stimB", B3 = "respB_stimA", B4 = "respB_stimC", B5 = "respB_stimD", B6 = "respB_Bead", B7 = "stimAlone_stimB",
               C1 = "respC_noStim", C2 = "respC_stimC", C3 = "respC_stimA", C4 = "respC_stimB", C5 = "respC_stimD", C6 = "respC_Bead", C7 = "stimAlone_stimC",
               D1 = "respD_noStim", D2 = "respD_stimD", D3 = "respD_stimA", D4 = "respD_stimB", D5 = "respD_stimC", D6 = "respD_Bead", D7 = "stimAlone_stimD")

FACS_dat = read.table(file = paste0(out_dir, "external_resources/FACS/08-Apr-2021_FACS_cellProliferationStatsTable.txt"), 
                      head = TRUE, fill = TRUE) %>% 
  mutate(Cells = ifelse(is.na(Cells), Statistic, Cells), # corrects the filling of columns
         percent_parent = ifelse(Statistic == Cells, 100, Statistic), 
         sample = stringr::str_split_fixed(string = Name, pattern = "/", n = 3)[,1], 
         sample_short = stringr::str_split_fixed(string = sample, pattern = "_", n = 3)[,2], 
         sample_recode = recode(sample_short, !!!sample_key))
FACS_dat$population = rep(c("All", "Lymphocytes", "Proliferating"), (nrow(FACS_dat)/3))
FACS_dat = select(FACS_dat, sample_recode, population, Cells, percent_parent)

FACS_counts_spread = FACS_dat %>% 
  dcast(sample_recode ~ population, value.var = "Cells") %>% 
  mutate(percent_pro = (Proliferating/All)*100, 
         percent_pro_lymph = (Proliferating/Lymphocytes)*100)

####
# collect active/proliferating T-cell data from sciPlexATAC
cell_counts_condition = data.frame(prj@cellColData) %>% 
  group_by(Replicate, Responder, Stimulator) %>% 
  summarise(total_cells = n()) %>% 
  mutate(condition = paste(Responder, Stimulator, Replicate, sep = "_"))

cell_counts_summary = data.frame(prj@cellColData) %>% 
  group_by(Replicate, Responder, Stimulator, cellType_broad) %>% 
  summarise(n_cells = n()) %>% 
  mutate(condition = paste(Responder, Stimulator, Replicate, sep = "_")) %>% 
  left_join(cell_counts_condition, by = "condition") %>% 
  select(
    Replicate = Replicate.x,
    Responder = Responder.x, 
    Stimulator = Stimulator.x, 
    cellType_broad, 
    n_cells, 
    total_cells) %>% 
  mutate(percent_cellType = (n_cells/total_cells)*100) 

Tcell_summary = cell_counts_summary %>% 
  filter(cellType_broad == "Activated_Tcell") %>% 
  mutate(resp_stim = paste0(Responder, "_", Stimulator)) %>% 
  group_by(resp_stim) %>% # sum replicates
  summarize(n_actTcells = sum(n_cells), 
            total_cells = sum(total_cells)) %>% 
  mutate(percent_actT = (n_actTcells/total_cells)*100)

###########################
# look at corr between FACS and ATAC T-cell recovery 

Tcell_dat = left_join(Tcell_summary, FACS_counts_spread, by = c("resp_stim" = "sample_recode")) %>% 
  mutate(resp = stringr::str_split_fixed(resp_stim, "_", 2)[,1])

# exclude stim alone:
Tcell_dat %>% filter(!grepl("stimAlone", resp_stim)) %>%
ggscatter(x = "percent_actT", y = "percent_pro",
          color = "black", size = 2, # Points color, shape and size
          xlab = "% Active T-cells recovered (ATAC)",
          ylab = "% Proliferating T-cells (FACS)",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"))+
  geom_abline(colour = "red")

ggsave(file = "scatterplot_PercentActTcell_FACSvsATAC_wBead.pdf", height = 3, width = 3)


# exclude bead and stim alone:
Tcell_dat %>% filter(!grepl("Bead", resp_stim), !grepl("stimAlone", resp_stim)) %>% 
  ggscatter(x = "percent_actT", y = "percent_pro",
            color = "black", size = 2, # Points color, shape and size
            xlab = "% Active T-cells recovered (ATAC)",
            ylab = "% Proliferating T-cells (FACS)",
            add = "reg.line",  # Add regressin line
            add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
            conf.int = TRUE, # Add confidence interval
            cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
            cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"))+
  geom_abline(colour = "red")


ggsave(file = "scatterplot_PercentActTcell_FACSvsATAC_noBead_orStimAlone.pdf", height = 3, width = 3)


