basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/191218_3Level_scichem_MLR_combinedRuns/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB4"))
setwd(paste0(out_dir, "results/NB4"))

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
  library(RColorBrewer)
  library(ggpubr)
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 

###################################### 
# Using Replicate information, plot cell type
# proportions for each condition. 
# want to look for significant changes 
# (i.e. activated Tcells after bead stimulation)
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

# make box plots of cell proportion faceted by responder. 
# prepare separate plots for each cell type 

# omit stim_alone condition from plots
ccs = filter(cell_counts_summary, Responder != "stimAlone")

celltypes = ccs %>% 
  distinct(cellType_broad)

for (ct in celltypes$cellType_broad){
  pdf(file = paste0("BoxPlot_PercentRecoveredCellType_byStimulator_", ct,"2.pdf"), 
      width = 3, height = 3)
  print(  
    ccs %>% 
      filter(cellType_broad == ct) %>% 
      ggplot(aes(x =  Stimulator, y = percent_cellType)) +
      geom_boxplot(fill = "gray") +
      facet_wrap("~Responder") +
      theme_bw() +
      #scale_fill_brewer(palette = "RdYlBu") +
      xlab("Stimulator") +
      ylab("cells") +
      ggtitle(paste0("Percent recovered ", ct)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
      stat_compare_means(label = "p.signif", method = "t.test",  
                         ref.group = "noStim", hide.ns = TRUE, 
                         label.y.npc = 0.85, col = "red")  # Pairwise comparison against reference
  )
  dev.off()
}
