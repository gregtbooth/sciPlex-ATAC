basepath = "/home/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB6"))
setwd(paste0(out_dir, "results/NB6"))

suppressPackageStartupMessages({
  library(ArchR)
 # library(Seurat)
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

prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered"), 
                       showLogo = FALSE) 

####################################################
# count the number and proportion of cells 
# within each cluster from each condition 
cdat = data.frame(prj@cellColData) %>% 
  tibble::rownames_to_column(var = "cellNames") 

# count treatment cells in clusters
cell_counts_summary = group_by(cdat, Drug, Relative_Dose, Stimulator, clusters_mon3) %>% 
  summarise(n_cells = n()) 

# stacked bar of proportion stimulators per cluster
cell_counts_summary %>%
  dplyr::select(Drug, Relative_Dose,Stimulator, clusters_mon3, n_cells) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = clusters_mon3, y =n_cells, fill = Relative_Dose), 
           color = "black", size = .25, stat = "identity", 
           position = "fill") +
  facet_wrap(~Drug+Stimulator) +
  #scale_fill_manual("Dose (uM)", values = clr_values) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  #scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Cluster") +
  ylab("Proportion") +
  ggsave("stackedbarplot_stimulator_DrugDose_Perclusters.pdf",
         width = 6, height = 4.5, unit = "in", dpi = 1200)

# stacked bar of proportion cluster per stimulator
cell_counts_summary %>% 
  dplyr::select(Drug, Relative_Dose, Stimulator, clusters_mon3, n_cells) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = Relative_Dose, y =n_cells, fill = clusters_mon3), 
           color = "black", size = .25, stat = "identity", 
           position = "fill") +
  facet_wrap(~Drug+Stimulator) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  xlab("Stimulator") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("stackedbarplot_clustersPerDrugDose.pdf",
       width = 6, height = 4.5, unit = "in", dpi = 1200)


# stacked bar of proportion cluster per stimulator 
# Vehical treated only 
cell_counts_summary %>% 
  filter(Relative_Dose == 0) %>% 
  group_by(clusters_mon3, Stimulator) %>% summarise(n_cells_new = sum(n_cells)) %>% 
  #dplyr::select(Drug, Relative_Dose, Stimulator, clusters_mon3, n_cells) %>%
  dplyr::select(Stimulator, clusters_mon3, n_cells_new) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = Stimulator, y =n_cells_new, fill = clusters_mon3), 
           color = "black", size = .25, stat = "identity", 
           position = "fill") +
  #facet_wrap(~Drug) +
  #scale_fill_manual("Dose (uM)", values = clr_values) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  #scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Stimulator") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("stackedbarplot_clustersPerStimulator_vehicle.pdf",
       width = 3, height = 2.25, unit = "in", dpi = 1200)

##########################################
# Remake stacked barplots with annotations
cell_counts_summary = group_by(cdat, Drug, Relative_Dose, Stimulator, cellType_broad) %>% 
  summarise(n_cells = n()) 

# stacked bar of proportion cluster per stimulator
cell_counts_summary %>%
  dplyr::select(Drug, Relative_Dose, Stimulator,cellType_broad, n_cells) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = Relative_Dose, y =n_cells, fill = cellType_broad), 
           color = "black", size = .25, stat = "identity", 
           position = "fill") +
  facet_wrap(~Drug+Stimulator) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  xlab("Stimulator") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("stackedbarplot_cellTypesPerDrugDose.pdf",
       width = 6, height = 4.5, unit = "in", dpi = 1200)


# Vehical treated only 
cell_counts_summary %>% 
  filter(Relative_Dose == 0) %>% 
  group_by(cellType_broad, Stimulator) %>% summarise(n_cells_new = sum(n_cells)) %>% 
  #dplyr::select(Drug, Relative_Dose, Stimulator, clusters_mon3, n_cells) %>%
  dplyr::select(Stimulator, cellType_broad, n_cells_new) %>%
  distinct() %>%
  ggplot() +
  geom_bar(aes(x = Stimulator, y =n_cells_new, fill = cellType_broad), 
           color = "black", size = .25, stat = "identity", 
           position = "fill") +
  #facet_wrap(~Drug) +
  #scale_fill_manual("Dose (uM)", values = clr_values) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +    
  #scale_x_continuous(breaks = seq(1,10,1))+
  xlab("Stimulator") +
  ylab("Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("stackedbarplot_cellTypesPerStimulator_vehicle.pdf",
       width = 3, height = 2.25, unit = "in", dpi = 1200)

###################################### 
# Using Replicate information, plot cell type
# proportions for each condition. 
# want to look for significant changes 
# (i.e. activated Tcells after bead stimulation)

# Vehicle treated cells only
cell_counts_condition = data.frame(prj@cellColData) %>% 
  filter(Relative_Dose == 0) %>% 
  group_by(Replicate, Stimulator, Drug) %>% 
  summarise(total_cells = n()) %>% 
  mutate(condition = paste(Drug, Stimulator, Replicate, sep = "_"))

cell_counts_summary = data.frame(prj@cellColData) %>% 
  filter(Relative_Dose == 0) %>% 
  group_by(Replicate, Drug, Stimulator, cellType_broad) %>% 
  summarise(n_cells = n()) %>% 
  mutate(condition = paste(Drug, Stimulator, Replicate, sep = "_")) %>% 
  left_join(cell_counts_condition, by = "condition") %>% 
  select(
    Replicate = Replicate.x,
    Drug = Drug.x, 
    Stimulator = Stimulator.x, 
    cellType_broad, 
    n_cells, 
    total_cells) %>% 
  mutate(percent_cellType = (n_cells/total_cells)*100)

# make box plots of cell proportion faceted by Drug 
# prepare separate plots for each cell type 

celltypes = cell_counts_summary %>% 
  distinct(cellType_broad)

for (ct in celltypes$cellType_broad){
  pdf(file = paste0("BoxPlot_PercentRecoveredCellType_VehicleOnly_byStimulator_", ct,".pdf"), 
      width = 6, height = 4.5)
  print(  
    cell_counts_summary %>% 
      filter(cellType_broad == ct) %>% 
      ggplot(aes(x =  Stimulator, y = percent_cellType)) +
      geom_boxplot() +
      facet_wrap("~Drug") +
      scale_color_brewer(palette='Set1') +
      xlab("Stimulator") +
      ylab("% recovered cells") +
      ggtitle(paste0("Percent recovered ", ct, " (vehicle only)")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      stat_compare_means(label = "p.signif", method = "t.test",  
                         ref.group = "StimD", hide.ns = TRUE, 
                         label.y.npc = 0.85, col = "red")  # Pairwise comparison against reference

  )
  dev.off()
}

## Use each vehicle treatement from separate drugs as replicates:
for (ct in celltypes$cellType_broad){
  pdf(file = paste0("BoxPlot_PercentRecoveredCellType_AllVehicleReps_byStimulator_", ct,".pdf"), 
      width = 3, height = 3)
  print(  
    cell_counts_summary %>% 
      filter(cellType_broad == ct) %>% 
      ggplot(aes(x =  Stimulator, y = percent_cellType)) +
      geom_boxplot() +
      #facet_wrap("~Drug") +
      scale_color_brewer(palette='Set1') +
      xlab("Stimulator") +
      ylab("% cells") +
      ggtitle(paste0("% Recovered ", ct)) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      stat_compare_means(label = "p.signif", method = "t.test",  
                         ref.group = "StimD", hide.ns = TRUE, 
                         label.y.npc = 0.85, col = "red")
    #theme_bw()
  )
  dev.off()
}

######################################
# Prepare heatmaps of % cell types recovered for each drug across doses
######################################

cell_counts_condition_dose = data.frame(prj@cellColData) %>% 
  #group_by(Stimulator, Drug, Relative_Dose) %>% 
  group_by(Drug, Relative_Dose) %>% 
  summarise(total_cells = n()) %>% 
 #mutate(condition = paste(Drug, Stimulator, Relative_Dose, sep = "_"))
  mutate(condition = paste(Drug,Relative_Dose, sep = "_"))

cell_counts_summary_dose = data.frame(prj@cellColData) %>% 
  #group_by(Drug, Stimulator, Relative_Dose, cellType_broad) %>% 
  group_by(Drug, Relative_Dose, cellType_broad) %>% 
  summarise(n_cells = n()) %>% 
  #mutate(condition = paste(Drug, Stimulator, Relative_Dose, sep = "_")) %>% 
  mutate(condition = paste(Drug, Relative_Dose, sep = "_")) %>% 
  left_join(cell_counts_condition_dose, by = "condition") %>% 
  select(
    Relative_dose = Relative_Dose.x,
    Drug = Drug.x, 
    #Stimulator = Stimulator.x, 
    cellType_broad, 
    n_cells, 
    total_cells) %>% 
  mutate(percent_cellType = (n_cells/total_cells)*100)#,
        # Stimulator = recode(Stimulator, StimA = "Allo", StimD = "Auto"))

celltypes = cell_counts_summary_dose %>% 
  distinct(cellType_broad)

for (ct in cell_counts_summary_dose$cellType_broad){
  cell_counts_summary_dose_f = filter(cell_counts_summary_dose, cellType_broad == ct)
  # plot heatmap
  pdf(file = paste0("HeatMap_percent_", ct, "_Recovered_acrossDoses_noStimFacet.pdf"), 
      width = 4, height = 3.5)
  print(  
    ggplot(cell_counts_summary_dose_f, aes(x = Relative_dose, y = Drug, fill = percent_cellType)) + 
      #facet_wrap(~Stimulator) +
      geom_tile(color = "gray") +
      scale_fill_gradient(low = "steelblue", high = "Firebrick") + 
      ylab("Drug") +
      xlab("Relative Dose") + 
      labs(fill = paste0("Percent ", ct))
  )
  dev.off()
}



# Prepare heatmaps of fold change in % recovered 
# cell types (StimA/StimD) for each drug across doses

stimA_cellCounts = filter(cell_counts_summary_dose, Stimulator == "StimA") %>% 
  mutate(condition = paste(cellType_broad, Drug, Relative_dose, sep = "_"))
stimD_cellCounts = filter(cell_counts_summary_dose, Stimulator == "StimD") %>% 
  mutate(condition = paste(cellType_broad, Drug, Relative_dose, sep = "_"))

df_compare = left_join(stimA_cellCounts, stimD_cellCounts, by = "condition") %>% 
  select(
    Relative_dose = Relative_dose.x,
    Drug = Drug.x, 
    cellType_broad = cellType_broad.x, 
    n_cells_A = n_cells.x, 
    total_cells_A = total_cells.x,
    percent_cellType_A = percent_cellType.x,
    n_cells_D = n_cells.y, 
    total_cells_D = total_cells.y,
    percent_cellType_D = percent_cellType.y) %>% 
  mutate(foldChangeA_D = percent_cellType_A/percent_cellType_D)

###
celltypes = df_compare %>% 
  distinct(cellType_broad)

for (ct in celltypes$cellType_broad){
  df_compare_f = filter(df_compare, cellType_broad == ct)
  # plot heatmap
  pdf(file = paste0("HeatMap_FCpercent_", ct, "_Recovered_stimA_stimD.pdf"), 
      width = 6, height = 3.5)
  print(  
  ggplot(df_compare_f, aes(x = Relative_dose, y = Drug, fill = foldChangeA_D)) + 
    geom_tile(color = "gray") +
    scale_fill_gradient(low = "steelblue", high = "Firebrick") + 
    ylab("Drug") +
    xlab("Relative Dose") + 
    labs(fill = paste0("FC Percent", ct,  "\n (stimA/stimD)"))
  )
  dev.off()
}


