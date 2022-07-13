basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/211103_3Level_scichem_MLR_Drugs/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB17"))
setwd(paste0(out_dir, "results/NB17"))

suppressPackageStartupMessages({
  library(ArchR)
  #library(Seurat)
  library(monocle3)
  library(cicero)
  library(devtools)
  #load_all('/net/trapnell/vol1/home/gtb7/git/github/monocle3_dev') # latest developer branch of monocle3
  library(dplyr)
  library(tidymodels)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
  library(RColorBrewer)
  #library(ggpubr)
  library(viridis)
  library(snowfall)
  library(furrr)
  library(gt)
  library(forcats)
  library(lmtest)
  plan(multicore)
})

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "MLR_drug_filtered_annotated"), 
                       showLogo = FALSE) 

# Load CDS with T-cell trajectory pseudotime information 
cds = readRDS(file = paste0(out_dir, "results/NB8/cds_NaiveActive_Tcells_PT"))

###############
# Look at how FRIP scores change across T-cell activation, faceted by drug and dose. 
cds_ff = cds[,colData(cds)$cellType_broad == "Activated_Tcell" | colData(cds)$cellType_broad == "Naive_Tcell"  & colData(cds)$Stimulator == "StimA"]

cdat = data.frame(colData(cds_ff))
drugs = dplyr::distinct(cdat, Drug)

for(d in drugs$Drug){
  cdat %>%
    filter(Drug == d) %>% 
    ggplot() +
    geom_boxplot(aes(x = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)),
                     y = FRIP, 
                     #fill = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)))) +
                     fill = factor(Relative_Dose, levels = c(0, 1, 10, 100)))) +
    #facet_wrap(~Drug) +
    theme_bw() +
    xlab("Pseudotime bin") +
    ylab("FRIP") +
    #ylim(4,15) +
    #scale_fill_viridis_d(option ='inferno', name = "Rel. Dose") +
    scale_fill_brewer(palette = "RdBu") +
    #theme(legend.position = "none") 
    ggsave(paste0("BoxPlot_FRIP_overPT_byDrugDose3_", d, ".pdf"),
           #ggsave(paste0("BoxPlot_FRIT_overPT_byDrugDose2_FacetDrugs.pdf"),
           height = 4 , width = 5, unit = "in")
}

###
# Model broad cell characters across T-cell trajectory
# ex) linear model of FRIP ~ ns::spline(Pseduotime, df = 3)*Drug 
# Get a predicted curve 
###

# function to generate linear model
full_model_PT_drugDose = function(drug = "SAHA"){
  cdat_f = cdat %>%
    filter(Drug == drug) 
  model_pt = glm(formula = "FRIP ~ splines::ns(Pseudotime, df = 3)*Relative_Dose", data = cdat_f, family = "gaussian") 
  return(model_pt)
}

reduced_model_PT = function(drug = "SAHA"){
  cdat_f = cdat %>%
    filter(Drug == drug) 
  reduced_model_pt = glm(formula = "FRIP ~ splines::ns(Pseudotime, df = 3)", data = cdat_f, family = "gaussian") 
  return(reduced_model_pt)
}

#function to fit model to PT_bin group averages: 
predict_from_PTbins = function(model_pt, drug = "SAHA"){
  cdat_f = cdat %>%
    filter(Drug == drug) 
  PT_bin_means = group_by(cdat_f, Pseudotime_bin, Relative_Dose) %>% 
    summarise(FRIP = mean(FRIP), 
              Pseudotime = mean(Pseudotime))
  
  PT_bin_means$mod_preds = predict(model_pt, PT_bin_means)
  return(PT_bin_means)
}

# function to plot model fits across pseudotime, separately for each dose
plot_preds  = function(model_preds, drug = "SAHA"){
  ggplot(model_preds, aes(x = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)),
                          y = mod_preds, 
                          group = factor(Relative_Dose, levels = c(0, 1, 10, 100)),
                          colour= factor(Relative_Dose, levels = c(0, 1, 10, 100)))) +
    geom_line() +
    theme_bw() +
    ggtitle(drug) +
    xlab("Pseudotime bin") +
    ylab("FRIP (pred)") +
    scale_colour_viridis_d(option ='inferno', name = "Rel. Dose") +
    ggsave(paste0("GLMfits_FRIP_overPT_byDrugDose2_", drug, ".pdf"),
           height = 2 , width = 3.5, unit = "in")
}

for(d in drugs$Drug){
  cat("modeling FRIP by PT*DrugDose for ", d, " treated cells" )
  mod_full = full_model_PT_drugDose(drug = d)
  mod_red = reduced_model_PT(drug = d)
  lrt = lrtest(mod_full, mod_red) # Likelihood ratio tests 
  print(lrt)
  preds = predict_from_PTbins(model_pt = mod_full, drug = d)
  plot_preds(model_preds = preds, drug = d)
}

####################
# generate dataframe for faceted plot
df_preds = data.frame(Pseudotime_bin = integer(), Relative_Dose = character(),  
                      FRIP = double(), Pseudotime = double(), mod_preds = double(), 
                      Treatment = character())
for(d in drugs$Drug){
  cat("modeling FRIP by PT*DrugDose for ", d, " treated cells" )
  mod_full = full_model_PT_drugDose(drug = d)
  preds = predict_from_PTbins(model_pt = mod_full, drug = d) 
  preds$treatment = d
  df_preds = rbind(df_preds, preds)
}

# make faceted plot 
ggplot(df_preds, aes(x = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)),
                        y = mod_preds, 
                        group = factor(Relative_Dose, levels = c(0, 1, 10, 100)),
                        colour= factor(Relative_Dose, levels = c(0, 1, 10, 100)))) +
  geom_line() +
  facet_wrap(~treatment) +
  theme_bw() +
  xlab("Pseudotime bin") +
  ylab("FRIP (pred)") +
  scale_color_brewer(palette = "RdBu", name = "Rel. Dose", direction = -1)
ggsave(paste0("GLMfits_FRIP_overPT_byDrugDose2_facet.pdf"),
       height = 4 , width = 5, unit = "in")

# make faceted boxplot: 
ggplot(cdat, aes(x = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)),
                 y = FRIP, 
                 #fill = factor(Pseudotime_bin, levels = c(1,2,3,4,5,6,7,8,9,10)))) +
                 fill = factor(Relative_Dose, levels = c(0, 1, 10, 100)))) +
  geom_boxplot()+
  facet_wrap(~Drug) +
  theme_bw() +
  xlab("Pseudotime bin") +
  ylab("FRIP") +
  #ylim(4,15) +
  scale_fill_brewer(palette = "RdBu", name = "Rel. Dose", direction = -1) +
ggsave(paste0("BoxPlot_FRIP_overPT_byDrugDose2_FacetDrugs.pdf"),
         height = 6 , width = 7.5, unit = "in")
