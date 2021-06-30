basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB8"))
setwd(paste0(out_dir, "results/NB8"))

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
  library(Gviz)
  library(biomaRt)
})

set.seed(2017) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

# load processed cds with processed naive to active T-cell trajectory 
cds = readRDS(file = paste0(out_dir, "results/NB7/cds_NaiveActive_Tcells_PT"))
cds <- detect_genes(cds)

# subset to peaks open in > 0.1% cells
cutoff = ncol(cds)*.001
cds = cds[rowData(cds)$num_cells_expressed > cutoff, ]

###########################################
###########################################
# this section takes a long time (~1-2 hrs) 
# Only need to run fit_models once: load the saved output in next section 

# fit regression models for each peak across pseudotime of T-cell activation 
peak_fits <- fit_models(
  cds = cds, 
  #model_formula_str = "~Pseudotime + nFrags", 
  model_formula_str = "~ splines::ns(Pseudotime, df = 3) + nFrags", 
  expression_family = "binomial")

# get coeficients from fit_models
peak_fit_coefs = coefficient_table(peak_fits)

peak_fit_coefs_f = dplyr::select(peak_fit_coefs, 
                                 peak, num_cells_expressed, status, term,
                                 estimate, std_err, test_val, p_value, 
                                 normalized_effect, model_component, q_value)

write.table(peak_fit_coefs_f, 
            file =  "Tcell_Activation_PT_DApeak_spline_coefs.txt", 
            row.names=F, 
            quote=F, 
            sep="\t")

###########################################
###########################################
peak_coefs = read.table("Tcell_Activation_PT_DApeak_spline_coefs.txt", sep = "\t",
                        head = TRUE, stringsAsFactors = FALSE)


peak_coefs_sig = filter(peak_coefs, grepl("Pseudotime", x = term), q_value < 0.05) %>%
  mutate(direction = ifelse(estimate > 0, "opening", "closing")) %>% 
  distinct(peak, .keep_all = TRUE) %>% 
  arrange(q_value)

######################
# get normalized pseudobulk bigwig tracks grouped by pseudotime bins
prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 

cdat = colData(cds) %>% 
  data.frame() %>% 
  tibble::rownames_to_column(var = "cellNames") %>% 
  mutate(PT_bin = paste0("PTbin_", Pseudotime_bin))

# count number of cells in each bin
group_by(cdat, PT_bin) %>% 
  summarise

# isolate cells used in the T-cell trajectory analysis
prj_f = prj[prj$cellNames %in% cdat$cellNames,]

# add PT bin data to temp project 
prj_f <- addCellColData(ArchRProj = prj_f, data = cdat$PT_bin,
                      cells = cdat$cellNames, name = "Pseudotime_bin", 
                      force = TRUE)

getGroupBW(
  ArchRProj = prj_f,
  groupBy = "Pseudotime_bin",
  normMethod = "ReadsInTSS", # creates a group scale factor = 10k/sum(reads in TSS) 
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)










