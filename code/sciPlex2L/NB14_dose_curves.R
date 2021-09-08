basepath = "/net/trapnell/vol1/home/gtb7/projects/scichem_ATAC/190521_scichem2_AllPlates/"
bin_directory = paste0(basepath, "analysis/bin/")
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB14"))
setwd(paste0(out_dir, "results/NB14"))

suppressPackageStartupMessages({
library(ArchR)
library(devtools)
library(monocle3)
library(tidymodels)
library(purrr)
library(dplyr)
library(drc)
library(glmnet)
})


################

count_cells = function(cds, normalization_terms=c("cell_type", "culture_plate"), 
                       normalization_formula_str=NULL, pseudodose=0){
  
  # Get a cell count for cells recovered from each well
  sci_chem_counts = 
    colData(cds) %>% 
    as.data.frame() %>%
    group_by(cell_type, sample, replicate, culture_plate, well_oligo, treatment, dose, vehicle, Relative_dose) %>% 
    summarize(cells = n()) %>% 
    ungroup()
  
  # Calculate the number of each cell type recovered
  sci_chem_counts = 
    sci_chem_counts %>% 
    group_by(cell_type) %>% 
    mutate(total_cells_for_sample = sum(cells)) %>% 
    ungroup()
  
  # Calculate a size factor that is the ratio of the total number of cells recovered over the mean per cell
  sci_chem_counts = 
    sci_chem_counts %>% 
    group_by(cell_type) %>% 
    mutate(cell_count_size_factor =  total_cells_for_sample / mean(total_cells_for_sample) ) %>% 
    ungroup()
  
  # Bring rownames (cell names)) into a column 
  sci_chem_counts = tibble::rownames_to_column(sci_chem_counts)
  
  # Build the model string from the normalization terms provided
  if (is.null(normalization_formula_str)){
    normalization_terms = lapply(normalization_terms, function(x) { 
      
      # Make sure that the normalization terms provided have more than one value
      if(nrow(unique(sci_chem_counts[,x])) > 1) {
        return(x)
      } else{ 
        return (NA)
      }
    })
    
    # Select those normalization terms that have non-NA values
    normalization_terms = normalization_terms[is.na(normalization_terms) == FALSE]
    
    # Build the model string for normalization
    normalization_formula_str = paste("~", paste(normalization_terms, collapse="*"))
    if (length(normalization_terms) > 0){
      normalization_formula = as.formula(normalization_formula_str)
    }
    else{
      normalization_formula = NULL
    } 
  }
  
  if (is.null(normalization_formula) == FALSE){
    
    # Use the cell counts and the model string to create a model matrix
    plate_model_matrix = 
      model.matrix(normalization_formula, data=sci_chem_counts)
    
    # Predict the contributions of the normalization terms to the cell count
    plate_effects_model = 
      glmnet(x = plate_model_matrix, 
             y = as.vector(sci_chem_counts[row.names(plate_model_matrix),]$cells), 
             family="poisson", 
             lambda=0)
    
    # Calculated expected recovery based on fitted values
    fitted_vals = predict(plate_effects_model, newx=plate_model_matrix, type="response")
    
    adj_r_squared = plate_effects_model$dev.ratio
    
    model_coefs = coef(plate_effects_model)
    model_coefs_mat = as.matrix(model_coefs)
    model_coefs_mat = as.data.frame(model_coefs_mat)
    
    cell_count_adjustment = fitted_vals - exp(model_coefs[1,1])
    
    sci_chem_counts$fitted_cells = NA
    sci_chem_counts$fitted_cells[sci_chem_counts$rowname %in% row.names(plate_model_matrix)] = fitted_vals
    
    sci_chem_counts$norm_cells = NA
    sci_chem_counts$norm_cells[sci_chem_counts$rowname %in% row.names(plate_model_matrix)] = 
      sci_chem_counts$cells[sci_chem_counts$rowname %in% row.names(plate_model_matrix)] - cell_count_adjustment
    
    sci_chem_counts$norm_cells[is.na(sci_chem_counts$norm_cells)] = 
      sci_chem_counts$cells[is.na(sci_chem_counts$norm_cells)]
    sci_chem_counts$norm_cells[sci_chem_counts$norm_cells < 0] = 0
  }else{
    message("After pruning un-needed normalization factors, no normalization needed.")
    sci_chem_counts$norm_cells = sci_chem_counts$cells
  }
  
  # Isolate vehicle treated cells
  vehicle_counts = 
    sci_chem_counts %>% 
    filter(vehicle) 
  
  
  sci_chem_counts = 
    sci_chem_counts %>% 
    dplyr::select(cell_type, sample, replicate, culture_plate, well_oligo, treatment, dose, norm_cells, cells, vehicle, Relative_dose)
  
  vehicle_counts = 
    vehicle_counts %>% 
    dplyr::select(cell_type, sample, replicate, culture_plate, well_oligo, norm_cells, cells, vehicle, Relative_dose)
  
  # For each treament add a vehicle treatment row where the dose is set to the pseudodose
  vehicle_count_rows = 
    unique(sci_chem_counts$treatment) %>% 
    purrr::map_dfr(.f = function(tmnt) {
      vehicle_count_rows = vehicle_counts
      vehicle_count_rows$treatment = tmnt
      vehicle_count_rows$dose = pseudodose
      vehicle_count_rows
    })
  
  
  
  # vehicle_count_rows = full_join(vehicle_counts, 
  #                                dplyr::select(sci_chem_counts, cell_type, sample, replicate, culture_plate, well_oligo,  treatment, vehicle) %>% distinct(),
  #                                by=c("cell_type", "sample", "replicate")) %>% distinct()
  # vehicle_count_rows$dose = pseudodose
  
  vehicle_count_rows  = 
    vehicle_count_rows %>% 
    dplyr::select(cell_type, sample, replicate, culture_plate, well_oligo, treatment, dose, vehicle, norm_cells, cells, Relative_dose)
  
  # FIXME: put this back!
  
  sci_chem_counts = rbind(sci_chem_counts, vehicle_count_rows)
  
  return(sci_chem_counts)
}

fit_drm_helper = function(df, verbose=FALSE){
  tryCatch({
    #drm(norm_cells ~ dose, data = df, fct=LL.4(), type="Poisson" )
    drm(norm_cells ~ dose, data = df, fct=LL.4())
  }, 
  error = function(e) { 
    if (verbose)
      print (e); 
    NA 
  })
}

fit_dose_response_models = function(cell_count_table, verbose=FALSE){
  sci_chem_models = cell_count_table %>% 
    nest(-cell_type, -treatment) %>% 
    mutate(
      drc_model = map(.f = fit_drm_helper, .x = data, verbose=verbose)
    ) 
  sci_chem_models
}

extract_terms_helper = function(model){
  if (class(model) == "drc"){
    coefs = coef(model)
    ed = tryCatch( { ED(model, 50, interval="delta", display=FALSE) } , error = function(e) { ed = rep(NA, 4) })
    # Isolate the individual parameters and store in a small tibble
    tmp = tibble(
      steepness = coefs[1],
      lower_response_limit = coefs[2],
      upper_response_limit = coefs[3],
      #ec_param = coefs[3],
      ec50 = ed[1],
      ec50_lower = ed[3],
      ec50_upper = ed[4]
    )
    tmp
  }else{
    return (tibble(
      steepness = NA,
      lower_response_limit = NA,
      upper_response_limit = NA,
      ec50 = NA,
      ec50_lower = NA,
      ec50_upper = NA
    ))
  }
  
}

extract_terms = function(sci_chem_models){
  sci_chem_models = sci_chem_models %>%
    mutate(drc.terms = purrr::map(.f = extract_terms_helper, .x = drc_model)) %>% 
    unnest(drc.terms)
  return(sci_chem_models)
}


test_reduction_helper = function(model){
  if (class(model) == "drc"){
    tryCatch( {  
      cfs = coef(model)
      #lower_response_limit = max(0, cfs[2])
      #upper_response_limit = max(0, cfs[3])
      lower_response_limit = predict(model, newdata=data.frame(dose=c(min(model$origData$dose))))
      upper_response_limit = predict(model, newdata=data.frame(dose=c(max(model$origData$dose))))
      fold_change = upper_response_limit / lower_response_limit
      #print (paste(upper_response_limit, lower_response_limit, fold_change))
      #print (noEffect(model))
      tbl = noEffect(model)
      tbl = c(fold_change, tbl)
      
      #print (class(tbl))
      names(tbl) = c("fold_change", 
                     "Chi-square test",
                     "Df", 
                     "p-value") 
      tbl = as_tibble(t(tbl))
      
      #print (tbl)
      return(tbl)
    }, error = function(e) { 
      print (e)
      tibble("fold_change" = NA,
             "Chi-square test" = NA,
             "Df" = NA, 
             "p-value" = NA) 
    })
  }else{
    tibble("fold_change" = NA,
           "Chi-square test" = NA,
           "Df" = NA, 
           "p-value" = NA) 
  }     
}

test_dose_curve = function(sci_chem_models){
  sci_chem_sig_response = sci_chem_models %>%
    mutate(fit_test_res = purrr::map(.f = test_reduction_helper, .x = drc_model)) %>% 
    unnest(fit_test_res) %>%
    dplyr::select(cell_type, treatment, "fold_change", "Chi-square test", "Df", "p-value")
  return(sci_chem_sig_response)
}

extract_dose_model_stats = function(sci_chem_counts, sci_chem_models){
  dose_range_df = sci_chem_counts %>% group_by(cell_type, treatment) %>% summarize(min_dose = min(dose), max_dose = max(dose))
  sci_chem_models = extract_terms(sci_chem_models)
  coefficient_df = sci_chem_models %>% dplyr::select(cell_type, treatment, ec50, ec50_lower, ec50_upper, upper_response_limit, lower_response_limit)
  coefficient_df$lower_response_limit = pmax(0, coefficient_df$lower_response_limit)
  coefficient_df$upper_response_limit = pmax(0, coefficient_df$upper_response_limit)
  coefficient_df = coefficient_df %>% mutate(fold_change = (upper_response_limit - lower_response_limit)/ upper_response_limit) 
  coefficient_df = left_join(coefficient_df, dose_range_df)
  coefficient_df$ec50[coefficient_df$ec50 > coefficient_df$max_dose] = NA
  coefficient_df$ec50_lower[coefficient_df$ec50_lower > coefficient_df$max_dose] = max(coefficient_df$max_dose) # FIXME: this isn't generic. Could be different max doses for different drugs
  coefficient_df$ec50_lower[coefficient_df$ec50_lower < coefficient_df$min_dose] = 0
  coefficient_df$ec50[coefficient_df$ec50 < coefficient_df$min_dose] = NA
  coefficient_df$ec50_upper[coefficient_df$ec50_upper > coefficient_df$max_dose] = max(coefficient_df$max_dose) # FIXME: this isn't generic. Could be different max doses for different drugs
  coefficient_df$ec50_upper[coefficient_df$ec50_upper < coefficient_df$min_dose] = 0
  return(coefficient_df)
}

gen_dose_range_helper = function(model){
  tryCatch({
    if (class(model) == "drc"){
      eps_dose = min(model$origData$dose[model$origData$dose>0], na.rm=T)
      exp(seq(log(min(model$origData$dose, na.rm=T) + eps_dose/2), log(max(model$origData$dose, na.rm=T)), length=100))
    }else{
      return (rep(NA, 100))
    }
  }, error = function(e) { print(e); return (rep(NA, 100)) }
  )
}

extract_curve_helper = function(model){
  tryCatch({
    if (class(model) == "drc"){
      eps_dose = min(model$origData$dose[model$origData$dose>0], na.rm=T)
      newdata =  expand.grid(dose=exp(seq(log(min(model$origData$dose, na.rm=T) + eps_dose/2), log(max(model$origData$dose, na.rm=T)), length=100)))
      #newdata = data.frame(dose=seq(min(model$origData$dose), max(model$origData$dose),length.out=100))
      res = predict(model, newdata=newdata)
    }else{
      return (rep(NA, 100)) 
    }
  }, 
  warning = function(e) { return (rep(NA, 100)) },
  error = function(e) { return (rep(NA, 100)) }
  )
}

extract_curves = function(sci_chem_models){
  sci_chem_curves = sci_chem_models %>%
    mutate(dose = purrr::map(.f = gen_dose_range_helper, .x = drc_model)) %>%
    mutate(fit = purrr::map(.f = extract_curve_helper, .x = drc_model)) %>%
    dplyr::select(cell_type, treatment,  dose, fit) %>% 
    unnest(dose, fit)
  return(sci_chem_curves)
}

plot_dose_response_curves = function(sci_chem_counts, sci_chem_models, color_by="replicate", ncol=NULL, point_size = 1){
  # sci_chem_data = sci_chem_models %>%
  #     #mutate(fit = purrr::map(.f = extract_curve, .x = drc_model)) %>%
  #     dplyr::select(cell_type, treatment, data) %>% 
  #     unnest(data)
  sci_chem_data = sci_chem_counts
  coefficient_df = extract_dose_model_stats(sci_chem_counts, sci_chem_models)
  
  sci_chem_curves = extract_curves(sci_chem_models)
  min_plot_dose = 0.05
  ggplot(aes(dose + min_plot_dose, norm_cells), data=sci_chem_data) + 
    geom_point(aes_string(color=color_by), size = point_size, stroke = 0, position="jitter") +
    geom_line(aes(dose, fit), data=sci_chem_curves) + 
    geom_rect(aes(x=ec50, y=0, xmin=ec50_lower + min_plot_dose, xmax=ec50_upper, ymin=0, ymax=Inf), 
              alpha=I(0.1), data=coefficient_df) + 
    geom_vline(aes(xintercept=ec50), linetype="dotted", data=coefficient_df) + 
    scale_x_log10() + 
    scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() + 
    facet_wrap(treatment~cell_type, scales="free_y", ncol=ncol)
}

#################################
## Greg usage below
## Create cds from ArchR project

prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered"), 
                       showLogo = FALSE) 

pMat = getMatrixFromArrow(ArrowFile = paste0(out_dir,"sc2_filtered/ArrowFiles/sc2.arrow"),
                          ArchRProj = prj, 
                          cellNames = prj$cellNames,
                          useMatrix = "PeakMatrix", 
                          binarize = TRUE)

# Create CDS from peak Matrix (SummarizedExperiment)
rd = data.frame(rowRanges(pMat)) %>% 
  rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end))

# peaks in this list are in a DIFFERENT ORDER 
# relative to pMat (rowRanges above is correct)
pd = data.frame(getPeakSet(prj)) %>% 
  rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end), 
         gene_short_name = ifelse(peakType == "Promoter" & distToGeneStart < 500, nearestGene, "NA")) %>% 
  dplyr::select(rowname, score, replicateScoreQuantile, groupScoreQuantile, 
                Reproducibility, GroupReplicate, distToGeneStart,
                nearestGene, peakType, distToTSS, nearestTSS,
                gene_short_name, GC, idx, peak)

a_rd = left_join(rd, pd, by = "peak") %>% 
  dplyr::select(rowname = rowname.x, idx = idx.x, score, replicateScoreQuantile, groupScoreQuantile, 
         Reproducibility, GroupReplicate, distToGeneStart,
         nearestGene, peakType, distToTSS, nearestTSS,
         gene_short_name, GC, peak)

row.names(a_rd) <- a_rd$peak
row.names(pMat) <- a_rd$peak
cds_atac = monocle3::new_cell_data_set(assays(pMat)$PeakMatrix, 
                                       cell_metadata = colData(pMat), 
                                       gene_metadata = a_rd)
#format colData of atac_cds

colData(cds_atac)$cell =
  row.names(colData(cds_atac))

colData(cds_atac)$log_nFrags =
  log(colData(cds_atac)$nFrags)

colData(cds_atac)$dose_character =
  factor(colData(cds_atac)$Relative_dose,
         levels = c("0", "0.1", "0.5", "1", "5", "10", "50", "100"))

colData(cds_atac)$new_treatment_label =
  sapply(colData(cds_atac)$treatment, function(x) {
    if (grepl("BMS", x))
      return("BMS345541")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
  })

colData(cds_atac)$solvent =
  sapply(pData(cds_atac)$treatment, function(x) {
    if (grepl("Dex", x))
      return("Ethanol")
    else
      return("DMSO")
  })



cds = cds_atac
#apply dose response functions to loaded cds
CellCounts = count_cells(cds = cds)
#CellCounts_filt$replicate = as.factor(CellCounts_filt$replicate)
modelFits = fit_dose_response_models(CellCounts)
ModelTerms = extract_terms(modelFits)
doseCurve = test_dose_curve(modelFits)
doseModelStats = extract_dose_model_stats(sci_chem_counts = CellCounts, sci_chem_models = modelFits)
doseCurves_F = extract_curves(modelFits)

#plot_dose_response_curves(sci_chem_counts = CellCounts, sci_chem_models = modelFits)


#manually plot 
CellCounts$dose_character <- as.character(CellCounts$Relative_dose)
CellCounts$dose_character <- factor(CellCounts$dose_character,levels = c("0", "0.1", "0.5", "1", "5","10", "50", "100")) 

min_plot_dose = 0
pdf("sciCHEM_ATAC_DoseResponse.pdf", width = 3.5, height = 3.5)
ggplot(aes(dose + min_plot_dose, norm_cells), data=CellCounts) + 
  geom_point(aes_string(color=CellCounts$dose_character), size = 1.5, stroke = 0, position="jitter") +
  geom_line(aes(dose, fit), data=doseCurves_F) +
  #geom_vline(aes(xintercept=ec50), linetype="dotted", data=ModelTerms) +
  xlab("Dose [uM]") +
  ylab("Cell counts") +
  #scale_x_log10(breaks = c(1, 10, 100, 1000, 10000)) + 
  scale_x_log10() + 
  facet_wrap(treatment~cell_type, scales="free", ncol=NULL) +
  scale_color_manual("Dose",values = c("0"="gray", "0.1"="#1D1147FF", "0.5"="#51127CFF", "1"="#822681FF", "5"="#B63679FF",
                                       "10"="#E65164FF", "50" = "#FB8861FF", "100"="#FEC287FF"))+
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  theme(legend.position = "none") +
  monocle3:::monocle_theme_opts() 
dev.off()


###################
Dex_counts = filter(CellCounts, treatment == "Dex")
dexModel = drm(norm_cells ~ dose, data = Dex_counts, fct=LL.4())



