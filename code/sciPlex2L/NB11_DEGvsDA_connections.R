basepath = "github/"
bin_directory = paste0(basepath, "analysis/bin/")
out_dir = paste0(basepath, "analysis/archr/sc2_integratedRNA_v2/")
dir.create(paste0(out_dir,"Plots/NB9/"))

setwd(paste0(out_dir,"Plots/NB9/"))

suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(cicero)
  library(ggplot2)
  library(ggridges)
  library(tidymodels)
  library(devtools)
  library(furrr)
  library(viridis)
  library(piano)
  library(UpSetR)
  library(snowfall)
  library(ggrastr)
  library(ggpubr)
  plan(multicore)
})

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e7)
set.seed(1989) 

prj = loadArchRProject(path = out_dir,
                       showLogo = FALSE)

# load DE Gene list: 
gene_model_res = read.table(file = paste0(out_dir, "Plots/NB8/sciPlex_RNA_A549_small_screen_DEGs.txt"), head = T)
dex_degs = gene_model_res %>%
  filter(grepl("dose_Dex", term)) %>%
  filter(q_value < 0.05) %>%
  dplyr::select(id, gene_short_name, term, estimate, q_value)

# load DA gene Score model list:
genescore_model_res = read.table(file = paste0(out_dir, "Plots/NB8/sciPlex_RNA_A549_small_screen_DAgeneScores.txt"), head = T)

# load DA site list: 
peak_model_res = read.table(file = paste0(out_dir, "Plots/NB8/sciPlex_RNA_A549_small_screen_DApeaks.txt"), head = T)
dex_das = peak_model_res %>%
  filter(grepl("dose_Dex", term)) %>%
  filter(q_value < 0.05) %>%
  dplyr::select(peak, term, estimate, q_value)

# retreive peaklist from archR project. 
peakSet = getPeakSet(prj)

# retreive peak to gene links from archR project. 
p2g <- getPeak2GeneLinks(
  ArchRProj = prj,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)

###############################################
p_info = data.frame(peakSet)
p_info = mutate(p_info, peak = paste0(seqnames, "_", start, "_", end))
Conn_df = data.frame(p2g)
P_df = data.frame(p2g@metadata$peakSet)
G_df = data.frame(p2g@metadata$geneSet)

# function returns all connections (and peak info) for a gene
getGeneConns = function(gene, peak_info = p_info, conn_df = Conn_df, peak_df = P_df, 
                        gene_df = G_df){
  valid_genes = gene_df$name
  if(!(gene %in% valid_genes)){message("gene not in ArchR project")}
  
  else{   
    g_idx = filter(gene_df, name == gene) %>% 
      select(idx) 
    p_idx = filter(conn_df, idxRNA == g_idx$idx) %>% 
      select(idxATAC, Correlation) 

    if(length(p_idx$idxATAC) > 0){
      gene_conns = peak_df[p_idx$idxATAC,]
      gene_conns$p2gCorr = p_idx$Correlation

      gene_conns = mutate(gene_conns, peak = paste0(seqnames, "_", start, "_", end)) %>% 
        select(peak, p2gCorr)
      gene_conns$connectedGene = gene
      gene_conns = dplyr::left_join(gene_conns, peak_info, by = "peak") %>% 
        select(peak, connectedGene, p2gCorr, seqnames, start, end,
               nearestGene, peakType, distToTSS, nearestTSS, 
               GC, idx)
    }
    else{
      #message("gene has no connections")
      gene_conns = data.frame(peak = NA, connectedGene = gene, p2gCorr = NA, 
                              seqnames =NA , start = NA, end = NA, 
                              nearestGene = NA, peakType = NA, 
                              distToTSS = NA, nearestTSS = NA,  
                              GC = NA, idx =NA)
      }
      return(gene_conns)
  }
}
############################
# apply above function to list of genes
genes = unique(gene_model_res$gene_short_name)

gene_conns = data.frame(peak = character(), connectedGene = character(), p2gCorr = numeric(), 
                        seqnames = character(), start = integer(), end = integer(), 
                        nearestGene = character(), peakType = character(), 
                        distToTSS = integer(), nearestTSS = character(),  
                        GC = numeric(), idx =integer())
for(gene in genes){
  gc = getGeneConns(gene)
  gene_conns = rbind(gene_conns, gc)
}

head(gene_conns)
############################
# add model results for each condition
gene_conns_fin = filter(gene_conns, !is.na(peak))
# SAHA
SAHA_deg_res = filter(gene_model_res, term == "dose_SAHA") %>% 
  select(gene_short_name, estimate, q_value) %>% 
  rename(SAHA_gene_estimate = estimate, SAHA_gene_q_val = q_value)
  
SAHA_daPeak_res = filter(peak_model_res, term == "dose_SAHA") %>% 
  select(peak, estimate, q_value) %>% 
  rename(SAHA_peak_estimate = estimate, SAHA_peak_q_val = q_value)

# Dex
Dex_deg_res = filter(gene_model_res, term == "dose_Dex") %>% 
  select(gene_short_name, estimate, q_value) %>% 
  rename(Dex_gene_estimate = estimate, Dex_gene_q_val = q_value)

Dex_daPeak_res = filter(peak_model_res, term == "dose_Dex") %>% 
  select(peak, estimate, q_value) %>% 
  rename(Dex_peak_estimate = estimate, Dex_peak_q_val = q_value)

# Nutlin
Nutlin_deg_res = filter(gene_model_res, term == "dose_Nutlin") %>% 
  select(gene_short_name, estimate, q_value) %>% 
  rename(Nutlin_gene_estimate = estimate, Nutlin_gene_q_val = q_value)

Nutlin_daPeak_res = filter(peak_model_res, term == "dose_Nutlin") %>% 
  select(peak, estimate, q_value) %>% 
  rename(Nutlin_peak_estimate = estimate, Nutlin_peak_q_val = q_value)

# BMS
BMS_deg_res = filter(gene_model_res, term == "dose_BMS") %>% 
  select(gene_short_name, estimate, q_value) %>% 
  rename(BMS_gene_estimate = estimate, BMS_gene_q_val = q_value)

BMS_daPeak_res = filter(peak_model_res, term == "dose_BMS") %>% 
  select(peak, estimate, q_value) %>% 
  rename(BMS_peak_estimate = estimate, BMS_peak_q_val = q_value)

# add all info to gene connection df
gene_conns_final = left_join(gene_conns_fin, SAHA_deg_res, by = c("connectedGene" = "gene_short_name"))
gene_conns_final = left_join(gene_conns_final, SAHA_daPeak_res, by = "peak")

gene_conns_final = left_join(gene_conns_final, Dex_deg_res, by = c("connectedGene" = "gene_short_name"))
gene_conns_final = left_join(gene_conns_final, Dex_daPeak_res, by = "peak")

gene_conns_final = left_join(gene_conns_final, Nutlin_deg_res, by = c("connectedGene" = "gene_short_name"))
gene_conns_final = left_join(gene_conns_final, Nutlin_daPeak_res, by = "peak")

gene_conns_final = left_join(gene_conns_final, BMS_deg_res, by = c("connectedGene" = "gene_short_name"))
gene_conns_final = left_join(gene_conns_final, BMS_daPeak_res, by = "peak")

write.table(gene_conns_final, file = paste0(out_dir, "Plots/NB9/g2p_conns_with_DEG_DAS_dat.txt"), 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)



#################################################################################
#################################################################################
# 6/26/2020 start here

# analyze relationship between sig changed genes and connected peaks. 
gene_conns_f = read.table(file = paste0(out_dir, "Plots/NB9/g2p_conns_with_DEG_DAS_dat.txt"), 
                              head = TRUE)


# Summarize and Reformat data for ggplotting 
count_conns = group_by(gene_conns_f, connectedGene) %>% 
  mutate(tot_conn_peaks = n(),
         tot_conn_peaks_pr = length(peakType[peakType == "Promoter"]),
         tot_conn_peaks_distal = length(peakType[peakType == "Distal"]),
         tot_conn_peaks_intronic = length(peakType[peakType == "Intronic"]),
         tot_conn_peaks_exonic = length(peakType[peakType == "Exonic"]),
         gene_direction_SAHA = ifelse(SAHA_gene_q_val <= 0.05 & SAHA_gene_estimate > 0, "up",
                                      ifelse(SAHA_gene_q_val <= 0.05 & SAHA_gene_estimate < 0, "down", "un")),
         gene_direction_Dex = ifelse(Dex_gene_q_val <= 0.05 & Dex_gene_estimate > 0, "up",
                                      ifelse(Dex_gene_q_val <= 0.05 & Dex_gene_estimate < 0, "down", "un")),
         gene_direction_BMS = ifelse(BMS_gene_q_val <= 0.05 & BMS_gene_estimate > 0, "up",
                                      ifelse(BMS_gene_q_val <= 0.05 & BMS_gene_estimate < 0, "down", "un")),
         gene_direction_Nutlin = ifelse(Nutlin_gene_q_val <= 0.05 & Nutlin_gene_estimate > 0, "up",
                                      ifelse(Nutlin_gene_q_val <= 0.05 & Nutlin_gene_estimate < 0, "down", "un")),
         SAHA_tot_sig_peaks = length(SAHA_peak_q_val[!is.na(SAHA_peak_q_val) & SAHA_peak_q_val <= 0.05]),
         SAHA_sig_peaks_up = length(SAHA_peak_q_val[!is.na(SAHA_peak_q_val) & SAHA_peak_q_val <= 0.05 & !is.na(SAHA_peak_estimate) & SAHA_peak_estimate > 0]),
         SAHA_sig_peaks_dn = length(SAHA_peak_q_val[!is.na(SAHA_peak_q_val) & SAHA_peak_q_val <= 0.05 & !is.na(SAHA_peak_estimate) & SAHA_peak_estimate < 0]),
         Dex_tot_sig_peaks = length(Dex_peak_q_val[!is.na(Dex_peak_q_val) & Dex_peak_q_val <= 0.05]),
         Dex_sig_peaks_up = length(Dex_peak_q_val[!is.na(Dex_peak_q_val) & Dex_peak_q_val <= 0.05 & !is.na(Dex_peak_estimate) & Dex_peak_estimate > 0]),
         Dex_sig_peaks_dn = length(Dex_peak_q_val[!is.na(Dex_peak_q_val) & Dex_peak_q_val <= 0.05 & !is.na(Dex_peak_estimate) & Dex_peak_estimate < 0]),
         BMS_tot_sig_peaks = length(BMS_peak_q_val[!is.na(BMS_peak_q_val) & BMS_peak_q_val <= 0.05]),
         BMS_sig_peaks_up = length(BMS_peak_q_val[!is.na(BMS_peak_q_val) & BMS_peak_q_val <= 0.05 & !is.na(BMS_peak_estimate) & BMS_peak_estimate > 0]),
         BMS_sig_peaks_dn = length(BMS_peak_q_val[!is.na(BMS_peak_q_val) & BMS_peak_q_val <= 0.05 & !is.na(BMS_peak_estimate) & BMS_peak_estimate < 0]),
         Nutlin_tot_sig_peaks = length(Nutlin_peak_q_val[!is.na(Nutlin_peak_q_val) & Nutlin_peak_q_val <= 0.05]),
         Nutlin_sig_peaks_up = length(Nutlin_peak_q_val[!is.na(Nutlin_peak_q_val) & Nutlin_peak_q_val <= 0.05 & !is.na(Nutlin_peak_estimate) & Nutlin_peak_estimate > 0]),
         Nutlin_sig_peaks_dn = length(Nutlin_peak_q_val[!is.na(Nutlin_peak_q_val) & Nutlin_peak_q_val <= 0.05 & !is.na(Nutlin_peak_estimate) & Nutlin_peak_estimate < 0])) %>%
  ungroup()


count_conns_clean = distinct(count_conns, connectedGene, .keep_all = TRUE) %>% 
  select(connectedGene, tot_conn_peaks, tot_conn_peaks_pr, tot_conn_peaks_distal,
         tot_conn_peaks_intronic, tot_conn_peaks_exonic, gene_direction_SAHA,
         gene_direction_Dex, gene_direction_BMS, gene_direction_Nutlin,
         SAHA_tot_sig_peaks, SAHA_sig_peaks_up, SAHA_sig_peaks_dn, 
         Dex_tot_sig_peaks, Dex_sig_peaks_up, Dex_sig_peaks_dn,
         BMS_tot_sig_peaks, BMS_sig_peaks_up, BMS_sig_peaks_dn,
         Nutlin_tot_sig_peaks, Nutlin_sig_peaks_up, Nutlin_sig_peaks_dn)


cc_melt_gene_dir = melt(count_conns_clean, 
                        id.vars= c("connectedGene", "tot_conn_peaks", "tot_conn_peaks_pr", "tot_conn_peaks_distal",
                                   "tot_conn_peaks_intronic", "tot_conn_peaks_exonic"), 
                        measure.vars = c("gene_direction_SAHA", "gene_direction_Dex", 
                                         "gene_direction_BMS", "gene_direction_Nutlin"), 
                        variable.name = "gene_direction_drug", value.name = "direction")
cc_melt_gene_dir$order_idx = row.names(cc_melt_gene_dir)

cc_melt_tot_sig_peaks = melt(count_conns_clean, 
                        id.vars= c("connectedGene"), 
                        measure.vars = c("SAHA_tot_sig_peaks", "Dex_tot_sig_peaks", 
                                         "BMS_tot_sig_peaks", "Nutlin_tot_sig_peaks"), 
                        variable.name = "tot_sig_peaks_drug", value.name = "n_sigPeaks")
cc_melt_tot_sig_peaks$order_idx = row.names(cc_melt_tot_sig_peaks)

cc_melt_sig_peaks_up = melt(count_conns_clean, 
                             id.vars= c("connectedGene"), 
                             measure.vars = c("SAHA_sig_peaks_up", "Dex_sig_peaks_up", 
                                              "BMS_sig_peaks_up", "Nutlin_sig_peaks_up"), 
                             variable.name = "sig_peaks_up_drug", value.name = "n_Peaks_up")
cc_melt_sig_peaks_up$order_idx = row.names(cc_melt_sig_peaks_up)

cc_melt_sig_peaks_dn = melt(count_conns_clean, 
                            id.vars= c("connectedGene"), 
                            measure.vars = c("SAHA_sig_peaks_dn", "Dex_sig_peaks_dn", 
                                             "BMS_sig_peaks_dn", "Nutlin_sig_peaks_dn"), 
                            variable.name = "sig_peaks_down_drug", value.name = "n_Peaks_dn")
cc_melt_sig_peaks_dn$order_idx = row.names(cc_melt_sig_peaks_dn)

# merge separately melted dfs 
final_merged_conn_data = inner_join(cc_melt_gene_dir, cc_melt_tot_sig_peaks, by = "order_idx")
final_merged_conn_data = inner_join(final_merged_conn_data, cc_melt_sig_peaks_up, by = "order_idx")
final_merged_conn_data = inner_join(final_merged_conn_data, cc_melt_sig_peaks_dn, by = "order_idx")

# clean up 
final_merged_conn_data = mutate(final_merged_conn_data, 
                                drug = ifelse(grepl("SAHA", gene_direction_drug), "SAHA",
                                              ifelse(grepl("Dex", gene_direction_drug), "Dex",
                                                     ifelse(grepl("Nutlin", gene_direction_drug), "Nutlin", "BMS")))) %>% 
                                  select(-c(connectedGene.y, connectedGene.x.x, connectedGene.y.y,
                                            gene_direction_drug, tot_sig_peaks_drug, sig_peaks_up_drug,
                                            sig_peaks_down_drug))

# eliminate NAs in gene directions
final_merged_conn_data = mutate(final_merged_conn_data, 
                                direction = ifelse(is.na(direction), "un", direction))

#####################################
# plots describing sig peak-gene connections

# number of peaks connected to up, down, unchanged genes 
pdf(paste0(out_dir,"Plots/NB9/Boxplot_Conns_per_Gene.pdf"), width = 5, height = 5)
ggplot(final_merged_conn_data, aes(x = direction, y = log10(tot_conn_peaks), fill = direction)) +
  geom_violin() + 
  geom_boxplot(width=0.1) + 
  facet_wrap(~drug) + 
  xlab("Drugged Gene Expression") + 
  ylab("log10(n connected peaks)")
dev.off()

# proportion of peak types connected to up, down, unchanged genes 

x = group_by(final_merged_conn_data, drug, direction) %>% 
  summarize(connected_pr = sum(tot_conn_peaks_pr),
            connected_intronic = sum(tot_conn_peaks_intronic), 
            connected_exonic = sum(tot_conn_peaks_exonic), 
            connected_distal = sum(tot_conn_peaks_distal))

# number of sig peaks connected to up, down, unchanged genes 
pdf(paste0(out_dir,"Plots/NB9/Boxplot_sigPeak_Conns_per_Gene.pdf"), width = 5, height = 5)
ggplot(final_merged_conn_data, aes(x = direction, y = log10(n_sigPeaks), fill = direction)) +
  geom_violin() + 
  geom_boxplot(width=0.1) + 
  facet_wrap(~drug) + 
  xlab("Drugged Gene Expression") + 
  ylab("log10(n sig changed peaks)")
dev.off()


# fraction of peaks sig up/down/unchanged connected to up, down, unchanged genes 
plot_df_1 = group_by(final_merged_conn_data, drug, direction) %>% 
  mutate(n_Peaks_un = tot_conn_peaks - n_sigPeaks) %>% 
  summarize(connected_up = sum(n_Peaks_up),
            connected_unch = sum(n_Peaks_un), 
            connected_down = sum(n_Peaks_dn)) %>% 
  melt(id.vars= c("drug", "direction"), 
       measure.vars = c("connected_up", "connected_unch", 
                        "connected_down"), 
       variable.name = "peak_direction", value.name = "nPeaks")
pdf(paste0(out_dir,"Plots/NB9/BarPlot_FracPeakDirections.pdf"), width = 5, height = 5)
ggplot(plot_df) +
  geom_col(aes(x = direction, y = nPeaks, fill = peak_direction), position = "fill") + 
  facet_wrap(~drug) + 
  scale_fill_viridis(discrete = T) +
  xlab("Drugged Gene Expression") + 
  ylab("fraction connected peaks")
dev.off()

# fraction of up/down/unchanged genes connected to at least one sig peak
plot_df_2 = mutate(final_merged_conn_data, 
                 has_sigPeak_conn = ifelse(n_sigPeaks > 0,  1, 0),
                 has_NosigPeak_conn = ifelse(n_sigPeaks == 0,  1, 0),
                 has_sigPeakUp_conn = ifelse(n_Peaks_up > 0,  1, 0),
                 has_sigPeakDown_conn = ifelse(n_Peaks_dn > 0,  1, 0)) %>% 
  group_by(drug, direction) %>% 
  summarise( genes_withNoSig_conn = sum(has_NosigPeak_conn),
             genes_withSigUp_conn = sum(has_sigPeakUp_conn), 
             genes_withSigDown_conn = sum(has_sigPeakDown_conn)) %>% 
  melt(id.vars= c("drug", "direction"), 
       measure.vars = c("genes_withNoSig_conn", "genes_withSigUp_conn", 
                        "genes_withSigDown_conn"), 
       variable.name = "gene_class", value.name = "nGenes")
# bars represented as fraction
pdf(paste0(out_dir,"Plots/NB9/BarPlot_FracSigGenes_with_sigConns.pdf"), width = 5, height = 5)
ggplot(plot_df_2) +
  geom_col(aes(x = direction, y = nGenes, fill = gene_class), position = "fill") + 
  facet_wrap(~drug) + 
  scale_fill_viridis(discrete = T) +
  xlab("Drugged Gene Expression") + 
  ylab("fraction connected peaks")
dev.off()   

# bars actual gene numbers
pdf(paste0(out_dir,"Plots/NB9/BarPlot_SigGenes_with_sigConns.pdf"), width = 5, height = 5)
ggplot(plot_df_2) +
  geom_col(aes(x = direction, y = nGenes, fill = gene_class)) + 
  facet_wrap(~drug) + 
  scale_fill_viridis(discrete = T) +
  xlab("Drugged Gene Expression") + 
  ylab("fraction connected peaks")
dev.off()     

###################################################
# Metrics from perspective of peaks. 
# i.e. what fraction of sig changed peaks connected to sig changed genes? 
peak_conns_df =  distinct(count_conns, peak, .keep_all = TRUE) %>% 
  select(peak, connectedGene, gene_direction_SAHA, gene_direction_Dex, 
         gene_direction_BMS, gene_direction_Nutlin, 
         SAHA_peak_estimate, SAHA_peak_q_val,
         Dex_peak_estimate, Dex_peak_q_val,
         BMS_peak_estimate,  BMS_peak_q_val,
         Nutlin_peak_estimate, Nutlin_peak_q_val) %>% 
  melt(id.vars= c("peak", "connectedGene"), 
       measure.vars = c("genes_withNoSig_conn", "genes_withSigUp_conn", 
                        "genes_withSigDown_conn"), 
       variable.name = "gene_class", value.name = "nGenes")

######################################################
######################################################
# 07/22/2020

# Simplified peak-gene analysis: 
# 1) isolate DA peaks 
# 2) identify any DA peaks which are promoters 
# 3) Are DA promoters concordant DEGs in same drug context?
# 4) identify any connected sites which include a DA site. 
# 5) Assess directionality and significance of transcripts associated with connected promoter

#########
# load DEG res
gene_model_res = read.table(file = paste0(out_dir, "Plots/NB8/sciPlex_RNA_A549_small_screen_DEGs.txt"), head = T)
# load DA site res: 
peak_model_res = read.table(file = paste0(out_dir, "Plots/NB8/sciPlex_RNA_A549_small_screen_DApeaks.txt"), head = T)
# load full peak set
peakSet = data.frame(getPeakSet(prj)) %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end), 
         gene_short_name = ifelse(peakType == "Promoter" & distToGeneStart < 500, nearestGene, "NA")) %>% 
  select(distToGeneStart, nearestGene, peakType, distToTSS, nearestTSS,  
         GC, gene_short_name, peak)

# 1) for drug, isolate DA peaks and DEGs
treatment = "dose_SAHA"
DA_sites = filter(peak_model_res, grepl("dose_SAHA", term)) %>%
  filter(q_value < 0.05)
DEGs = filter(gene_model_res, grepl("dose_SAHA", term)) %>%
  filter(q_value < 0.05)

# get DA promoters only. 
DA_proms = filter(DA_sites, peakType == "Promoter", !is.na(gene_short_name)) %>% 
  select(distToGeneStart, nearestGene, peakType, distToTSS, nearestTSS, 
         gene_short_name, GC, idx, peak, num_cells_expressed, term, estimate,
         std_err, q_value)
cat("number of DA promoter peaks = ", dim(DA_proms)[1], "\n")

# 2) isloate sig promoters (ATAC) with sig DEG (RNA)
DA_proms_DEGs = inner_join(DA_proms, DEGs, by = "gene_short_name")
cat("number of DA prometers for DE genes = ", dim(DA_proms_DEGs)[1], "\n")

# 3) Are DA promoters concordant DEGs in same drug context?
DA_proms_DEGs$concordant = sign(DA_proms_DEGs$estimate.x) == sign(DA_proms_DEGs$estimate.y)
cat("number of concordant DA prometers and DE genes = ", dim(DA_proms_DEGs[DA_proms_DEGs$concordant == TRUE,])[1], "\n")

# 4) identify any connected sites which include a DA site. 
DA_peaks  = DA_sites %>% 
  select(distToGeneStart, nearestGene, peakType, distToTSS, nearestTSS, 
         gene_short_name, GC, idx, peak, num_cells_expressed, term, estimate,
         std_err, q_value)

Pconns <- getCoAccessibility(
  ArchRProj = prj,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
Pconns_df = data.frame(Pconns)
peaks = data.frame(Pconns@metadata$peakSet) %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end))
peaks$idx = as.numeric(row.names(peaks))

# simplify peak index map
peak_idx_map = select(peaks, peak, idx)

# add query peaks to connection df
Pconns_df_1 = inner_join(Pconns_df, peak_idx_map, by = c("queryHits" = "idx"))
Pconns_df_1 = rename(Pconns_df_1, q_peak = peak)

# add subject peaks to connection df
Pconns_df_2 = inner_join(Pconns_df_1, peak_idx_map, by = c("subjectHits" = "idx"))
Pconns_df_2 = rename(Pconns_df_2, s_peak = peak)

# get any connection involving a DA site 
Pconns_DA = dplyr::filter(Pconns_df_2, q_peak %in% DA_sites$peak)

# get any DA site connected to a promoter 
promoter_peaks = filter(peakSet, !(gene_short_name=="NA"))
#
Pconns_DA_pr = inner_join(Pconns_DA, promoter_peaks, by = c("s_peak" = "peak"))


#########################################################
#########################################################
# correlate Drug dose coefficients from Transcript abundance
# and Gene scores at DE genes 

# get all sig genes for all drugs
DEgene_model_res = 
  gene_model_res %>% 
  filter(grepl("dose", term)) %>% 
  select(id, gene_short_name, num_cells_expressed, 
         term, estimate, p_value, q_value) %>% 
  mutate(gene_and_term = paste(gene_short_name, term, sep ="_")) %>% 
  filter(q_value < 0.05)

# isolate gene score estimates for joining to DEGene info
gs_model_res = 
  genescore_model_res %>% 
  filter(grepl("dose", term)) %>% 
  mutate(gene_and_term = paste(gene_short_name, term, sep ="_")) %>% 
  select(gene_and_term, estimate, q_value) %>% 
  rename(estimate_gs = estimate, q_value_gs = q_value)
  
RNAvsGenescore_estim = 
  dplyr::left_join(DEgene_model_res, gs_model_res, by = "gene_and_term")

## scatterplots of drug specific coefficients
ggplot(RNAvsGenescore_estim, aes(x =  estimate, y = estimate_gs)) +
  geom_point_rast(size = 2, alpha = 0.5) +
  facet_wrap(~term) +
  xlab("RNA Coefficient") +
  ylab("ATAC Coefficient") +
  stat_cor(method = "pearson", aes(label = ..r.label..),
           label.x = -0.5, label.y = 0.3) +
  geom_smooth(method=lm, se = TRUE) + 
  theme_bw() 
ggsave(filename = paste0(out_dir,"Plots/NB8/CoeficientScatter_RNAvsGeneScore_DEGs.pdf"),
       width = 4, height = 4)







