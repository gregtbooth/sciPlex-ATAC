# This bit of code compares our sciplex-ATAC peaks with GR binding sites. 
# GR binding sites we're defined by Vockley et. al. 2016 using ChIP seq. 

basepath = "github/"
bin_directory = paste0(basepath, "analysis/bin/")
ref_path = "genomes/HG19/annotations/downloaded_data/"
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB13"))
setwd(paste0(out_dir, "results/NB13"))

suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(cicero)
  library(ggplot2)
  library(ggridges)
  library(tidymodels)
  library(ggrastr)
  library(ggpubr)
  library(glmnet)
})

# load processed ArchR project
prj = loadArchRProject(path = paste0(out_dir, "sc2_filtered/"),
                       showLogo = FALSE)

# create cds from peak x cell matrix
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
  select(rowname = rowname.x, idx = idx.x, score, replicateScoreQuantile, groupScoreQuantile, 
         Reproducibility, GroupReplicate, distToGeneStart,
         nearestGene, peakType, distToTSS, nearestTSS,
         gene_short_name, GC, peak)

row.names(a_rd) <- a_rd$peak
row.names(pMat) <- a_rd$peak
cds_archr = monocle3::new_cell_data_set(assays(pMat)$PeakMatrix, 
                                       cell_metadata = colData(pMat), 
                                       gene_metadata = a_rd)

##########################################
# identify peaks which overlap with GR binding sites 

grbs = read.table(file = paste0(ref_path, "GR_BindingSites.bed"), head = TRUE) %>% 
  mutate(type = "GR_ChIP") %>% 
  select(Chromosome, Start, end, type, merged_reads_per_peak)

cds_archr = annotate_cds_by_site(cds_archr, grbs)

peak_dat = data.frame(rowData(cds_archr)) %>% 
  select(site_name = peak, peakType, gene_short_name, type)

########################################
# add Dex-dose response info (from sciPlex ATAC) 
# for promoters peaks, add Dex-dose response info (from sciPlex RNA) 
d = "Dex"

DA.res = read.table(paste0(out_dir, "results/NB5/sciPlex_ATAC_A549_small_screen_DApeaks.txt"), head = T)

RNA.res = read.table(file = paste0(out_dir, "results/NB5/sciPlex_RNA_A549_small_screen_DEGs.txt"), head = T) %>% 
  dplyr::filter(term == paste0("dose_", d)) %>% 
  dplyr::select(tss.id = id, gene_short_name, term, expr.estimate = estimate, pval = p_value, 
                qval = q_value, num.cells.expr = num_cells_expressed) 

drug_effects = filter(DA.res, term == paste0("dose_", d)) %>% 
  rename(site_name = peak) %>% 
  mutate(gene_short_name = ifelse(peakType == "Promoter" & distToTSS < 500, 
                                  nearestGene, 'NA'), 
         peak.type = ifelse(peakType == "Promoter" & distToTSS < 500, 
                            "Promoter", 'Distal'), 
         promoter_gene = gene_short_name, 
         direction = ifelse(q_value < 0.05 & estimate > 0, "Opening", 
                            ifelse(q_value < 0.05 & estimate < 0, "Closing", "Static"))) %>% 
  left_join(RNA.res, by = "gene_short_name") %>% 
  select( c(nearestGene, peak.type, distToTSS, 
            nearestTSS, GC, site_name, gene_short_name, 
            promoter_gene, tss.id, direction, 
            atac.estimate = estimate, atac.qval = q_value, 
            rna.estimate = expr.estimate, rna.qval = qval)) %>% 
  distinct(site_name, .keep_all = TRUE)


# final peak data ONLY includes peaks present within > 1% of dex-treated or vehicle cells
pDat = left_join(select(peak_dat, -gene_short_name), drug_effects, by = "site_name") %>% 
  select(-peakType, -promoter_gene) %>% 
  filter(!is.na(direction)) # only sites in > 1% of dex-treated cells were evaluated in DA analysis

########################################
# Collect stats for GR binding sites

# fraction of GRBS in my peaks. 
tot_grbs = dim(grbs)[1]
grbs_peaks = filter(pDat, !is.na(type)) %>% 
  summarise(grbs_peaks = n())
cat(paste0("total number of grbs = ", tot_grbs, "\n"))
cat(paste0("number of grbs contained in peaks found in >1% of dex/vehicle cells = ", grbs_peaks, "\n"))
cat(paste0("fraction of grbs contained in peaks = ", grbs_peaks/tot_grbs, "\n"))

# number of opening closing and static sites with or without GR peaks 
no_grbs_trends = filter(pDat, is.na(type)) %>% 
  group_by(direction) %>% 
  summarise(no_grbs = n())
yes_grbs_trends = filter(pDat, !is.na(type)) %>% 
  group_by(direction) %>% 
  summarise(yes_grbs = n())

trend_summary = left_join(no_grbs_trends, yes_grbs_trends, by = "direction") %>% 
  mutate(percent_grbs = (yes_grbs/(yes_grbs+no_grbs))*100)


# make stacked barplot of grbs containing sites
m_trend_summary = reshape2::melt(data = trend_summary, id.vars = "direction", 
                                 measure.vars = c("no_grbs", "yes_grbs")) %>% 
  ggplot() +
  geom_bar(aes(x = direction, y = value, fill = variable ), 
         color = "black", size = .25, 
         position="fill", stat = "identity") +
  #scale_fill_manual("sites", values = clr_values) +
  monocle3:::monocle_theme_opts() +
  theme(legend.position = "right", 
        text = element_text(size = 6),        
        legend.key.width = unit(0.15,"line"), 
        legend.key.height = unit(0.1,"line")) +  
  xlab("Dex-Affected Sites") +
  ylab("proportion GRBS") +
  ggsave("Stacked_Barplot_proportion_GRBS.pdf",
         width = 2.25, height = 1.5, unit = "in", dpi = 1200)

#######################################
# how many dex responsive gene promoters are connected 
# (through cicero) to GR binding sites?

conns = read.table(paste0(out_dir, "results/NB10/cicero_peak_conns.txt"), head = TRUE)

get_tss_conn_stats = function(pDat, cicero_conns, conn_threshold = 0.1){
  # cicero connections 
  cons.info = cicero_conns %>% 
      mutate(p1_center = (as.numeric(stringr::str_split_fixed(Peak1, pattern = "_", 3)[,2]) + 
                          as.numeric(stringr::str_split_fixed(Peak1, pattern = "_", 3)[,3]))/2, 
           p2_center = (as.numeric(stringr::str_split_fixed(Peak2, pattern = "_", 3)[,2]) + 
                          as.numeric(stringr::str_split_fixed(Peak2, pattern = "_", 3)[,3]))/2, 
          dist = abs(p2_center - p1_center))

  cons.info = select(cons.info, Peak1, Peak2, dist, coaccess)

  # set thresholds and define cicero-connections (graphical lasso)
  glasso.edges = filter(cons.info, dist >= 1000,  coaccess >= conn_threshold)

  ##########
  tss.distal.peaks = filter(pDat, !is.na(tss.id)) %>% 
    select(site_name, tss.id, atac.estimate, atac.qval, 
           direction, rna.estimate, rna.qval) %>%
    inner_join(glasso.edges, by = c("site_name" = "Peak1")) %>%
    group_by(tss.id, Peak2) %>% 
    summarize(coaccess = max(coaccess)) %>%
    # add connected distal peak info
    inner_join(select(pDat, site_name, d.type = type, d.atac.estimate = atac.estimate, 
                      d.atac.qval = atac.qval,  d.acc.trend = direction),
               by = c("Peak2" = "site_name")) 

  # for each promoter, summarize stats for distal peak connections 
  # including how many distal GR binding sites connected 
  tss.distal.peak.stats = group_by(tss.distal.peaks, tss.id) %>%
    summarize(
      n.dp = n(),
      n.o.dp = sum(d.acc.trend == "Opening"),
      n.c.dp = sum(d.acc.trend == "Closing"),
      n.grbs.dp = sum(d.type == "GR_ChIP"), # number of connected GR sites
      n.grbs.o.dp = sum(d.type == "GR_ChIP" & d.acc.trend == "Opening"), # number of connected GR sites that significantly open
      n.grbs.c.dp = sum(d.type == "GR_ChIP" & d.acc.trend == "Closing"), # number of connected GR sites that significantly close
      best.o.dp = max(ifelse(d.acc.trend == "Opening", coaccess, 0)),
      best.c.dp = max(ifelse(d.acc.trend == "Closing", coaccess, 0))) %>% 
    mutate(n.grbs.dp = ifelse(is.na(n.grbs.dp), 0, n.grbs.dp),
           n.grbs.o.dp = ifelse(is.na(n.grbs.o.dp), 0, n.grbs.o.dp),
           n.grbs.c.dp = ifelse(is.na(n.grbs.c.dp), 0, n.grbs.c.dp))

  # add each promoters DEG info to table
  tss.stats = data.frame(right_join(tss.distal.peak.stats, 
                                     select(filter(pDat, !is.na(tss.id)), 
                                            tss.id, gene_short_name, 
                                            rna.estimate, rna.qval), 
                                     by = "tss.id")) %>%
    distinct(tss.id, .keep_all = TRUE) %>% 
    mutate(rna.trend = ifelse(rna.qval <= 0.05 & rna.estimate > 0, "up",
                              ifelse(rna.qval <= 0.05 & rna.estimate < 0, "down", "no_effect")))

  # change any NAs to 0
  for (i in 2:9) {
    tss.stats[, i] = ifelse(is.na(tss.stats[, i]), 0, tss.stats[, i])
  }

  message("total number of promoters in analysis = ", dim(tss.stats)[1], "\n")
  message("total number of dex affected genes in analysis = ", dim(tss.stats[tss.stats$rna.qval < 0.05,])[1], "\n")
  
  return(tss.stats)
}

t.stats = get_tss_conn_stats(pDat = pDat,  cicero_conns = conns, conn_threshold = 0.1)

# summarize number of promoters with distal connections 
summary.stats = t.stats %>% 
  group_by(rna.trend) %>% 
  summarize(tot_pr = n(),
            pr_with_conns = sum(n.dp > 0), 
            pr_with_GRBS_conn = sum(n.grbs.dp >0))


# histogram of distal connections per promoter
ggplot(t.stats) +
  geom_histogram(aes(x = n.dp)) + # y = ..density..)) +
  facet_wrap(~rna.trend) + 
  xlab("distal connections") +
  ylab("fraction of promoters") +
  #ylim(0, 0.2) +
  theme_bw() 
ggsave(filename ="histogram_distal_peak_connections.pdf", width = 7.5, height = 2.5)

##############################  

