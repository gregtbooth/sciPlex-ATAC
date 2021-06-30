
# 04/15/21 
#### UNFINISHED #######

basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/")
dir.create(paste0(out_dir, "results/NB12"))
setwd(paste0(out_dir, "results/NB12"))

suppressPackageStartupMessages({
  library(ArchR)
  library(Seurat)
  library(monocle3)
  library(cicero)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
  library(Gviz)
  library(biomaRt)
  library(RColorBrewer)
  library(viridis)
})

options(dplyr.summarise.inform = FALSE) 

set.seed(2017)

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

prj = loadArchRProject(path = paste0(out_dir, "mlr_filtered_annotated"), 
                       showLogo = FALSE) 


##############################################
# 1) Get Cicero connections
##############################################
# only need to run once (output saved)
# make a peak x cell cds object for running cicero. 

# load PeakMatrix into memory
pMat = getMatrixFromArrow(
  ArrowFile = paste0(out_dir, "mlr_filtered_annotated/ArrowFiles/MLR.arrow"),
  ArchRProj = prj, 
  cellNames = prj$cellNames,
  useMatrix = "PeakMatrix",
  binarize = TRUE)

# Create CDS from peak Matrix (SummarizedExperiment)
rd = data.frame(rowRanges(pMat)) %>% 
  tibble::rownames_to_column() %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end))

# peaks in this list are in a DIFFERENT ORDER 
# relative to pMat (rowRanges above is correct)
pd = data.frame(getPeakSet(prj)) %>% 
  tibble::rownames_to_column() %>% 
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

cds_a = monocle3::new_cell_data_set(assays(pMat)$PeakMatrix, 
                                    cell_metadata = colData(pMat), 
                                    gene_metadata = a_rd)

#find connections from ALL cells 
# preprocess cds for cicero: 
message("Getting Cicero connections using All cells\n")

cds_archr = cds_a

set.seed(2017)
cds_archr <- detect_genes(cds_archr)
cds_archr = cds_archr[rowData(cds_archr)$num_cells_expressed > 0,]
cds_archr <- estimate_size_factors(cds_archr)
cds_archr = preprocess_cds(cds_archr, method = "LSI", num_dimensions=50)
cds_archr = reduce_dimension(cds_archr, reduction_method = 'UMAP', preprocess_method = "LSI")

umap_coords <- reducedDims(cds_archr)$UMAP
cicero_cds <- make_cicero_cds(cds_archr, reduced_coordinates = umap_coords)

# run cicero on preprocessed data 
data("human.hg19.genome")
conns <- run_cicero(cicero_cds, 
                    human.hg19.genome, 
                    window = 5e+05,
                    silent = FALSE,
                    sample_num = 100)

filter(conns, coaccess > 0) %>% 
  write.table(file = "cicero_peak_conns.txt", 
              sep = "\t",
              quote = FALSE,
              row.names = FALSE, 
              col.names = TRUE)

##############################################
##############################################
# Identify promoter peaks connected to sites which 
# differ between Allo and Bead stimulated T-cells

peak_type = data.frame(getPeakSet(prj)) %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end)) %>% 
  dplyr::select(peak, peakType, distToGeneStart)

# load DA coefs. 
coefs_PT = read.table(paste0(out_dir, "results/NB10/Tcell_Pseudotime_DApeak_coefs_directions.txt"), head = TRUE) %>% 
  filter(term == "Pseudotime") %>% 
  left_join(peak_type, by = "peak") %>% 
  mutate(gene_short_name = ifelse(peakType == "Promoter" & distToGeneStart < 500, nearestGene, "NA")) 

# collect data on all promoter peaks
coefs_promoters = filter(coefs_PT, gene_short_name != "NA")
# filter for DA promoters
coefs_sig_promoters = filter(coefs_promoters, p_value <= 0.05)


# find sig DA sites connected to promoters
cons.info = 
  read.table("cicero_peak_conns.txt", head = TRUE) %>% 
  mutate(p1_center = (as.numeric(stringr::str_split_fixed(Peak1, pattern = "_", 3)[,2]) + 
                        as.numeric(stringr::str_split_fixed(Peak1, pattern = "_", 3)[,3]))/2, 
         p2_center = (as.numeric(stringr::str_split_fixed(Peak2, pattern = "_", 3)[,2]) + 
                        as.numeric(stringr::str_split_fixed(Peak2, pattern = "_", 3)[,3]))/2, 
         dist = abs(p2_center - p1_center)) 

# filter connections based on threshold
cicero_threshold = 0.1

cons.info_f = 
  cons.info[, c("Peak1", "Peak2", "dist", "coaccess")] %>% 
  filter(!is.na(dist), dist >= 1000, !is.na(coaccess), coaccess > cicero_threshold)

peak_info= data.frame(getPeakSet(prj)) %>% 
  mutate(peak = paste0(seqnames, "_", start, "_", end), 
         gene_short_name = ifelse(peakType == "Promoter" & distToGeneStart < 500, nearestGene, "NA")) %>% 
  dplyr::select(peak, distToGeneStart,
                nearestGene, peakType, distToGeneStart, nearestTSS,
                gene_short_name, GC)

# isolate connections to/from significant sites and attach peak_info
DA_coefs = filter(coefs_PT, p_value <= 0.05) 
DA_site_conns = cons.info_f %>% 
  left_join(DA_coefs, by = c("Peak1" = "peak")) %>%  #attach info for peak1
  left_join(DA_coefs, by = c("Peak2" = "peak")) %>% # attach info for peak2
  dplyr::filter(!is.na(term.x)) %>% 
  dplyr::select(Peak1, Peak2, dist, coaccess, num_cells_expressed, 
         term, estimate, p_value, q_value, direction, peakType_p1 = peakType.x, distToGeneStart_p1 = distToGeneStart.x,
         gene_short_name_p1 = gene_short_name.x, peakType_p2 = peakType.y, distToGeneStart_p2 = distToGeneStart.y,
         gene_short_name_p2 = gene_short_name.y)

# filter for Sig sites which are either promoters themselves, or connected to promoters. 
DA_pr_conns = filter(DA_site_conns, gene_short_name_p1 != "NA" | gene_short_name_p2 != "NA")

# isolate non promoter DA sites that are connected to promoters (possibly non_DA genes)
DA_pr_conns_dist = filter(DA_pr_conns, gene_short_name_p1 == "NA") %>% 
  dplyr::select(gene_short_name = gene_short_name_p2, peakType = peakType_p1, estimate, p_value, direction)






# isolate list of genes (peak1 or peak2 based) for GSEA 
p1_genes = dplyr::select(DA_pr_conns, gene = gene_short_name_p1, direction)
p2_genes = dplyr::select(DA_pr_conns, gene = gene_short_name_p2, direction)

genes = rbind(p1_genes, p2_genes) %>% 
  dplyr::filter(gene != "NA") %>% 
  group_by(gene, direction) %>% 
  summarise(n_conn_dir = n())

genes_dir_spread = spread(genes, direction, n_conn_dir, fill = 0) %>% 
  tibble::column_to_rownames(var = "gene")

# Find the directional bias, meaning: of promoter connected DA sites, what fraction are changing in the same direction)  
genes_dir_spread$bias<-apply(X=data.frame(genes_dir_spread), MARGIN=1, FUN=function(x){max(x)/sum(x)})
genes_dir_spread$max_dir = colnames(genes_dir_spread)[max.col(genes_dir_spread, ties.method = "first")]


# filter for genes that are confidently changing in only one direction: 
genes_dir_final = dplyr::filter(genes_dir_spread, bias > 0.5) %>% 
  tibble::rownames_to_column(var = "gene") %>% 
  mutate(nconns = Closing + Opening + Static +Transient) %>% 
  dplyr::select(gene, nconns, bias, max_dir)

# run Gene Set Enrichment analysis: 



  
                       
                       
                       
                       
                       
                       
        




