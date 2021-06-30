#set working environment in R

basepath = "github/"
out_dir = paste0(basepath, "analysis/archr/sc2_integratedRNA_v2/")
dir.create(paste0(out_dir,"Plots/NB6/"))
setwd(out_dir)

# load the developer version (master has issues with coAccessibility functions)
library(devtools)
load_all("~/scripts/github/ArchR") 

suppressPackageStartupMessages({
  #library(ArchR)
  library(monocle3)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
})

set.seed(1989) # ensures reproducibility of previous random number generation

DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)

addArchRGenome("hg19")

# load ArchR project 
prj = loadArchRProject(path = out_dir,
                       showLogo = FALSE)


#############################################################################
# get co-accessibility
prj <- addCoAccessibility(
  ArchRProj = prj,
  reducedDims = "IterativeLSI"
)

# retreive co-accessibility:
# note that all connection data stored in metadata(prj@peakSet)$CoAccessibility
cA <- getCoAccessibility(
  ArchRProj = prj,
  corCutOff = 0.5,
  resolution = 1000,
  returnLoops = TRUE
)

cA[[1]]

# demo plot of connections
markerGenes  <- c("COPS7A")

#isolate SAHA treated cells
SAHAonly = prj[grepl("SAHA", prj@cellColData$treatment_hashbased),]

p <- plotBrowserTrack(
  ArchRProj = SAHAonly, 
  groupBy = "treatment_hashbased", 
  geneSymbol = markerGenes, 
  upstream = 10000,
  downstream = 10000,
  loops = getCoAccessibility(prj)
)

plotPDF(plotList = p, 
        name = "/NB6/COPS7A-Browser-CoAccessibility_SAHAonly.pdf", 
        ArchRProj = prj, 
        addDOC = FALSE, width = 5, height = 5)

######################################################################
# get peak to gene links
prj <- addPeak2GeneLinks(
  ArchRProj = prj,
  reducedDims = "IterativeLSI"
)

# retrieve peak to gene links:
# note that all connection data stored in metadata(prj@peakSet)$
p2g <- getPeak2GeneLinks(
  ArchRProj = prj,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)

head(p2g)

# summary plots of numbers of connected peaks/gene and genes/peak

peakSet = getPeakSet(prj)
connPeakInfo = data.frame(peakSet)[p2g$idxATAC,]
connGeneInfo = data.frame(p2g@metadata$geneSet)[p2g$idxRNA,]
connPeakInfo$connGene = connGeneInfo$name
connPeakInfo$idxATAC = p2g$idxATAC
connPeakInfo$idxRNA = p2g$idxRNA
connPeakInfo = mutate(connPeakInfo, 
                      peak = paste(seqnames, start, end, sep = "_"))
# peaks per gene
ppg = group_by(connPeakInfo, connGene) %>% 
  summarize(nConnPeaks = n())

pdf(paste0(out_dir,"Plots/NB6/Hist_peaksPerGene.pdf"), width = 3, height = 3)
ggplot(ppg, aes(x = nConnPeaks)) +
  geom_histogram() + 
  labs(title = "number connected peaks", 
       x = "connected peaks per gene", 
       y = "number of genes") 
dev.off()

# genes per peak
cpp = group_by(connPeakInfo, peak, peakType, nearestGene) %>% 
  summarize(nConnGenes = n())

pdf(paste0(out_dir,"Plots/NB6/Hist_GenesPerPeak_facetPeakType.pdf"), width = 3, height = 3)
ggplot(cpp, aes(x = nConnGenes)) +
  geom_histogram(bins = 10) + 
  facet_wrap(~peakType) + 
  labs(title = "number connected genes", 
       x = "connected genes per peak", 
       y = "number of peaks") 
dev.off()

# boxplots of genes per peak by peak type(promoter, distal etc.)
pdf(paste0(out_dir,"Plots/NB6/Boxplot_GenesPerPeak_byType.pdf"), width = 4, height = 4)
ggplot(cpp, aes(x = peakType, y = log10(nConnGenes), fill = peakType)) +
  geom_violin() + 
  geom_boxplot(width=0.1) + 
  labs(title = "number connected genes", 
       x = "peak type", 
       y = "log10(number of conn genes)")
dev.off()

# demo plot of connections
markerGenes  <- c("COPS7A")

#isolate SAHA treated cells
SAHAonly = prj[grepl("SAHA", prj@cellColData$treatment_hashbased),]

p3 <- plotBrowserTrack(
  ArchRProj = SAHAonly, 
  groupBy = "treatment_hashbased", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(prj)
)

plotPDF(plotList = p3, 
        name = "/NB6/COPS7A-Browser-peak2Gene_loops_SAHAonly.pdf", 
        ArchRProj = prj, 
        addDOC = FALSE, width = 5, height = 5)


#########################################################
# save project with connection info
saveArchRProject(ArchRProj = prj, 
                 load = FALSE)



###########################################
# compare with cicero
#library(cicero)
# reduce dimensions for cicero aggrigation
#cds_atac <- detect_genes(cds_atac)
#cds_atac = cds_atac[rowData(cds_atac)$num_cells_expressed > 0,]
#cds_p = align_cds(counts(cds_p), residual_model_formula_str = "~umi")
#cds_atac <- estimate_size_factors(cds_atac)
#cds_atac = preprocess_cds(cds_atac, method = "LSI", num_dimensions=50)
#cds_atac = reduce_dimension(cds_atac, reduction_method = 'UMAP', preprocess_method = "LSI")

# custom plot
#colData(cds_atac)$UMAP_1 <- cds_atac@reducedDims$UMAP[,1]
#colData(cds_atac)$UMAP_2 <- cds_atac@reducedDims$UMAP[,2]
#TCdat_pr = data.frame(colData(cds_atac))

#pdf(paste0(out_dir,"Plots/NB8/UMAP_preCiceroAgg.pdf"), width = 3, height = 2.5)
#ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = treatment)) +
#  geom_point(size = 0.4, stroke = 0) +
#  theme(legend.position = "right", text = element_text(size = 12),  
#        legend.key.width = unit(0.5,"line"), legend.key.height = unit(0.5,"line")) + 
#  xlab("Component 1") +
#  ylab("Component 2") +
#  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
#  monocle3:::monocle_theme_opts()
#dev.off()

# create cicero aggregate cds (based on reduced dims above)
#umap_coords <- reducedDims(cds_atac)$UMAP
#cicero_cds <- make_cicero_cds(cds_atac, reduced_coordinates = umap_coords)
#saveRDS(cicero_cds, file = "cicero_cds")

#cicero_cds = readRDS("cicero_cds")

#data(human.hg19.genome)
# identify cicero connections from agg cds. 
#cicero_conns <- run_cicero(cds = cicero_cds,
#                           genomic_coords = human.hg19.genome,
#                           window = 5e+05,
#                           silent = FALSE,
#                           sample_num = 100)
##########################################

