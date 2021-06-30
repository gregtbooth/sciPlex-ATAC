# 
#load the following modules 
#module load MACS/2.2.7.1 
#
# start R (4.0.0)

#set working environment in R
basepath = "github/"
out_dir =paste0(basepath, "analysis/archr_revised/")
dir.create(paste0(out_dir, "results/NB4"))
setwd(paste0(out_dir, "results/NB4"))

suppressPackageStartupMessages(c(
  library(ArchR),
  library(dplyr),
  library(Gviz),
  library(biomaRt)
))

addArchRGenome("hg19")
set.seed(1)

# load processed ArchR project from NB2
prj = loadArchRProject(path = paste0(out_dir,"sc2_filtered/"),
                       showLogo = FALSE)

###############################
# make pseduo bulk replicates for monocle3-defined clusters
prj = addGroupCoverages(
  ArchRProj = prj,
  groupBy=  "clusters_mon3",
  useLabels = TRUE,
  minCells= 50,
  maxCells = 500,
  minReplicates = 2,
  maxReplicates = 5,
  sampleRatio = 0.8,  
  threads =  getArchRThreads(),
  verbose =  TRUE
)

# call peaks. 
prj <- addReproduciblePeakSet(
  ArchRProj = prj, 
  groupBy = "clusters_mon3", 
  peakMethod = "Macs2",
  reproducibility = "2",
  peaksPerCell = 500,
  maxPeaks = 150000,
  minCells = 25,
  excludeChr = c("chrM", "chrY"), 
  extendSummits = 250,
  promoterRegion = c(2000, 100),
  genomeAnnotation = getGenomeAnnotation(prj),
  geneAnnotation = getGeneAnnotation(prj),
  plot = TRUE,
  threads = getArchRThreads(),
  force = TRUE
)

# add peak matrix to arrow file
prj <- addPeakMatrix(prj, 
                     binarize = FALSE, # don't binarize here (can be binarized when retrieved)
                     force = TRUE)

# add peak motif information
projHeme5 <- addMotifAnnotations(ArchRProj = prj, 
                                 motifSet = "cisbp", 
                                 name = "Motif")


# save output
saveArchRProject(ArchRProj = prj, 
                 load = FALSE)

######################################################


# save genome ranges and .bed files for called peaks
peaks = getPeakSet(prj)
dir.create("PeakCalls")
saveRDS(peaks, file=paste0("PeakCalls/ArchrPeakSet"))
# make bed version
peaks_df = data.frame(peaks) %>% 
  select(chr = seqnames, start, end, score, peakType)
write.table(peaks_df, file = "PeakCalls/ArchrPeaks.bed", 
            col.names = TRUE, row.names = FALSE, 
            quote = FALSE, sep = "\t")

#Plot FRIP values by dose
cd = data.frame(prj@cellColData) %>% 
  mutate(Relative_dose_factor = 
  factor(Relative_dose, levels = c("0", "0.1", "0.5", "1",
                                      "5", "10", "50", "100")),
  new_treatment_label =
  sapply(treatment, function(x) {
    if (grepl("BMS", x))
      return("BMS345541")
    if (grepl("Dex", x))
      return("Dex")
    if (grepl("Nutlin", x))
      return("Nutlin3A")
    if (grepl("SAHA", x))
      return("SAHA")
  }))

ggplot(cd) +
  geom_boxplot(aes(x =  Relative_dose_factor, y = FRIP)) +
  facet_wrap("~new_treatment_label") +
  monocle3:::monocle_theme_opts() +
  scale_color_brewer(palette='Set1') +
  xlab("Relative Dose") +
  ylab("FRIP") + 
  theme_bw()
ggsave(filename = paste0("FRIPbyDose_drugFacet.png", sep = ""),
       width = 4, height = 3)



###################################
# create pseudobulk bigwig files for 
# all treatment conditions (drug and dose)

prj$hash_condition = paste0(prj$treatment, "_", prj$dose)

getGroupBW(
  ArchRProj = prj,
  groupBy = "hash_condition",
  #normMethod = "ReadsInTSS", # creates a group scale factor = 10k/sum(reads in TSS) 
  normMethod = "nFrags", # creates a group scale factor = 10k/sum(nFrags)
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)


##########################################################
##########################################################

# Use GVIZ for plotting pseudobulk data from bigwigs

plot_pb_browsertracks = function(archrPrj = prj, drug = "SAHA", bw_list = SAHA_bw_list,
                                 chrom = "chr12", start_window = 6810000, 
                                 end_window = 6840000, locus = "locus", Ylim = c(0,20))
{
  dir.create("browser_tracks/")
  ## Load BW files to data tracks
  dose_0_Track <- DataTrack(bw_list$dose_0, type = "histogram", fill = "gray",
                            col.histogram = "gray", name = paste0(drug, "_0"), 
                            ylim = Ylim)
  dose_0.1_Track <- DataTrack(bw_list$dose_0.1, type = "histogram", fill = "#1D1147FF", 
                              col.histogram = "#1D1147FF", name = paste0(drug, "_0.1"), 
                              ylim = Ylim)
  dose_0.5_Track <- DataTrack(bw_list$dose_0.5, type = "histogram", fill = "#51127CFF",
                              col.histogram = "#51127CFF", name = paste0(drug, "_0.5"), 
                              ylim = Ylim)
  dose_1_Track <- DataTrack(bw_list$dose_1, type = "histogram", fill = "#822681FF",
                            col.histogram = "#822681FF", name = paste0(drug, "_1"), 
                            ylim = Ylim)
  dose_5_Track <- DataTrack(bw_list$dose_5, type = "histogram", fill = "#B63679FF", 
                            col.histogram = "#B63679FF", name = paste0(drug, "_5"), 
                            ylim = Ylim)
  dose_10_Track <- DataTrack(bw_list$dose_10, type = "histogram", fill = "#E65164FF", 
                             col.histogram = "#E65164FF", name = paste0(drug, "_10"), 
                             ylim = Ylim)
  dose_50_Track <- DataTrack(bw_list$dose_50, type = "histogram", fill = "#FB8861FF", 
                             col.histogram = "#FB8861FF", name = paste0(drug, "_50"), 
                             ylim = Ylim)
  dose_100_Track <- DataTrack(bw_list$dose_100, type = "histogram", fill = "#FEC287FF", 
                              col.histogram = "#FEC287FF", name = paste0(drug, "_100"), 
                              ylim = Ylim)
  
  # prepare extra feature tracks to display
  bm <- useMart(host = "grch37.ensembl.org", 
                biomart = "ENSEMBL_MART_ENSEMBL", 
                dataset = "hsapiens_gene_ensembl")
  
  bmt <- BiomartGeneRegionTrack(genome = "hg19", chromosome = chrom,
                                start = start_window, end = end_window, 
                                filter = list(with_ox_refseq_mrna = TRUE),
                                stacking = "dense")
  
  ideoTrack <- IdeogramTrack(genome="hg19", chromosome= chrom)
  axisTrack <- GenomeAxisTrack()
  
  allpeaks = getPeakSet(archrPrj)
  # must restrict GRanges object to specific chromosome
  peaks = allpeaks[allpeaks@seqnames == chrom,]
  pTrack <- AnnotationTrack(peaks, genome= "hg19", name = "peaks", 
                            col = "darkblue", fill = "darkblue")
  
  #plot all tracks together
  pdf(paste0("browser_tracks/gviz_", drug, "_dose_", locus, ".pdf"), width = 3, height = 6)
  print(
    plotTracks(c(ideoTrack, 
                 axisTrack, 
                 dose_0_Track,
                 dose_0.1_Track, 
                 dose_0.5_Track, 
                 dose_1_Track, 
                 dose_5_Track, 
                 dose_10_Track, 
                 dose_50_Track, 
                 dose_100_Track, 
                 pTrack,
                 bmt),
               sizes = c(0.5, 1, 1,1,1,1,1,1,1,1, 0.5, 0.5), 
               from = start_window,
               to = end_window, 
               chromosome = chrom, 
               type = "histogram")
  )
  dev.off()
}

# plot tracks using above function:

# SAHA tracks
SAHA_bw_list = list(
  dose_0 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/SAHA_0-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_0.1 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/SAHA_1-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_0.5 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/SAHA_5-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_1 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/SAHA_10-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_5 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/SAHA_50-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_10 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/SAHA_100-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_50 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/SAHA_500-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_100 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/SAHA_1000-TileSize-100-normMethod-ReadsInTSS-ArchR.bw")
)

# Dex tracks
Dex_bw_list = list(
  dose_0 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Dex_0-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_0.1 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Dex_0.5-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_0.5 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Dex_2.5-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_1 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Dex_5-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_5 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Dex_25-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_10 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Dex_50-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_50 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Dex_250-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_100 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Dex_500-TileSize-100-normMethod-ReadsInTSS-ArchR.bw")
)

# BMS tracks
BMS_bw_list = list(
  dose_0 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/BMS_0-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_0.1 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/BMS_1-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_0.5 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/BMS_5-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_1 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/BMS_10-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_5 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/BMS_50-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_10 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/BMS_100-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_50 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/BMS_500-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_100 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/BMS_1000-TileSize-100-normMethod-ReadsInTSS-ArchR.bw")
)

# Nutlin tracks
Nutlin_bw_list = list(
  dose_0 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Nutlin_0-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_0.1 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Nutlin_0.25-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_0.5 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Nutlin_1.25-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_1 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Nutlin_2.5-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_5 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Nutlin_12.5-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_10 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Nutlin_25-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_50 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Nutlin_125-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"),
  dose_100 = paste0(out_dir,"sc2_filtered/GroupBigWigs/bigWigs_TSSreads_norm/Nutlin_250-TileSize-100-normMethod-ReadsInTSS-ArchR.bw")
)


# SAHA track (COPS7A)
plot_pb_browsertracks(archrPrj = prj, drug = "SAHA",
                      bw_list = SAHA_bw_list, chrom = "chr3",
                      start_window = 171520000, end_window = 171560000, 
                      locus = "PLD1", Ylim = c(0, 20))

# Dex track (COPS7A)
plot_pb_browsertracks(archrPrj = prj, drug = "Dex",
                      bw_list = Dex_bw_list, chrom = "chr1",
                      start_window = 42753353, end_window = 42803854, 
                      locus = "FOXJ3", Ylim = c(0, 20))

# BMS track (COPS7A)
plot_pb_browsertracks(archrPrj = prj, drug = "BMS",
                      bw_list = BMS_bw_list, chrom = "chr17",
                      start_window = 78995342, end_window = 79045843, 
                      locus = "BAIAP2", Ylim = c(0, 20))

# Nutlin track (COPS7A)
plot_pb_browsertracks(archrPrj = prj, drug = "Nutlin",
                      bw_list = Nutlin_bw_list, chrom = "chr2",
                      start_window = 30102949, end_window = 30200000, 
                      locus = "ALK", Ylim = c(0, 20))







