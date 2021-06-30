options(stringsAsFactors = FALSE)
library(Matrix)
library(Seurat)
library(ggplot2)
library(irlba)
library(patchwork)
library(plyr)
library(dplyr)
library(stringr)
library(SnapATAC)
library(GenomicRanges)
library(cisTopic)

########################################
# Utility functions
########################################
# Loads 10x dataset into a sparse matrix object
# Args:
#   matrix_fn (string): name of mtx or mtx.gz file coontaining MM formatted matrix
#   peaks_fn (string): name of peaks BED file (matches rows of matrix)
#   barcodes_fn (string): name of file containing cell barcodes (matches columns of matrix)
# Returns:
#   sparse matrix: matrix represented by files above with row and column names set (chr_start_stop format used for rownames)
load_tenx_atac = function(matrix_fn, peaks_fn, barcodes_fn) {
  atac_matrix = readMM(matrix_fn)
  colnames(atac_matrix) = read.delim(barcodes_fn, header=FALSE)$V1
  peaks = read.delim(peaks_fn, header=FALSE)
  peaks = paste(peaks$V1, peaks$V2, peaks$V3, sep = '_')
  rownames(atac_matrix) = peaks
  return(atac_matrix)
}

# Allows filtering of sites measured as non-zero in less than a given number of cells
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   cells (int): filter sites if they have less than this number of cells above zero
# Returns:
#   sparse matrix: filtered sparse matrix
filter_features = function(bmat, cells=10) {
  bmat = bmat[Matrix::rowSums(bmat) >= cells, ]
  return(bmat)
}

# Allows filtering of cells with below a given number of non-zero features
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   features_above_zero (int): filter cells if they have less than this number of features above zero
# Returns:
#   sparse matrix: filtered sparse matrix
filter_cells = function(bmat, features_above_zero=100) {
  bmat = bmat[, Matrix::colSums(bmat > 0) >= features_above_zero]
  return(bmat)
}

# Takes sparse matrix object and downsamples to a given fraction of entries remaining.
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
#   fraction_remaining (float): float (0, 1) that indicates fraction of non-zero entries to retain
#   cells_per_site_min (int): min cells a site must be measured in to retain the site in matrix
#   sites_per_cell_min (int): min sites a cell must have non-zero entries in to retain the cell in matrix
# Returns:
#   sparse matrix: downsampled sparse matrix
downsample_matrix = function(bmat, fraction_remaining=0.5, cells_per_site_min=1, sites_per_cell_min=1) {
  set.seed(2019)
  non_zero_entries = which(bmat@x > 0)
  indices_to_zero = sample(non_zero_entries, size=ceiling(length(non_zero_entries) * (1 - fraction_remaining)))
  bmat@x[indices_to_zero] = 0
  
  # Make sure to get rid of stuff that has gone to ~ 0 after downsampling
  bmat = filter_features(bmat, cells=cells_per_site_min)
  bmat = filter_cells(bmat, features_above_zero=sites_per_cell_min)
  return(bmat)
}

########################################
# Functions for LSI
########################################
# Helper function to do fast version of row scaling of sparse TF matrix by IDF vector.
# Exploits the way that data is stored within sparseMatrix object. Seems to be much more memory efficient than tf * idf and faster than DelayedArray.
# Args:
#   tf (sparse matrix): term frequency matrix
#   idf (vector): IDF vector
# Returns:
#   sparse matrix: TF-IDF matrix
safe_tfidf_multiply = function(tf, idf) {
   tf = t(tf)
   tf@x <- tf@x * rep.int(idf, diff(tf@p))
   tf = t(tf)
   return(tf)
}

# Perform TF-IDF on binary matrix
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
#   frequencies (bool): divide bmat by colSums (if FALSE simply use bmat for TF matrix)
#   log_scale_tf (bool): log scale TF matrix if TRUE
#   scale_factor (float): multiply terms in TF matrix by scale_factor prior to log1p. Equivalent to adding small pseudocount but doesn't cast to dense matrix at any point.
# Returns:
#   sparse matrix: TF-IDF matrix
tfidf = function(bmat, frequencies=TRUE, log_scale_tf=TRUE, scale_factor=100000) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  
  # TF-IDF
  tf_idf_counts = safe_tfidf_multiply(tf, idf)
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  return(tf_idf_counts)
}

# Perform current version of TF-IDF used by 10x on binary matrix
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
# Returns:
#   sparse matrix: TF-IDF matrix
tenx_tfidf = function(bmat) {
  idf = log(ncol(bmat) + 1) - log(1 + Matrix::rowSums(bmat))
  tf_idf_counts = safe_tfidf_multiply(bmat, idf)
  
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  tf_idf_counts = as(tf_idf_counts, "sparseMatrix")
  return(tf_idf_counts)
}

# Perform fast PCA (irlba) on matrix, retaining observation names
# Args:
#   mat (sparse matrix): matrix to use for PCA (no further scaling or centering done)
#   dims (int): number of PCs to calculate
# Returns:
#   sparse matrix: TF-IDF matrix
do_pca = function(mat, dims=50) {
  pca.results = irlba(t(mat), nv=dims)
  final_result = pca.results$u %*% diag(pca.results$d)
  rownames(final_result) = colnames(mat)
  colnames(final_result) = paste0('PC_', 1:dims)
  return(final_result)
}

########################################
# Helper functions for dim reduction
########################################
# Wrapper for performing further dim reduction (tSNE/UMAP) and clustering given PCA space via Seurat.
# Args:
#   atac_matrix (sparse matrix): matrix to store in Seurat object (not used in computations)
#   cell_embeddings (matrix): typically PCA coordinates of cells but could be any set of reduced dim coordinates
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   metadata (dataframe): dataframe of metadata (rowonames are cell names) to add to Seurat object
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
# Returns:
#   Seurat object: seurat object
run_dim_reduction = function(atac_matrix, cell_embeddings, dims, metadata=NULL, reduction='pca.l2') {
  if (is.null(metadata)) {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix)
  } else {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix, meta.data = metadata)
  }

  seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=cell_embeddings, key='PC_', assay='RNA')
  seurat_obj = seurat_obj %>%
                Seurat::L2Dim(reduction='pca') %>%
                Seurat::RunUMAP(reduction = reduction, dims = dims) %>%
#                Seurat::RunTSNE(reduction = reduction, dims = dims) %>%
                Seurat::FindNeighbors(reduction=reduction, nn.eps=0.25, dims=dims)
  return(seurat_obj)
}

# Helper function for plotting Spearman correlations of a given metadata column with all dimensions in a reduced space.
# Args:
#   seurat_obj (seurat object): Seurat object to use
#   reduction (string): name of reduction to use
#   column (string): name of column in metadata to use
# Returns:
#   ggplot object: plot object
plot_pc_correlation = function(seurat_obj, reduction, column='nCount_RNA') {
  coords = Seurat::Embeddings(seurat_obj, reduction=reduction)
  column_value = seurat_obj@meta.data[, column]
  correlations = apply(coords, 2, function(x) {cor(x, column_value, method='spearman')})
  correlations_df = data.frame(correlation=correlations, PC=1:ncol(coords))
  
  plot_obj = ggplot(correlations_df, aes(PC, correlation)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 0, linetype='dashed', color='red')
    
  return(plot_obj)
}

########################################
# Wrapper functions for workflows
########################################
# Wrapper for full LSI workflow (TF-IDF and PCA + clustering + further dim reduction)
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   metadata: dataframe of metadata (rowonames are cell names) to add to Seurat object
#   log_scale_tf (bool): log scale TF matrix if TRUE
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on PCA results from TF-IDF matrix.
lsi_workflow = function(bmat, dims, metadata=NULL, log_scale_tf=TRUE, reduction='pca.l2', resolution=0.3) {
  tfidf_mat = tfidf(bmat, log_scale_tf=log_scale_tf)
  pca_mat = do_pca(tfidf_mat, dims=max(dims))
  
  seurat_obj = run_dim_reduction(bmat, pca_mat, dims, metadata) %>%
                Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)
  return(seurat_obj)
}

# Wrapper for 10x version of full LSI workflow (TF-IDF and PCA + clustering + further dim reduction). Only TF-IDF step is modified.
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   metadata: dataframe of metadata (rownames are cell names) to add to Seurat object
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on PCA results from TF-IDF matrix.
tenx_lsi_workflow = function(bmat, dims, metadata=NULL, reduction='pca.l2', resolution=0.3) {
  tfidf_mat = tenx_tfidf(bmat)
  pca_mat = do_pca(tfidf_mat, dims=max(dims))
  
  seurat_obj = run_dim_reduction(bmat, pca_mat, dims, metadata) %>%
                Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)
  return(seurat_obj)
}

# Runs cisTopic on binary matrix.
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   topic (vector of int): topic numbers to generate models for
# Returns:
#   cisTopic object: cisTopic object with models generated
cistopic_workflow = function(bmat, topic=c(40, 50)) {
  coords = str_split_fixed(rownames(bmat), '_', 3)
  new_coords = paste0(coords[, 1], ':', coords[, 2], '-', coords[, 3])
  rownames(bmat) = new_coords

  cisTopicObject = cisTopic::createcisTopicObject(bmat, project.name='mouse_atac')
  cisTopicObject = cisTopic::runModels(cisTopicObject, topic, seed=2019, nCores=1, burnin = 250, iterations = 500)
  return(cisTopicObject)
}

# Wrapper for SnapATAC workflow up until PCA step.
# Args:
#   snap_file (string): path to snap file
#   promooter.df (dataframe): dataframe with promoter definitions as shown in SnapATAC tutorials (dataframe not file name; see examples in markdown)
#   blacklist.df (dataframe): dataframe with blacklist region definitions as shown in SnapATAC tutorials (dataframe not file name; see examples in markdown)
#   fragment_number_threshold (int): threshold for number of unique fragments per cell
#   promoter_ratio_range (c(float, float)): vector with lower and upper bound of acceptable fraction of reads in promoter regions for cell filtering as used in SnapATAC tutorial.
#   window_z_range (c(float, float)): vector with lower and upper bound of acceptable window z-scores for non-zero entries for site filtering as used in SnapATAC tutorial.
#   sample_name (string): sample_name provided to SnapATAC
#   pc.num (int): total PCs to compute
# Returns:
#   snap object: SnapATAC object
# Notes:
#   This uses single core because multithreaded implementation interferes with Knitr. In running this code, any do.par=FALSE and num.cores=1 could be changed as needed.
snapatac_workflow = function(snap_file, promoter.df=NULL, blacklist.df=NULL, fragment_number_threshold=500, promoter_ratio_range=c(0.2, 0.8), window_z_range=c(-1.5, 1.5), sample_name='default', pc.num=50) {
  x.sp = createSnap(
    file=snap_file,
    sample=sample_name,
    do.par=FALSE,
    num.cores=1
  )
  
  plotBarcode(
    obj=x.sp, 
    pdf.file.name=NULL, 
    pdf.width=7, 
    pdf.height=7, 
    col="grey",
    border="grey",
    breaks=50
  )
  
  x.sp = filterCells(
    obj=x.sp, 
    subset.names=c("fragment.num", "UMI"),
    low.thresholds=c(fragment_number_threshold, fragment_number_threshold),
    high.thresholds=c(Inf, Inf)
  )
  
  x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=1)

  # Optionally filter cells based on ratio of reads in promoters
  if (!is.null(promoter.df)) {
    promoter.gr = GRanges(promoter.df[,1], IRanges(promoter.df[,2], promoter.df[,3]))
    ov = findOverlaps(x.sp@feature, promoter.gr)
    idy = queryHits(ov)
    promoter_ratio = SnapATAC::rowSums(x.sp[,idy, mat="bmat"], mat="bmat") / SnapATAC::rowSums(x.sp, mat="bmat")
    plot(log(SnapATAC::rowSums(x.sp, mat="bmat") + 1,10), promoter_ratio, cex=0.5, col="grey", xlab="log(count)", ylab="FIP Ratio", ylim=c(0,1 ))
    idx = which(promoter_ratio > promoter_ratio_range[1] & promoter_ratio < promoter_ratio_range[2])
    x.sp = x.sp[idx,]
  }
  
  x.sp = makeBinary(x.sp, mat="bmat");
  
  # Filter out non-standard contigs if present
  idy2 = grep("chrM|random", x.sp@feature)
  
  if (!is.null(blacklist.df)) {
  black_list.gr = GRanges(
    blacklist.df[,1], 
    IRanges(blacklist.df[,2], blacklist.df[,3])
  )
    idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr))
    
  } else {
    # No blacklist provided, so just ignore
    idy1 = c()
  }
  
  idy = unique(c(idy1, idy2))
  x.sp = x.sp[,-idy, mat="bmat"]
  
  # Filter based on frequency
  x.sp = filterBins(
    x.sp,
    low.threshold=window_z_range[1],
    high.threshold=window_z_range[2],
    mat="bmat"
  )
  
  plotBinCoverage(
    x.sp,
    pdf.file.name=NULL,
    col="grey",
    border="grey",
    breaks=10,
    xlim=c(-6,6)
  )
  
  x.sp = runJaccard(
    obj = x.sp,
    tmp.folder=tempdir(),
    mat = "bmat",
    max.var=2000,
    ncell.chunk=1000,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )

  x.sp = runNormJaccard(
    obj = x.sp,
    tmp.folder=tempdir(),
    ncell.chunk=1000,
    method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )

  x.sp = runDimReduct(
    x.sp,
    pc.num=pc.num,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=10
  )

  rownames(x.sp@bmat) = x.sp@barcode
  colnames(x.sp@bmat) = as.character(1:ncol(x.sp@bmat))

  return(x.sp)
}

# Reperform Jaccard + PCA on SnapATAC object (used for redoing after modifying matrix)
# Args:
#   snao_obj (snap object): snap object
#   pc.num (int): total PCs to compute
# Returns:
#   snap object: SnapATAC object
# Notes:
#   This uses single core because multithreaded implementation interferes with Knitr. In running this code, any do.par=FALSE and num.cores=1 could be changed as needed.
snapatac_rerun_jaccard = function(snap_obj, pc.num=50) {
  
  snap_obj = runJaccard(
    obj = snap_obj,
    tmp.folder=tempdir(),
    mat = "bmat",
    max.var=2000,
    ncell.chunk=1000,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )

  snap_obj = runNormJaccard(
    obj = snap_obj,
    tmp.folder=tempdir(),
    ncell.chunk=1000,
    method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )

  snap_obj = runDimReduct(
    snap_obj,
    pc.num=pc.num,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=10
  )
}

# Wrapper for cisTopic workflow (choose from models that have been run + clustering + further dim reduction)
# Args:
#   cistopicObject (cistopic object): cisTopic object with runModels having already been run
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   method (string): method argument to modelMatSelection function in cisTopic
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on topic matrix from cisTopic
seurat_workflow_on_cistopic = function(cistopicObject, method='Z-score', reduction='pca', resolution=0.3) {
  cistopicObject = cisTopic::selectModel(cistopicObject)

  cistopicObject.reduced_space = t(cisTopic::modelMatSelection(cistopicObject, target='cell', method=method))
  colnames(cistopicObject.reduced_space) = paste0('PC_', 1:ncol(cistopicObject.reduced_space))
  dimensions = ncol(cistopicObject.reduced_space)
  
  cistopicObject.seurat = run_dim_reduction(cistopicObject@binary.count.matrix, cistopicObject.reduced_space, dims=1:dimensions, reduction='pca')
  
  cistopicObject.seurat = cistopicObject.seurat %>% 
    Seurat::FindNeighbors(reduction=reduction, nn.eps=0.25, dims=1:dimensions) %>%
    Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)
  return(cistopicObject.seurat)
}

# Wrapper for running downstream Seurat workflow (clustering + further dim reduction) on PCA from Jaccard matrix generated by SnapATAC
# Args:
#   snao_obj (snap object): snap object with runDimReduct already run
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on PCA matrix from SnapATAC (note PCA is weighted by variance explained)
seurat_workflow_on_jaccard_pca = function(snap_obj, dims, reduction='pca', resolution=0.3) {
  pca_results.jaccard = snap_obj@smat@dmat %*% diag(snap_obj@smat@sdev)
  colnames(pca_results.jaccard) = paste0('PC_', 1:ncol(pca_results.jaccard))
  rownames(pca_results.jaccard) = snap_obj@barcode

  seurat_obj.jaccard = run_dim_reduction(t(snap_obj@bmat), pca_results.jaccard, dims, reduction=reduction)
  seurat_obj.jaccard = seurat_obj.jaccard %>%
    Seurat::FindNeighbors(nn.eps=0.25, dims=dims, reduction=reduction) %>%
    Seurat::FindClusters(n.start=20, resolution=resolution, dims=dims, reduction=reduction)
}

# Helper function to plot embeddings corresponding to same set of cells with one another's cluster assignments for comparison.
# Args:
#   seurat_obj1 (snap object): snap object with runDimReduct already run
#   seurat_obj2 (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   reduction (string): reduction to use for plot (can be tsne or umap)
#   description1 (string): title for seurat_obj1 (used in plot)
#   description2 (string): title for seurat_obj1 (used in plot)
#   cluster_column1 (string): column from metadata of seurat_obj1 to use for coloring plot
#   cluster_column2 (string): column from metadata of seurat_obj2 to use for coloring plot
# Returns:
#   ggolot object: ggplot object where each embedding is plotted colored by its own clusters and then again with the opposite object's clusters assigned for comparison. Four total panels shown.
plot_clustering_comparison = function(seurat_obj1, seurat_obj2, reduction, description1='', description2 = '', cluster_column1='RNA_snn_res.0.3', cluster_column2='RNA_snn_res.0.3') {
  # Clusters as called on each dataset
  seurat_obj1 = SetIdent(seurat_obj1, value=cluster_column1)
  seurat_obj2 = SetIdent(seurat_obj2, value=cluster_column2)
  
  p1 = DimPlot(seurat_obj1, reduction = 'tsne', pt.size=0.25) +
            ggtitle(description1)
  
  p2 = DimPlot(seurat_obj2, reduction = 'tsne', pt.size=0.25) +
            ggtitle(description2)
  
  # Now swap the labels to verify they are finding the same groups
  seurat_obj1@meta.data$cluster_seurat_obj2 = seurat_obj2@meta.data[, cluster_column2]
  seurat_obj2@meta.data$cluster_seurat_obj1 = seurat_obj1@meta.data[, cluster_column1]
  
  seurat_obj1 = SetIdent(seurat_obj1, value='cluster_seurat_obj2')
  seurat_obj2 = SetIdent(seurat_obj2, value='cluster_seurat_obj1')
  
  p3 = DimPlot(seurat_obj1, reduction = reduction, pt.size=0.25) +
            ggtitle(paste0(description1, ' (', description2, ' clusters)'))
  
  p4 = DimPlot(seurat_obj2, reduction = reduction, pt.size=0.25) +
            ggtitle(paste0(description2, ' (', description1, ' clusters)'))
  
  (p1 + p3) / (p2 + p4)
}