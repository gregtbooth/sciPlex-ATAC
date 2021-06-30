library(Matrix)
library(irlba)
library(Seurat)
library(dplyr)
library(plyr)
library(stringr)
library(GenomicRanges)
options(stringsAsFactors = FALSE)

.get_aux_files = function(mtx_file) {
  base_name = str_replace(mtx_file, '[.]gz$', '')
  base_name = str_replace(base_name, '[.]mtx$', '')
  
  features_file = paste0(base_name, '.rows.txt')
  cells_file = paste0(base_name, '.columns.txt')
  
  return(list(features=features_file, cells=cells_file))
}

load_mtx_file = function(mtx_file) {
  mat = readMM(mtx_file)
  dim_files = .get_aux_files(mtx_file)
  
  if (! file.exists(dim_files$features)) {
    stop(paste0(dim_files$features, ' file not found when loading ', mtx_file))
  }
  
  if (! file.exists(dim_files$cells)) {
    stop(paste0(dim_files$cells, ' file not found when loading ', mtx_file))
  }
  
  rownames(mat) = read.delim(dim_files$features, header=FALSE)$V1
  colnames(mat) = read.delim(dim_files$cells, header=FALSE)$V1
  return(as(mat, "dgCMatrix"))
}

########################################
# Utility functions
########################################
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

filter_features_z = function(bmat, lower_bound=-1.5, upper_bound=1.5, downsample=NULL) {
  feature_totals = log10(Matrix::rowSums(bmat) + 1)
  avg = mean(feature_totals)
  stdev = sd(feature_totals)

  feature_z = (feature_totals-avg)/stdev

  # Deal with corner case for missing or entirely non-zero windows (latter probably not a thing)
  feature_z[feature_totals == 0] = -Inf
  feature_z[feature_totals == ncol(bmat)] = Inf

  feature_totals = feature_totals[feature_z > lower_bound & feature_z < upper_bound]
  
  if (!is.null(downsample)) {
    probabilities = pnorm(feature_totals, mean=avg, sd=stdev)
    sampled_features = sample(names(feature_totals), prob=probabilities, size=downsample, replace=FALSE)
    feature_totals = feature_totals[sampled_features]
  }
  return(bmat[names(feature_totals),])
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

# Slightly better binarization from SnapATAC
binarize_matrix = function(mat, outlier.filter=1e-3) {
    if (max(mat@x) == 1) {
      return(mat)
    }

		# identify the cutoff using outlier.filter
		count_cutoff = max(1, quantile(mat@x, 1 - outlier.filter))
		mat@x[mat@x > count_cutoff] = 0
		mat@x[mat@x > 0] = 1
    return(mat)
}

safe_tfidf_multiply = function(tf, idf) {
  tf = t(tf)
  tf@x <- tf@x * rep.int(idf, diff(tf@p))
  tf = t(tf)
  return(tf)
}

safe_column_scale = function(bmat, scale_vector) {
  bmat@x <- bmat@x / rep.int(scale_vector, diff(bmat@p))
  return(bmat)
}

safe_column_multiply = function(bmat, scale_vector) {
  bmat@x <- bmat@x * rep.int(scale_vector, diff(bmat@p))
  return(bmat)
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

sparse_prcomp_irlba <- function(x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, ...) {
  a <- names(as.list(match.call()))
  ans <- list(scale=scale.)
  
  args <- list(A=x, nv=n, scale.=scale., center=center)
  if (!missing(...)) args <- c(args, list(...))
  s <- do.call(irlba, args=args)
  
  ans$sdev <- s$d / sqrt(max(1, nrow(x) - 1))
  ans$rotation <- s$v
  colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  ans$center <- args$center
  if (retx)
  {
    ans <- c(ans, list(d=s$d, x = s$u %*% diag(s$d)))
    colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), sep="")
  }
  class(ans) <- c("irlba_prcomp", "prcomp")
  ans
}

# Perform fast PCA (irlba) on matrix, retaining observation names
# Args:
#   mat (sparse matrix): matrix to use for PCA (no further scaling or centering done)
#   dims (int): number of PCs to calculate
# Returns:
#   sparse matrix: TF-IDF matrix
do_pca = function(mat, dims=50) {
  pca.results = sparse_prcomp_irlba(t(mat), n=dims, center=FALSE, scale.=FALSE)
  final_result = pca.results$x
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
    #Seurat::RunTSNE(reduction = reduction, dims = dims) %>%
    Seurat::FindNeighbors(reduction=reduction, nn.eps=0.5, dims=dims)
  return(seurat_obj)
}

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

atac_scrublet = function(bmat, k=NULL, fraction_sim_doublets=2, estimated_doublet_rate=0.1, dims=2:50) {
  # Determine KNN parameters
  if (is.null(k)) {
    k = round(0.5 * sqrt(ncol(bmat)))
  }
  kadj = round(k * (1+fraction_sim_doublets))
  
  # Perform TFIDF on original dataset
  message('[scrublet atac] Performing LSI-logTF on dataset...')
  
  # TF-IDF on original dataset
  bmat_colsums = Matrix::colSums(bmat)
  tf_original = safe_column_scale(bmat, bmat_colsums)
  tf_original@x = log1p(tf_original@x * 100000)
  idf_original = log(1 + ncol(bmat) / Matrix::rowSums(bmat))

  bmat.tfidf = safe_tfidf_multiply(tf_original, idf_original)
  rm(tf_original)

  bmat.pca = sparse_prcomp_irlba(t(bmat.tfidf), n=max(dims), center=FALSE, scale.=FALSE)
  rm(bmat.tfidf)
  
  # Make simulated doublets
  message('[scrublet atac] Simulating doublets...')
  set.seed(2019)
  doublet_half1 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
  doublet_half2 = sample(colnames(bmat), floor(ncol(bmat) * fraction_sim_doublets), replace=TRUE)
  doublets.bmat = bmat[, doublet_half1] + bmat[, doublet_half2]
  doublets.bmat@x[doublets.bmat@x > 1] = 1
  colnames(doublets.bmat) = paste0('DOUBLET_', 1:ncol(doublets.bmat))
  
  # Perform TF-IDF on doublets using IDF from original dataset
  doublet_colsums = bmat_colsums[doublet_half1] + bmat_colsums[doublet_half2] ## approximate (because of binarization after sum) to save time, but could recalculate
  tf_doublets = safe_column_scale(doublets.bmat, doublet_colsums)
  rm(doublets.bmat)
  tf_doublets@x = log1p(tf_doublets@x * 100000)
  doublets.tfidf = safe_tfidf_multiply(tf_doublets, idf_original)
  rm(tf_doublets)

  # Project doublets into PCA space and weight by variance explained
  message('[scrublet atac] Projecting doublets into PCA space...')
  doublets.pca = t(doublets.tfidf) %*% bmat.pca$rotation
  rm(doublets.tfidf)
  doublets.pca.weighted = doublets.pca
  
  # Jam it all into a Seurat object for some of the downstream stepst
  message('[scrublet atac] Making Seurat object...')
  rownames(bmat.pca$x) = colnames(bmat)
  combined.metadata = rbind(data.frame(cell=colnames(bmat), doublet=FALSE), data.frame(cell=rownames(doublets.pca), doublet=TRUE))
  rownames(combined.metadata) = combined.metadata$cell
  
  combined.pca = rbind(bmat.pca$x, doublets.pca.weighted)
  
  rownames(combined.pca) = rownames(combined.metadata)
  combined.seurat = CreateSeuratObject(counts=t(combined.pca), meta.data = combined.metadata)
  combined.seurat[['pca']] = Seurat::CreateDimReducObject(embeddings=as.matrix(combined.pca), key='PC_', assay='RNA')
  
  message('[scrublet atac] Finding KNN...')
  combined.seurat = combined.seurat %>%
    Seurat::L2Dim(reduction='pca') %>%
    Seurat::FindNeighbors(reduction='pca.l2', dims=dims, k=kadj, compute.SNN = FALSE) # nn.eps = 0.25    
  
  # From KNN, calculate doublet likelihood as defined in Scrublet paper
  message('[scrublet atac] Calculating doublet neighbors...')
  doublet_mask = ifelse(combined.seurat@meta.data$doublet, 1, 0)
  doublet_neighbors = Matrix::rowSums(safe_column_multiply(combined.seurat@graphs$RNA_nn, doublet_mask))

  message('[scrublet atac] Finalizing doublet likelihoods...')
  doublet_score = doublet_neighbors / kadj
  q = (doublet_neighbors + 1)/ (kadj + 2)
  doublet_likelihood = q * estimated_doublet_rate / fraction_sim_doublets / (1 - estimated_doublet_rate - q * (1 - estimated_doublet_rate - estimated_doublet_rate / fraction_sim_doublets))
  
  # Return Seurat object with doublet likelihood as an extra column
  result = data.frame(cell=rownames(combined.seurat@meta.data), 'doublet_score'=doublet_score, 'doublet_likelihood'=doublet_likelihood, 'simulated_doublet'=combined.seurat@meta.data$doublet)
  rownames(result) = result$cell
  return(result)
}

# Takes a set of features chr_start_stop and a dataframe of chr, start, end and returns only features that don't overlap regions in the dataframe
filter_regions = function(features, blacklist_df) {
  column_names = c('chrom', 'start', 'end')

  if (! all(column_names %in% colnames(blacklist_df))) {
    stop('chrom, start, and end must be columns in blacklist df.')
  }

  features_df = as.data.frame(str_split_fixed(features, '_', n=3))
  colnames(features_df) = column_names
  features_df$start = as.numeric(features_df$start)
  features_df$end = as.numeric(features_df$end)

  features_df.gr = GRanges(
    features_df$chrom,
    IRanges(features_df$start, features_df$end)
  )

  black_list.gr = GRanges(
    blacklist_df$chrom,
    IRanges(blacklist_df$start, blacklist_df$end)
  )

  matching_hits = queryHits(findOverlaps(features_df.gr, black_list.gr))
  return(features[-matching_hits])
}

downsample_clusters = function(mat, clusters, max_cells_per_cluster=Inf) {
  downsampled_cells = c()
  for (id in unique(clusters)) {
    cells_in_cluster = colnames(mat)[clusters == id]

    if (length(x = cells_in_cluster) > max_cells_per_cluster) {
      cells_in_cluster = sample(x = cells_in_cluster, size = max_cells_per_cluster, replace=FALSE)
    }

    downsampled_cells = c(downsampled_cells, cells_in_cluster)
  }
  return(mat[, downsampled_cells])
}
