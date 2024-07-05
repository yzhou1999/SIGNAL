###################################################################
##### SIGNAL core functions
###################################################################
# Written by Y. Zhou on May 2, 2024.

#' Group-centralized principal component analysis
#'
#' @param X Input data matrix.
#' @param meta Input metadata (data.frame).
#' @param g_factor Group variable (s).
#' @param b_factor Batch variable (s).
#' @param lambda Tuning parameter.
#' @param npcs How many dimensions to reduce.
#' @param excluded.cells Metadata about which cells should not be used.
#' @param output.all Whether to return all calculated values.
#' @param do.scale Whether to perform global scaling.
#' @param do.cosine Whether to cosine normalize the integration result.
#' @param block.sizeGB Memory settings for block computing.
#' @param block.size Block size.
#'
#' @import Matrix
#' @import bigstatsr
#' @import RSpectra
#' @import RcppEigen
#' @importFrom dplyr `%>%`
#' @importFrom Rcpp evalCpp
#' @useDynLib SIGNAL
#'
#' @export
Run.gcPCA <- function(X, meta, g_factor, b_factor, lambda = 50, npcs = 30, excluded.cells = NULL,
                      output.all = FALSE, do.scale = F, do.cosine = F,
                      block.sizeGB = 2, block.size = 4000) {
  if (do.scale) X <- scale_data(X, do.center = F)
  if (!is.matrix(X)) X <- as.matrix(X)
  options(bigstatsr.block.sizeGB = block.sizeGB)

  X <- bigstatsr::big_copy(X, block.size = block.size)$save()

  message("Run gcPCA!")
  N = length(g_factor)
  if (N == 1) {
    G <- group_onehot(meta, g_factor)
    BG <- group_onehot(meta, c(g_factor, b_factor))
    if (!is.null(excluded.cells)) {
      G[excluded.cells, ] <- 0
      BG[excluded.cells, ] <- 0
    }

    XXt <- bigstatsr::big_tcrossprodSelf(X, big_scale(center = FALSE, scale = FALSE), block.size = block.size)[]
    total_mmt <- group_mmt(X, matrix(1, nrow = ncol(X), ncol = 1) )
    g_mmt <- group_mmt(X, G)
    bg_mmt <- group_mmt(X, BG)
    V <- RSpectra::eigs_sym(XXt - total_mmt - lambda*(bg_mmt - g_mmt), k = npcs, which = "LA")[["vectors"]] %>% as.matrix()
    gcpca_res <- t(bigstatsr::big_cprodMat(X, V))
    if (do.cosine) gcpca_res <- cosineNorm(gcpca_res)
    message("gcPCA done!")
    if (output.all) gcpca_res = list('gcpca_res' = gcpca_res,
                                     'XXt' = XXt,
                                     'total_mmt' = total_mmt,
                                     'g_mmt' = g_mmt,
                                     'bg_mmt' = bg_mmt,
                                     'V' = V)
  } else if (N == 2) {
    G_list = lapply(g_factor, function(g) group_onehot(meta, g))
    BG_list = lapply(g_factor, function(g) group_onehot(meta, c(g, b_factor)))
    cross_G = group_onehot(meta, g_factor)
    cross_BG = group_onehot(meta, c(g_factor, b_factor))
    if (!is.null(excluded.cells)) {
      G_list[[1]][excluded.cells, ] <- 0
      G_list[[2]][excluded.cells, ] <- 0
      BG_list[[1]][excluded.cells, ] <- 0
      BG_list[[2]][excluded.cells, ] <- 0
      cross_G[excluded.cells, ] <- 0
      cross_BG[excluded.cells, ] <- 0
    }

    XXt <- bigstatsr::big_tcrossprodSelf(X, big_scale(center = FALSE, scale = FALSE), block.size = block.size)[]
    total_mmt <- group_mmt(X, matrix(1, nrow = ncol(X), ncol = 1) )
    g_mmt <- lapply(G_list, function(G) group_mmt(X, G))
    bg_mmt <- lapply(BG_list, function(BG) group_mmt(X, BG))
    cross_g_mmt <- group_mmt(X, cross_G)
    cross_bg_mmt <- group_mmt(X, cross_BG)
    V <- RSpectra::eigs_sym(XXt - total_mmt -
                              lambda*(Reduce('+', bg_mmt) - Reduce('+', g_mmt) -
                                        cross_bg_mmt + cross_g_mmt), k = npcs, which = "LA")[["vectors"]] %>% as.matrix()
    gcpca_res <- t(bigstatsr::big_cprodMat(X, V))
    if (do.cosine) gcpca_res <- cosineNorm(gcpca_res)
    message("gcPCA done!")
    if (output.all) gcpca_res = list('gcpca_res' = gcpca_res,
                                     'XXt' = XXt,
                                     'total_mmt' = total_mmt,
                                     'g_mmt' = g_mmt,
                                     'bg_mmt' = bg_mmt,
                                     'V' = V)
  }
  return(gcpca_res)
}

#' Predict cell labels using group technial effects.
#'
#' @param ref_data Reference data matrix.
#' @param ref_group Reference labels.
#' @param query_data Query data matrix.
#' @param do.cosine Whether to cosine normalize the data matrices.
#' @param do.scale Whether to perform scaling.
#'
#' @import Matrix
#' @importFrom dplyr `%>%`
#'
#' @export
Run.PGTV.Single <- function(ref_data, ref_group, query_data, do.cosine = T, do.scale = F) {
  if (do.cosine == TRUE) {
    ref_data <- cosineNorm(ref_data)
    query_data <- cosineNorm(query_data) %>% as("dgCMatrix")
  }
  if (do.scale == TRUE) {
    ref_data <- scale_data(ref_data, do.center = F)
    query_data <- scale_data(query_data, do.center = F)
  }
  groups_K <- unique(ref_group)
  K <- length(groups_K)
  N_k <- as.matrix(table(ref_group))[groups_K,]
  x_mean_k <- sapply(1:K, function(k) {
    if (length(which(ref_group == groups_K[k])) > 1) {
      rowMeans(ref_data[, which(ref_group == groups_K[k])])
    } else {
      ref_data[, which(ref_group == groups_K[k])]
    }
  })
  colnames(x_mean_k) <- groups_K
  x_mean_square_kK <- sapply(1:K, function(k) N_k[k]*sum(x_mean_k[, k]^2))
  x_mean_kK <- sapply(1:K, function(k) N_k[k]*x_mean_k[, k])
  pred_query_group <- groups_K[pred_k(x_mean_square_kK, query_data, x_mean_kK, as.numeric(N_k))]
  return(pred_query_group)
}

#' Predict cell labels for a data matrix using its given cell labels.
#'
#' @param data Data matrix.
#' @param group Given labels.
#' @param do.cosine Whether to cosine normalize the data matrices.
#' @param do.scale Whether to perform scaling.
#' @param max.iters The maximum number of iterations.
#' @param thresh Whether to perform scaling.
#' @param do.scale Stopping threshold for classification agreement metrics.
#'
#' @import Matrix
#' @import mclust
#' @importFrom dplyr `%>%`
#'
#' @export
Run.PGTV.Single.Cycle <- function(data, group, do.cosine = T, do.scale = F, max.iters = 20, thresh = 1e-5) {
  if (is.null(colnames(data)) ) colnames(data) <- paste0("cell", 1:ncol(data))
  names(group) <- colnames(data)
  ref_data <- data
  if (do.cosine == TRUE) ref_data <- cosineNorm(ref_data) %>% as("dgCMatrix")
  if (do.scale == TRUE) ref_data <- scale_data(ref_data, do.center = F)
  ref_group <- group
  query_data <- ref_data
  query_group <- ref_group

  iters <- 1
  ari_value0 <- 0
  ari_diff <- 1
  aris <- c()
  ref_index <- colnames(ref_data)
  while (iters < max.iters & ari_diff > thresh) {
    message(paste0("Perform prediction round ", iters))
    groups_K <- unique(ref_group)
    K <- length(groups_K)
    N_k <- as.matrix(table(ref_group))[groups_K,]
    x_mean_k <- sapply(1:K, function(k) {
      if (length(which(ref_group == groups_K[k])) > 1) {
        rowMeans(ref_data[, which(ref_group == groups_K[k])])
      } else {
        ref_data[, which(ref_group == groups_K[k])]
      }
    })
    colnames(x_mean_k) <- groups_K
    x_mean_square_kK <- sapply(1:K, function(k) N_k[k]*sum(x_mean_k[, k]^2))
    x_mean_kK <- sapply(1:K, function(k) N_k[k]*x_mean_k[, k])
    pred_query_group <- groups_K[pred_k(x_mean_square_kK, query_data, x_mean_kK, as.numeric(N_k))]
    names(pred_query_group) <- colnames(data)
    ari_value <- mclust::adjustedRandIndex(pred_query_group[ref_index], ref_group)
    ari_diff <- abs(ari_value - ari_value0) / ari_value
    aris <- c(aris, ari_value)
    message(paste0("ARI value = ", ari_value))
    message(paste0("ARI diff = ", ari_diff))
    ref_index <- ref_index[which(ref_group == pred_query_group[ref_index])]
    ref_data <- query_data[, which(colnames(query_data) %in% ref_index)]
    ref_group <- query_group[which(colnames(query_data) %in% ref_index)]
    iters <- iters + 1
    ari_value0 <- ari_value
  }
  print(aris)
  return(pred_query_group)
}

#' Predict cell labels of query data using cell labels of reference data.
#'
#' @param ref_data Reference data matrix.
#' @param ref_group Reference labels.
#' @param query_data Query data matrix.
#' @param do.cosine Whether to cosine normalize the data matrices.
#' @param do.scale Whether to perform scaling.
#' @param conf.score Whether to calculate confidence scores.
#' @param proj.dims Projection dimension of query data to reference data.
#' @param eta The proportion of cells that are considered reliably projected.
#' @param smooth.k Using kNNs to smooth the scores.
#'
#' @import Matrix
#' @import BiocNeighbors
#' @import irlba
#' @importFrom sparseMatrixStats rowSds
#' @importFrom dplyr `%>%`
#'
#' @export
Run.LabelTransfer.Single <- function(ref_data, ref_group, query_data, do.cosine = T, do.scale = F,
                                     conf.score = T, proj.dims = 30, eta = 0.25, smooth.k = 5) {
  if (is.list(ref_data)) {
    ref_meta <- data.frame("Batch" = Reduce('c', lapply(1:length(ref_data), function(i) rep(paste0("Batch_", i), ncol(ref_data[[i]])))),
                           "Group" = ref_group)
    rownames(ref_meta) <- Reduce('c', lapply(ref_data, colnames))
    ref_data <- do.call(cbind, ref_data)
  }
  if (do.cosine == TRUE) {
    ref_data <- cosineNorm(ref_data) %>% as("dgCMatrix")
    query_data <- cosineNorm(query_data) %>% as("dgCMatrix")
  }
  if (conf.score == TRUE) query_sds <- sparseMatrixStats::rowSds(query_data)
  if (do.scale == TRUE) {
    ref_sds <- sparseMatrixStats::rowSds(ref_data)
    ref_data <- scale_data(ref_data, do.center = F, row.sds = ref_sds)
    if (conf.score == TRUE) {
      query_data <- scale_data(query_data, do.center = F, row.sds = query_sds)
    } else {
      query_data <- scale_data(query_data, do.center = F)
    }
  }
  groups_K <- unique(ref_group)
  K <- length(groups_K)
  N_k <- as.matrix(table(ref_group))[groups_K,]
  x_mean_k <- sapply(1:K, function(k) {
    if (length(which(ref_group == groups_K[k])) > 1) {
      rowMeans(ref_data[, which(ref_group == groups_K[k])])
    } else {
      ref_data[, which(ref_group == groups_K[k])]
    }
  })
  colnames(x_mean_k) <- groups_K
  x_mean_square_kK <- sapply(1:K, function(k) N_k[k]*sum(x_mean_k[, k]^2))
  x_mean_kK <- sapply(1:K, function(k) N_k[k]*x_mean_k[, k])
  pred_query_labels <- groups_K[pred_k(x_mean_square_kK, query_data, x_mean_kK, as.numeric(N_k))]
  if (!conf.score) return(pred_query_labels)

  query_mean <- rowMeans(query_data)
  scaled_query <- scale_data(query_data, row.means = query_mean, do.scale = F)
  query_projV <- irlba(t(scaled_query), nv=proj.dims)[["v"]]
  if (do.scale == TRUE) {
    scaled_ref_means <- scale_data(x_mean_k[, unique(pred_query_labels)], row.means = query_mean, row.sds = query_sds/ref_sds)
  } else {
    scaled_ref_means <- scale_data(x_mean_k[, unique(pred_query_labels)], row.means = query_mean, row.sds = query_sds)
  }
  query_refdata <- crossprod(scaled_ref_means, query_projV) %>% t() %>% cosineNorm()
  query_querydata <- crossprod(scaled_query, query_projV) %>% t() %>% cosineNorm()
  sim_query_ref <- crossprod(query_refdata, query_querydata)
  rownames(sim_query_ref) <- unique(pred_query_labels)
  confidence_scores <- sapply(1:length(pred_query_labels), function(i) sim_query_ref[pred_query_labels[i], i])
  # Adjust the score
  query_knns <- BiocNeighbors::findKNN(t(query_querydata), k = smooth.k)$index
  knn_scores <- matrix(confidence_scores[query_knns], ncol = smooth.k)
  confidence_scores <- sapply(1:length(pred_query_labels), function(i) mean(knn_scores[i, ]))

  max.score <- quantile(confidence_scores, 1-eta)
  min.score <- min(quantile(confidence_scores, 0.01), 0)
  confidence_scores <- (confidence_scores - min.score) / (max.score - min.score)
  confidence_scores[confidence_scores > 1] <- 1
  confidence_scores[confidence_scores < 0] <- 0

  return(transfer_res = data.frame(row.names = colnames(query_data),
                                   "Prediction" = pred_query_labels,
                                   "Confidence" = confidence_scores))
}


#' Compute the sum of the outer product of the scaled group means
#'
#' @param X Input FBM object.
#' @param G Group one-hot matrix.
#' @param output.all Whether to output all results.
#'
#' @import Matrix
#' @import bigstatsr
group_mmt <- function(X, G, output.all = FALSE) {
  K <- ncol(G)
  Ns <- colSums(G)
  tG <- Matrix::t(as(G, "dgCMatrix"))
  K_indexs <- tG@i + 1
  indexs <- lapply(1:K, function(k) which(G[, k] > 0))

  mmt <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  if (output.all) {
    group_means <- matrix(0, nrow = nrow(X), ncol = K)
    for (k in 1:K) {
      if (length(indexs[[k]]) == 1) {
        group_means[, k] <- X[, indexs[[k]] ]
        mmt <- mmt + group_means[, k] %o% group_means[, k] * Ns[k]
      } else if (length(indexs[[k]]) > 1) {
        group_means[, k] <- Matrix::rowMeans(X[, indexs[[k]] ])
        mmt <- mmt + group_means[, k] %o% group_means[, k] * Ns[k]
      }
    }

    return(list(group_means = group_means,
                Ns = Ns,
                mmt = mmt))
  } else {
    for (k in 1:K) {
      if (length(indexs[[k]]) == 1) {
        x <- X[, indexs[[k]] ]*sqrt(Ns[k])
        mmt <- mmt + x %o% x
      } else if (length(indexs[[k]]) > 1) {
        x <- Matrix::rowMeans(X[, indexs[[k]] ])*sqrt(Ns[k])
        mmt <- mmt + x %o% x
      }
    }

    return(mmt)
  }
}

#' Compute covariance matrix
#'
#' @param X Input FBM object.
#' @param G Group ont-hot matrix.
#' @param block.size Block size.
#'
#' @import Matrix
#' @import bigstatsr
compute_cov <- function(X, G, block.size) {
  X0 <- big_copy(X)
  group_center(X0, G)
  bigstatsr::big_tcrossprodSelf(X0, big_scale(center = FALSE, scale = FALSE), block.size = block.size)[]
}

#' Compute one-hot matrix for given data frame and variable (s)
#'
#' @param x Input data frame.
#' @param ivar Variable (s) for one-hot computation.
#'
#' @importFrom stats model.matrix
group_onehot <- function(x, ivar) {
  if (length(unique(x[, ivar])) == 1) {
    matrix(1, nrow = length(x[, ivar]), ncol = 1)
  } else {
    x <- data.frame(ivar = x[, ivar])
    x <- Reduce(paste0, x)
    model.matrix(~ 0 + x)
  }
}

#' Scale data matrix
#'
#' @param data.x Input data matrix.
#' @param do.center Whether center the row values. (default TRUE)
#' @param do.scale Whether scale the row values. (default TRUE)
#' @param row.means The provided row means to center. (default NULL)
#' @param row.sds The provided row standard deviations to scale. (default NULL)
#'
#' @import Matrix
scale_data <- function(data.x,
                       do.center = T,
                       do.scale = T,
                       row.means = NULL,
                       row.sds = NULL) {
  if (do.center) {
    if (is.null(row.means)) {
      data_mean <- Matrix::rowMeans(data.x)
    } else {
      data_mean <- row.means
    }
    data.x <- data.x - sapply(1:ncol(data.x), function(i) data_mean)
  }
  if (do.scale) {
    if (is.null(row.sds)) {
      data_stddev <- matrixStats::rowSds(as.matrix(data.x))
    } else {
      data_stddev <- row.sds
    }
    index <- which(data_stddev > 0)
    data.x[index, ] <- data.x[index, ] / sapply(1:ncol(data.x), function(i) data_stddev[index])
  }

  data.x
}

#' Cosine normalize data matrix
#' @param x Input data matrix.
#' 
#' @import Matrix
cosineNorm <- function(x) {
  l2 <- sqrt(Matrix::colSums(x^2))
  l2 <- pmax(1e-8, l2)
  mat <- scale(x, center = F, scale = l2)
  mat
}
