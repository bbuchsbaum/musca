
#' @noRd
createSimFun <- function(S) {
  # Extract the labels from the similarity matrix
  label_names <- rownames(S)
  
  # The function to return
  function(labels) {
    # Find indices of input labels in the label_names
    indices <- match(labels, label_names)
    
    # Create the similarity matrix by indexing S
    sim_matrix <- S[indices, indices]
    
    return(sim_matrix)
  }
}



#' @noRd
compute_R <- function(X) {
  svd_X <- svd(t(X))
  V1 <- svd_X$v[, svd_X$d > 1]
  S1_inv_sq <- Matrix::Diagonal(x=1 / (svd_X$d[svd_X$d > 1]^2))
  R <- V1 %*% (Matrix::Diagonal(ncol(V1)) - S1_inv_sq) %*% t(V1)
  return(R)
}


#' @noRd
compute_rank_matrices <- function(strata) {
  Sl <- purrr::map(strata, function(x) {
    R <- compute_R(x$x)
  })
  
  Sl
  
}

#' @examples
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=1, run=rep(1:5, 2)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=2, run=rep(1:5, 2)))
#' d3 <- multidesign(matrix(rnorm(10*20), 10, 20), data.frame(y=1:10, subject=3, run=rep(1:5, 2)))
#' S <- matrix(runif(10*10), 10, 10)
#' S <- abs(cor(S))
#' row.names(S) <- colnames(S) <- 1:10
#' simfun <- createSimFun(S)
#' hd <- hyperdesign(list(d1,d2,d3))
lowrank_align.hyperdesign <- function(data, y, 
                                     preproc=center(), 
                                     ncomp=2,
                                     simfun,
                                     u=.5,
                                     lambda=.0001) {
  
  y <- rlang::enquo(y)
  
  labels <- factor(unlist(purrr::map(data, function(x) x$design %>% select(!!y) %>% pull(!!y))))
  label_set <- levels(labels)
  
  ninstances <- length(labels)
  nsets <- length(data)
  
  pdata <- multivarious::init_transform(data, preproc) 
  proclist <- attr(pdata, "preproc")
  
  names(proclist) <- names(pdata)
  
  block_indices <- block_indices(pdata)
  proc <- multivarious::concat_pre_processors(proclist, block_indices)
  names(block_indices) <- names(pdata)
  
  #browser()
  Rs <- compute_rank_matrices(pdata)
  R_block <- Matrix::bdiag(Rs)
  C_block <- simfun(labels)
  C_block[C_block < 0] <- 0
  C_block <- Matrix(C_block, sparse=TRUE)
  C_block <- (C_block + t(C_block))/2
  
  
  diag(C_block) <- 0
  #C_block <- neighborweights::normalize_adjacency(C_block)
  
  #M <- (diag(N) - R) %*% (diag(N) - R)
  
  R_i <- (Matrix::Diagonal(nrow(R_block)) - R_block)
  M <- crossprod(R_i, R_i)
  
  Ds <- Matrix::Diagonal(x=Matrix::rowSums(C_block))
  Ls <- Ds - C_block
  s_m <- PRIMME::eigs_sym(M, NEig=1,  which="LA",method='PRIMME_DEFAULT_MIN_MATVECS')
  s_l <- PRIMME::eigs_sym(Ls, NEig=1, which="LA",method='PRIMME_DEFAULT_MIN_MATVECS')
  e_ratio <- s_l$values[1]/ s_m$values[1] 
  M <- M * e_ratio
  
  Z <- (1-u) * M + (2 * u * Ls)

  # Solve the eigenvalue problem
  eigen_result <- PRIMME::eigs_sym(Z, NEig=ncomp, which="SA",method='PRIMME_DEFAULT_MIN_MATVECS')
  #eigen_result2 <- RSpectra::eigs_sym(Z, 2, which="SM")
  
  v <- eigen_result$vectors
  Xp <- Matrix::bdiag(lapply(pdata, function(x) x$x))
  
  if (is.null(lambda)) {
    cvfit <- cv.glmnet(Xp, v, family = "mgaussian", alpha=0, intercept=FALSE)
    lambda <- cvfit$lambda.min
  } 
  rfit <- glmnet(Xp, v, family = "mgaussian", alpha=0, lambda=lambda, intercept=FALSE)
  cfs <- coef(rfit)
  cfs <- do.call(cbind, cfs)[-1,,drop=FALSE]
  
  multivarious::multiblock_biprojector(
    v=cfs,
    s=v,
    sdev=apply(v,2,sd),
    preproc=proc,
    block_indices=block_indices,
    labels=labels,
    u=u,
    classes="lowrank_align"
  )
  
}