#' @keywords internal
compute_sim_mat <- function(blocks, FUN, ...) {
  pairs <- combn(length(blocks),2)
  M <- matrix(0, length(blocks), length(blocks))
  
  for (i in 1:ncol(pairs)) {
    p1 <- pairs[1,i]
    p2 <- pairs[2,i]
    sim <- FUN(blocks[[p1]], blocks[[p2]], ...)
    M[p1,p2] <- sim
    M[p2,p1] <- sim
  }
  
  M
}

#' @keywords internal
normalization_factors <- function(blocks, type=c("MFA", "RV", "RV2", "None", "Frob")) {
  type <- match.arg(type)
  message("normalization type:", type)
  alpha <- if (type == "MFA") {
    unlist(lapply(blocks, function(X) 1/(multivarious::svd_wrapper(X, ncomp=1, method="svds")$sdev[1]^2)))
  } else if (type == "RV" && length(blocks) > 2) {
    smat <- compute_sim_mat(blocks, function(x1,x2) MatrixCorrelation::RV2(x1,x2))
    diag(smat) <- 1
    abs(svd_wrapper(smat, ncomp=1, method="svds")$u[,1])
  } else if (type == "RV2" && length(blocks) > 2) {
    smat <- compute_sim_mat(blocks, function(x1,x2) MatrixCorrelation::RV(x1,x2))
    diag(smat) <- 1
    abs(svd_wrapper(smat, ncomp=1, method="svds")$u[,1])
  } else if (type == "Frob") {
    unlist(lapply(as.list(blocks), function(X) sum(X^2)))
  } else {
    rep(1, length(blocks))
  }
}




#' multiple factor analysis
#' 
#' principal component analysis for multiple blocks of data collected over th same instances
#' 
#' 
#' @param X a \code{block_matrix} object
#' @param ncomp the number of components to estimate
#' @param preproc a preprocessing pipeline, default is `center()`
#' @param normalization the normalization method: MFA, RV, RV-MFA, or None (see details).
#' @param A when `normalization` is `custom`, a weight (sparse or dense) matrix for the columns
#' @param M when `normalization` is `custom`, a weight (sparse or dense) matrix for the rows
#' @export
#' 
#' @references
#' Abdi, H., Williams, L. J., & Valentin, D. (2013). Multiple factor analysis: principal component analysis for multitable and multiblock data sets. 
#' \emph{Wiley Interdisciplinary reviews: computational statistics}, 5(2), 149-179.
#' @export
#' 
#' @examples 
#' 
#' X <- replicate(5, { matrix(rnorm(10*10), 10, 10) }, simplify=FALSE)
#' res <- mfa(X, ncomp=3, normalization="MFA")
#' p <- project_block(res, X[[1]], 1)
#' stopifnot(ncol(scores(res)) == 3)
#' 
#' labs <- letters[1:10]
#' cfier <- classifier(res, new_data=do.call(cbind, X), labels=labs)
#' pred <- predict(cfier, X[1:2,])
#' cfier2 <- classifier(res, new_data=X[[2]], labels=labs, colind=res$block_indices[[2]])
#' pred2 <- predict(cfier2, X[1:2,res$block_indices[[2]]])
mfa.list <- function(data, preproc=center(), ncomp=2,
                normalization=c("MFA", "RV", "None", "Frob", "custom"), 
                M=NULL, A=NULL, ...) {
  
  
  chk::chk_true(length(data) > 1)
  for (i in 1:length(data)) {
    chk::chkor(chk::chk_matrix(data[[i]]), chk::chk_s4_class(data[[i]], "Matrix"))
  }
  
  nrs <- sapply(data, nrow)
  chk::chk_true(all(nrs == nrs[1]))
  nr <- nrs[1]
  
  normalization <- match.arg(normalization)
  
  if (normalization == "custom") {
    chk::chkor(chk::chk_not_null(A), chk::chk_not_null(M))
  }
  
  if (is.null(names(data))) {
    names(data) <- paste0("B", 1:length(data))
  }
  
  ## set up pre-processors
  proclist <- lapply(seq_along(data), function(dat) {
    multivarious:::fresh(preproc) %>% prep()
  })
  
  names(proclist) <- names(data)
  
  
  ## future possibility, for large datasets.
  ## if p >> n for each block, then could project to n components, then concatenate.
  ## first must check if A exists and has any weights that span across blocks.
  
  
  ## pre-process blocks
  strata <- seq_along(data) %>% purrr::map(function(i) {
    p <- proclist[[i]]
    Xi <- data[[i]]
    Xout <- init_transform(p, Xi)
    Xout
  })
  
  ## calculate block normalization factors
  if (normalization != "custom") {
    alpha <- normalization_factors(strata, type=normalization)
    A <- rep(alpha, sapply(strata, ncol))
  } else {
    alpha <- rep(1, length(strata))
  }
  
  ## compute block indicees
  block_indices <- list()
  ind <- 1
  for (i in 1:length(strata)) {
    block_indices[[i]] <- seq(ind, ind+ncol(strata[[i]])-1)
    ind <- ind + ncol(strata[[i]])
  }
  
  proc <- multivarious::concat_pre_processors(proclist, block_indices)

  ## fit genpca
  Xp <- do.call(cbind, strata)
  fit <- genpca::genpca(Xp, 
                preproc=pass(),
                A=A, 
                M=M,
                ncomp=ncomp,
                ...)
  
  fit[["block_indices"]] <- block_indices
  fit[["alpha"]] <- alpha
  fit[["normalization"]] <- normalization
  fit[["names"]] <- names(data)
  
  ## this is awkward...
  ## instead, we need a "delegation" mechanism, where a multiblock projector simply wraps a projector
  ## here, we rely on the fact that we use "pass()" pre-processing for inner genpca fit
  fit[["preproc"]] <- proc

  class(fit) <- c("mfa", "multiblock_biprojector", "multiblock_projector", class(fit))
  fit
 
  
}
