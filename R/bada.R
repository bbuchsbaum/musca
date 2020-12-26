#collapse_obs <- function(dat) {
#  do.call(rbind, dat %>% pull(.obs) %>% purrr::map( ~ .x()))
#}

between_class_scatter <- function(X, Y, mu) {
  p <- ncol(X)
  Y <- droplevels(Y)
  levs <- levels(Y)
  
  gmeans <- group_means(Y,X)
  gmeans <- sweep(gmeans, 2, mu, "-")
  
  n <- tabulate(Y)
  
  res <- lapply(seq_along(levs), function(i) {
    n[i] * tcrossprod(gmeans[i,], gmeans[i,])
  })
  
  Reduce("+", res)
  
}


pooled_scatter <- function(X, Y) {
  ina <- as.integer(droplevels(Y))
  s <- crossprod(X)
  ni <- sqrt(tabulate(ina))
  mi <- rowsum(X, ina)/ni
  k <- length(ni)
  denom <- dim(X)[1] - k
  for (i in 1:k) s <- s - tcrossprod(mi[i, ])
  s
}

within_class_scatter <- function(X, Y) {
  pooled_scatter(X,Y)
}

## add a version that uses Sw as within subject covariance matrix

#' Barycentric Discriminant Analysis
#' 
#' A component technique that maximizes the between group variance over a set of variables. 
#' 
#' @import chk
#' @param formula the fixed effects formula
#' @param subject the subject formula
#' @param data the `multidesign` containing the X and Y data.
#'
#' @param ncomp number of components to estimate
#' @param preproc pre-processing function, defaults to \code{center}
#' @param A the column constraints
#' @param M the row constraints
#' @param ... arguments to pass through
#' 
#' @details 
#' 
#' The \code{subject} argument can be used to model multilevel structure.
#' 
#' @references
#' Abdi, H., Williams, L. J., & Bera, M. (2017). Barycentric discriminant analysis. \emph{Encyclopedia of Social Network Analysis and Mining}, 1-20.
#' 
#' 
#' @export
#' @examples 
#' 
#' X <- matrix(rnorm(50*100), 50, 100)
#' Y <- tibble(condition=rep(letters[1:5], 10), subject=rep(1:5, each=10))
#' 
#' des <- multivarious::multidesign(X,Y)
#' bada(condition, subject=subject, data=des)
bada <- function(y, subject, data, ncomp=2, preproc=center(), resdim=20, rescomp=2, ...) {
  y <- rlang::enquo(y)
  subject <- rlang::enquo(subject)
  
  labels <- factor(data$design %>% select(!!y) %>% pull(!!y))
  label_set <- levels(labels)
  
  subjects <- factor(data$design %>% select(!!subject) %>% pull(!!subject))
  subject_set <- levels(subjects)
  
  condition_means <- data %>% multivarious:::summarize_by.multidesign(sfun=colMeans, rlang::quos(y))
  
  ## condition means
  Xw <- condition_means$x
  
  ## grand mean
  Xg <- colMeans(Xw)
  
  ## data split by subject
  sdat <- split(data, subject)
  
  
  #global_preproc <- fresh(preproc)(center=Xg)
  
  ## pre-processors
  proclist <- lapply(seq_along(sdat), function(sd) {
    multivarious:::fresh(preproc) %>% prep()
  })
  
  names(proclist) <- as.character(subject_set)
  
  ## subject-split and preprocessed data
  strata <- seq_along(sdat) %>% purrr::map(function(i) {
    p <- proclist[[i]]
    Xi <- sdat[[i]]$x
    Xout <- init_transform(p, Xi)
    multidesign(Xout, sdat[[i]]$design)
  })
  
  block_indices <- list()
  ind <- 1
  for (i in 1:length(strata)) {
    block_indices[[i]] <- seq(ind, ind+ncol(strata[[i]]$x)-1)
    ind <- ind + ncol(strata[[i]]$x)
  }
  
  names(block_indices) <- as.character(subject_set)
  
  #block_indices <- lapply()
  
  ## group designs
  Dc <- strata %>% purrr::map(function(s) {
    multivarious:::summarize_by.multidesign(s, !!y)
  })
  
  ## group barycenters
  Xc <-Reduce("+", lapply(Dc, "[[", "x"))/length(Dc)
  ## dangerous, requires consistent ordering. Better to extract from design
  row.names(Xc) <- label_set
  ncomp <- min(ncomp, nrow(Xc))

  ## group pca
  pca_group <- pca(Xc, ncomp=ncomp, preproc=pass())
  
  residual_strata <- strata %>% purrr::map(function(s) {
    levs <- s$design %>% pull(!!y)
    s$x <- s$x - Xc[levs,,drop=FALSE]
    s
  })
  

  Xresid <- do.call(rbind, residual_strata %>% purrr::map( ~ .x$x))
  
  pca_resid <- pca(Xresid, ncomp=resdim, method="irlba")
  Xpca_resid <- scores(pca_resid)
  
  Sw <- within_class_scatter(Xpca_resid, interaction(subjects, labels))
  Sb <- between_class_scatter(Xpca_resid, interaction(subjects, labels), 
                              colMeans(Xpca_resid))
  
  eigout <- eigen(solve(Sw, Sb))
  #scores <- Xpca_resid %*% eigout$vectors
  resid_v <- pca_resid$v %*% eigout$vectors
  v <- cbind(pca_group$v, resid_v)
  vq <- qr(v)
  v <- qr.Q(vq)
  
  ## compute projections, one subject at a time.
  s <- do.call(rbind, strata %>% purrr::map(function(s) {
    s$x %*% v
  }))
  
  #s <- Xall %*% v

  proc <- multivarious:::concat_pre_processors(proclist, block_indices)
  
  multivarious:::discriminant_projector(v=v, s=s, apply(s,2,sd), 
                                        preproc = proc,
                                        proclist = proclist,
                                        labels=labels, 
                                        subjects=subjects,
                                        barycenters=Xc,
                                        block_indices=block_indices,
                                        classes="bada")
 
}


#' @export
reprocess.bada <- function(x, new_data, colind=NULL, subject=NULL) {
  if (is.null(colind) && is.null(subject)) {
    chk::chk_equal(ncol(new_data), shape(x)[1])
    ## pre-process every which way...
    Reduce("+", lapply(seq_along(x$block_indices), function(i) {
      apply_transform(x$preproc, new_data, colind=x$block_indices[[i]])
    }))/length(x$block_indices)
    
  } else if (!is.null(subject)) {
    ## pre-process along one block
    chk::chk_character(subject)
    chk::chk_equal(ncol(new_data), shape(x)[1])
    sind <- x$block_indices[[subject]]
    if (!is.null(colind)) {
      ## relative subset using colind
      cind <- sind[colind]
    }
    apply_transform(x$preproc, new_data, colind=cind)
  } else {
    ## colind not null. pre-process every which way using colind per block
    Reduce("+", lapply(seq_along(x$block_indices), function(i) {
      apply_transform(x$preproc, new_data, colind=colind)
    }))/length(x$block_indices)
  }
  
}

#' project 
#' 
#' @param subject the `character` id for the subject to project on. 
#' This is required for subject-specific pre-processing. 
#' If missing, the group average pre-processing parameters are used.
#' 
#' @export
#' @rdname project
project.bada <- function(x, new_data,subject) {
  if (missing(subject)) {
    NextMethod(x,new_data)
  } else {
    #Xp <- multivarious::apply_transform(preproc, new_data)
    reprocess(x, new_data, subject=subject) %*% coef(x)
  }
}

