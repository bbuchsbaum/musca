#collapse_obs <- function(dat) {
#  do.call(rbind, dat %>% pull(.obs) %>% purrr::map( ~ .x()))
#}





#' @import chk
#' @param resdim pca dimensionality for residual analysis (only relevant if `rescomp` > 0)
#' @param rescomp number of final residual components (default = 0, no residual aanalysis)
#' @param ... arguments to pass through
#' @rdname bada
#' @export
bada.multidesign <- function(data, y, subject, preproc=center(), ncomp=2,resdim=20, rescomp=0, ...) {
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
  
  ## requires consistent ordering. Better to extract from design
  row.names(Xc) <- label_set
  ncomp <- min(ncomp, nrow(Xc))

  ## group pca
  pca_group <- pca(Xc, ncomp=ncomp, preproc=pass())
  
  
  ## residual analysis
  if (rescomp > 0) {
    chk::chk_true(resdim > 0)
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
  } else {
    resdim <- 0
    rescomp <- 0
    v <- pca_group$v
  }
  
  ## compute projections, one subject at a time.
  s <- do.call(rbind, strata %>% purrr::map(function(s) {
    s$x %*% v
  }))
  
  proc <- multivarious:::concat_pre_processors(proclist, block_indices)
  
  multivarious:::discriminant_projector(v=v, s=s, 
                                        sdev=apply(s,2,sd), 
                                        preproc = proc,
                                        proclist = proclist,
                                        labels=labels, 
                                        resdim=resdim,
                                        rescomp=rescomp,
                                        subjects=subjects,
                                        barycenters=Xc,
                                        block_indices=block_indices,
                                        classes="bada")
 
}


#' @export
#' 
reprocess.bada <- function(x, new_data, colind=NULL, block=NULL) {
  if (is.null(colind) && is.null(block)) {
    ## how to pre-process wheen you don't know the subject?
    ## we pre-process every way and average.
    chk::chk_equal(ncol(new_data), shape(x)[1])
    
    ## pre-process every which way...
    Reduce("+", lapply(seq_along(x$block_indices), function(i) {
      apply_transform(x$preproc, new_data, colind=x$block_indices[[i]])
    }))/length(x$block_indices)
    
  } else if (!is.null(block)) {
    ## pre-process along one block
    chk::chk_character(block)
    chk::chk_equal(ncol(new_data), shape(x)[1])
    sind <- x$block_indices[[block]]
    if (!is.null(colind)) {
      ## relative subset using colind
      sind <- sind[colind]
    }
    apply_transform(x$preproc, new_data, colind=sind)
  } else {
    ## colind not null. pre-process every which way using colind per block
    Reduce("+", lapply(seq_along(x$block_indices), function(i) {
      apply_transform(x$preproc, new_data, colind=colind)
    }))/length(x$block_indices)
  }
  
}

#' new sample projection
#' 
#' Project one or more samples of onto a subspace.
#' 
#' @inheritParams multivarious::project
#' 
#' @param block the `character` id for the subject/block to project on. 
#' This is required for block-specific pre-processing. 
#' If missing, the group average pre-processing parameters are used.
#' 
#' @export
project.bada <- function(x, new_data, block) {
  if (missing(block)) {
    NextMethod(x,new_data)
  } else {
    #Xp <- multivarious::apply_transform(preproc, new_data)
    reprocess(x, new_data, block=block) %*% coef(x)
  }
}

