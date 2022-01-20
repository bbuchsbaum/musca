
#  ## use normalized
#  A <- Wsn - dweight*Wrn
#  decomp <- PRIMME::eigs_sym(u*Wn + (1-u)*A, NEig=ncomp)
#  Y <- decomp$vectors[,1:ncomp]
#}


rescale <- function(z) {
  t(apply(z,1, function(x) x/sqrt(sum(x^2))))
}

#' @keywords internal
coskern <- function() {
  rval <- function(x, y = NULL) {
    if (!is(x, "vector")) 
      stop("x must be a vector")
    if (!is(y, "vector") && !is.null(y)) 
      stop("y must be a vector")
    if (is(x, "vector") && is.null(y)) {
      coop::cosine(x)
    }
    if (is(x, "vector") && is(y, "vector")) {
      if (!length(x) == length(y)) 
        stop("number of dimension must be the same on both data points")
      coop::cosine(x, y)
    }
  }
  
  return(new("kernel", .Data = rval, kpar = list()))
}

#' @keywords internal
stratified_subsample <- function(labels, nperlabel) {
  nlabs <- length(unique(labels))
  sp <- split(seq_along(labels), labels)
  unlist(lapply(sp, function(idx) sample(idx, nperlabel)))
}

#' @keywords internal
normalize_laplacian <- function(L, D) {
  dvals <- diag(D)
  dvals[dvals == 0] <- 1e-4
  Dinvsq <-  Matrix::Diagonal(x=1/sqrt(dvals))
  Dinvsq %*% L %*% Dinvsq
}

#' @keywords internal
normalize_adjacency <- function(A, D) {
  Dinvsq <-  Matrix::Diagonal(x=1/sqrt(diag(D)))
  Dinvsq %*% A %*% Dinvsq
}


# medoids <- function(X, n=20, sfrac=1) {
#   d <- kmed::distNumeric(X,X)
#   res <- kmed::fastkmed(d,n, iterate=50)
#   res$medoid
# }

# 
# kmns <- function(X, k) {
#   res <- kmeans(X, round(k))
#   tmp <- FNN::get.knnx(X, res$centers, k=1)
#   tmp$nn.index[,1]
# }

#' @keywords internal
kcentroids <- function(X, k, sfrac=.5) {
  chk::chk_true(sfrac > 0 && sfrac <= 1)
  #res <- kmeans(X, round(k))
  #res$centers
  if (sfrac < 1) {
    n <- round(sfrac * ncol(X))
    sidx <- sort(sample(1:ncol(X), n))
    res <- cluster::pam(X[,sidx,drop=FALSE], k, metric="manhattan")
  } else{
    res <- cluster::pam(X, k, metric="manhattan")
  }
  
  sort(res$id.med)
 
}

#' @keywords internal
class_medoids <- function(X, L) {
  L <- as.factor(L)
  sl <- split(1:length(L), L)
  ids <- sapply(sl, function(sidx) {
    cp <- cluster::pam(X[sidx,], k=1, metric="manhattan")
    sidx[cp$id.med]
  })
  
  ids
}




compute_local_similarity <- function(strata, y, knn, weight_mode, type, sigma) {
  y <- rlang::enquo(y)
  Sl <- purrr::map(strata, function(x) {
    labs <- x$design %>% select(!!y) %>% pull(!!y)
    g <- neighborweights::graph_weights(x$x,
                                        weight_mode=weight_mode,
                                        neighbor_mode="knn",
                                        k=knn,
                                        type=type,
                                        sigma=sigma)
    
    #labs <- labels[x$design$.index]
    cg <- neighborweights::class_graph(labs)
    r <- neighborweights::repulsion_graph(g, cg, method="weighted")
    list(G=adjacency(g), R=adjacency(r))
  })
  
  Sl
  
}

compute_laplacians <- function(Ws, Wr, W, Wd, normalize=FALSE) {
  Ds <- Matrix::Diagonal(x=Matrix::rowSums(Ws))
  Ls <- Ds - Ws
 
    
  Dr <- Matrix::Diagonal(x=Matrix::rowSums(Wr))
  Lr <- Dr - Wr
  
    
  D <- Matrix::Diagonal(x=Matrix::rowSums(W))
  L <- D - W
  
  Dd <- Matrix::Diagonal(x=Matrix::rowSums(Wd))
  Ld <- Dd - Wd
  
  if (normalize) {
    list(Ls=normalize_laplacian(Ls, Ds),
         Lr=normalize_laplacian(Lr, Dr),
         L=normalize_laplacian(L, D),
         Ld=normalize_laplacian(Ld, Dd))
  } else {
    list(Ls=Ls,
         Lr=Lr,
         L=L,
         Ld=Ld)
  }
  
}

compute_between_graph <- function(strata, y) {
  y <- rlang::enquo(y)
  
  dlabels <- lapply(strata, function(s) {
    labs <- s$design %>% select(!!y) %>% pull(!!y)
    meds <- sort(class_medoids(s$x, labs))
    medlabels <- rep(NA, length(labs))
    medlabels[meds] <- names(meds)
    medlabels
  })
  
  neighborweights::binary_label_matrix(unlist(dlabels), unlist(dlabels), type="d")
  
}

compute_kernels <- function(strata, kernel, sample_frac) {
  Ks <- if (sample_frac == 1) {
    purrr::map(strata, function(x) {
      k <- kernlab::kernelMatrix(kernel, rescale(x$x))
      k/nrow(k)
    })
  } else {
    sidx <- vector(length(strata), mode="list")
    purrr::map(seq_along(strata), function(i) {
      x <- rescale(strata[[i]]$x)
      n <- round(nrow(x) * sample_frac)
      ci <- kcentroids(x, n)
      sidx[[i]] <- ci
      k <- t(kernlab::kernelMatrix(kernel, x, x[ci,]))
      k/nrow(k)
    })
  }
  
  Ks
}

normalize_graphs <- function(Sl, Ws, Wd) {
  diag(Ws) <- 0
  
  
  W <- Matrix::bdiag(lapply(Sl, "[[", "G"))
  Wr <- Matrix::bdiag(lapply(Sl, "[[", "R"))
  
  Sws <- sum(Ws)
  Sw <- sum(W)
  
  Ws <- Ws/Sws * Sw
  
  Swr <- sum(Wr)
  Wr <- Wr/Swr * Sw
  
  Swd <- sum(Wd)
  Wd <- Wd/Swd * Sw
  
  list(W=W, Wr=Wr, Ws=Ws, Wd=Wd)
}
 

## features by instances (d by n)
# Xlist <- lapply(1:5, function(i) matrix(rnorm(50*20), 50, 20))
# Y <-  lapply(1:5, function(i) sample(letters[1:3],50, replace=TRUE))
#' @import Matrix
kema.multidesign <- function(data, y, 
                             subject, 
                             preproc=center(), 
                             ncomp=2, 
                             knn=5, 
                             sigma=.73, u=.5, 
                             kernel=coskern(), 
                             sample_frac=1,
                             use_laplacian=TRUE, 
                             specreg=TRUE,
                             dweight=.1) {
  
  subject <- rlang::enquo(subject)
  y <- rlang::enquo(y)
  subjects <- factor(data$design %>% select(!!subject) %>% pull(!!subject))
  subject_set <- levels(subjects)
  
  strata <- multidesign::hyperdesign(split(data, subject))
  kema(strata, !!y, preproc, ncomp, knn, sigma, u, kernel, sample_frac, use_laplacian,specreg, dweight)

 
}

#' @importFrom multivarious init_transform
kema.hyperdesign <- function(data, y, 
                             preproc=center(), 
                             ncomp=2, 
                             knn=5, 
                             sigma=.73, 
                             u=.5, 
                             kernel=coskern(), 
                             sample_frac=1,
                             use_laplacian=TRUE, 
                             specreg=TRUE,
                             dweight=.1) {
  
  chk::chk_number(ncomp)
  chk::chk_range(sample_frac, c(0,1))
  chk::chk_logical(specreg)
  chk::chk_number(dweight)
  chk::chk_number(knn)
  chk::chk_range(u, c(0,1))
  
  
  y <- rlang::enquo(y)
  
  #browser()
  labels <- factor(unlist(purrr::map(data, function(x) x$design %>% select(!!y) %>% pull(!!y))))
  label_set <- levels(labels)
  
  ninstances <- length(labels)
  nsets <- length(data)
  
  pdata <- init_transform(data, preproc) 
  proclist <- attr(pdata, "preproc")
 
  names(proclist) <- names(pdata)
  
  block_indices <- block_indices(pdata)
  
  proc <- multivarious::concat_pre_processors(proclist, block_indices)
  names(block_indices) <- names(pdata)
  
  kema_fit(pdata, proc, ncomp, knn, sigma, u, !!y, labels, kernel, sample_frac, specreg, dweight, block_indices)
  
}

kema_fit <- function(strata, proc, ncomp, knn, sigma, u, y, labels, kernel, sample_frac, specreg, dweight, block_indices) {
  chk::chk_number(ncomp)
  chk::chk_range(sample_frac, c(0,1))
  chk::chk_logical(specreg)
  chk::chk_number(dweight)
  chk::chk_range(u, c(0,1))
  
  y <- rlang::enquo(y)
  
  ## data similarity
  Sl <- compute_local_similarity(strata, !!y, knn, 
                                 weight_mode="normalized", 
                                 type="normal",  
                                 sigma=sigma)
  
  ## class pull
  Ws <- neighborweights::binary_label_matrix(labels, labels)
  
  ## class push
  Wd <- compute_between_graph(strata, !!y)
  
  ## reweight graphs
  G <- normalize_graphs(Sl,Ws, Wd)
  
  ## compute full or subsampled kernels
  Ks <- compute_kernels(strata, kernel, sample_frac)
  kernel_indices <- get_block_indices(Ks, byrow=TRUE)
  Z <- Matrix::bdiag(Ks)
  
  ## compute laplacians
  Lap <- compute_laplacians(G$Ws,G$Wr,G$W,G$Wd, specreg)
  
  kemfit <- kema_solve(strata, Z, Ks, Lap, kernel_indices, specreg, ncomp, u, dweight, sample_frac)

  multivarious::multiblock_biprojector(
    v=kemfit$coef,
    s=kemfit$scores,
    sdev=apply(kemfit$scores,2,sd),
    preproc=proc,
    alpha=kemfit$vectors,
    block_indices=block_indices,
    Ks=Ks,
    sample_frac,
    dweight=dweight,
    labels=labels,
    classes="kema"
  )
}


#' @import glmnet
kema_solve <- function(strata, Z, Ks, Lap, kernel_indices, specreg, ncomp, u, dweight, sample_frac, lambda=.0001) {
  if (specreg) {
    A <- Lap$Ls - dweight*(Lap$Lr + Lap$Ld)
    decomp <- PRIMME::eigs_sym(u*Lap$L + (1-u)*A, NEig=ncomp+1, which="SA")
    Y <- decomp$vectors[,1:ncomp]
    
    if (sample_frac < 1) {
      rfit <- glmnet(t(Z), Y, family = "mgaussian", alpha=0, lambda=lambda)
    } else {
      rfit <- glmnet(Z, Y, family = "mgaussian", alpha=0, lambda=lambda)
    }
    cfs <- coef(rfit)
    cfs <- do.call(cbind, cfs)[-1,,drop=FALSE]
    vectors <- cfs
  } else {
    
    Zl <- Z %*% (u*Lap$L + (1-u)*Lap$Ls) %*% t(Z)
    Zr <- Z %*% (Lap$Lr + Lap$Ld) %*% t(Z)
    decomp <- trace_ratio(Zr, Zl, ncomp=ncomp+1)
    vectors <- decomp$vectors[,1:ncol(decomp$vectors),drop=FALSE]
  }
  
  scores <- if (sample_frac == 1) {
    ## S <- K * E 
    Z %*% vectors
  } else {
    Matrix::t(Z) %*% vectors
  }
  
  v <- do.call(rbind, lapply(1:length(strata), function(i) {
    xi <- strata[[i]]$x
    kind <- seq(kernel_indices[i,1], kernel_indices[i,2])
    alpha_i = vectors[kind,,drop=FALSE]
    
    if (sample_frac < 1) {
      v_i = t(xi) %*% t(Ks[[i]]) %*% alpha_i
    } else {
      v_i = t(xi) %*% alpha_i
    }
    
    ## vi_a (project from feature to latent space)
    v_i
    #t(xi) %*% scores[kind,] %*% diag(x=1/decomp$values)
  }))
  
  
  list(coef=v, scores=scores, vectors=vectors)
}


  
