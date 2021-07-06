
##Xs <- lapply(1:5, function(i) matrix(rnorm(10*10), 10, 10))
##Xs <- lapply(Xs, cor)

double_center <- function(x) {
  1/2 * t(scale(t(scale(x, scale=FALSE)), scale=FALSE))
}

mat_inner_prod <- function(S1, S2) {
  sum(as.vector(S1) * as.vector(S2))
}

norm_crossprod <- function(S) {
  S/sqrt(sum(S^2))
}

compute_prodmat <- function(S) {
  C <- matrix(0, length(S), length(S))
  
  for (i in 1:length(S)) {
    for (j in 1:length(S)) {
      if (i == j) {
        C[i,j] <- 1
      } else if (i > j) {
        next
      } else {
        rv <- mat_inner_prod(S[[i]], S[[j]])
        C[i,j] <- rv
        C[j,i] <- rv
      }
    }
  }
  
  C
}


#' @export
covstatis.list <- function(data, ncomp=2, normalize=TRUE, dcenter=TRUE, labels=NULL) {
  
  nr <- sapply(data, nrow)
  nc <- sapply(data, ncol)
  
  assertthat::assert_that(all(nr == nr[1]))
  assertthat::assert_that(all(nc == nc[1]))
  assertthat::assert_that(all(nc[1] == nr[1]))
  
  block_labels <- names(data)
  
  if (is.null(block_labels)) {
    block_labels <- paste0("B_", 1:length(data))
  }
  
  if (is.null(labels)) {
    labels <- row.names(data[[1]])
    if (is.null(labels)) {
      labels <- paste0("Obs_", 1:nr[1])
    }
  } else {
    assertthat::assert_that(length(labels) == nr[1])
  }
  
 
  
  S <- data
  
  if (dcenter) {
    S <- lapply(S, double_center)
  }
  
  if (normalize) {
    S <- lapply(S, norm_crossprod)
  }
  
  C <- compute_prodmat(S)

  alpha <- abs(eigen(C)$vectors[,1])
  alpha <- alpha/(sum(alpha))
  
  Sall <- Reduce("+", lapply(1:length(alpha), function(i) S[[i]]*alpha[i]))
  fit <- eigen(Sall)
  keep <- which(1:length(fit$values) <= ncomp & fit$values > 1e-8)
  
  scores <- fit$vectors[,keep,drop=FALSE] %*% diag(sqrt(fit$values[keep]))
  projmat <- fit$vectors[,keep,drop=FALSE] %*% diag(1/sqrt(fit$values[keep]))
  
  ret <- multivarious::projector(fit$vectors[,keep,drop=FALSE], classes="covstatis", 
                                 s=scores,
                                 projmat=projmat,
                                 normalize=normalize, dcenter=dcenter, alpha=alpha, C=C,
                                 block_labels=block_labels, labels=labels)
  ret
  
}

#' @export
project_cov.covstatis <- function(x, new_data) {
  new_data <- as.matrix(new_data)
  assertthat::assert_that(nrow(new_data) == nrow(x$v), msg="`new_data` must be symmetric")
  assertthat::assert_that(isSymmetric(new_data), msg="`new_data` must be symmetric")
  
  if (x$double_center) {
    new_data <- double_center(new_data)
  }
  
  if (x$normalize) {
    new_data <- norm_crossprod(new_data)
  } 
  
  new_data %*% x$projmat
}

