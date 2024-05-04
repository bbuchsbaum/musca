

#' @noRd
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

#' @noRd
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

#' @noRd
within_class_scatter <- function(X, Y) {
  pooled_scatter(X,Y)
}