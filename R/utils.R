#' @keywords internal
get_block_indices <- function(Xlist, byrow=FALSE) {
  lens <- if (byrow) {
    sapply(Xlist, nrow)
  } else {
    sapply(Xlist, ncol)
  }
  
  csum <- cumsum(lens)
  csum1 <- c(0, csum[-length(csum)])
  m <- as.matrix(cbind(csum1+1, csum))
  colnames(m) <- c("start", "end")
  m
}


#' trace_ratio optimization
#' 
#' @param A the numerator matrix
#' @param B the denominator matrix
#' 
#' @keywords internal
#' 
#' 
#' @export
#' @examples
#' data(iris)
#' X <- as.matrix(iris[,1:4])
#' Y <- iris[,5]
#' A <- between_class_scatter(X, Y, colMeans(X))
#' B <- within_class_scatter(X, Y)
#' ret <- trace_ratio(A,B, ncomp=3)
trace_ratio <- function(A, B, ncomp=2, eps=1e-6, maxiter=100) {
  ## in the the language
  n = nrow(A)
  p = ncomp
  
  ## prepare the initializer
  Vold = qr.Q(qr(matrix(rnorm(n*p),ncol=p)))
  #Vold <- matrix(rnorm(n*n), n, n)
  rhoold = 0
  for (i in 1:maxiter){
    message("iter: ", i)
    if (i > 10) {
      #Vnew   = PRIMME::eigs_sym(A-rhoold*B,p, x0=Vold)$vectors
      Vnew   = PRIMME::eigs_sym(A-rhoold*B,p, x0=Vold)$vectors
      #Vnew   = PRIMME::eigs_sym(A-rhoold*B,p)$vectors
      #Vnew   = RSpectra::eigs_sym(A-rhoold*B,k=p, initvec=Vold[,1])$vectors
      #Vnew <- rsvd::rsvd(A-rhoold*B, k=p)$v
      #Vnew   = irlba::partial_eigen(A-rhoold*B,n=p)$vectors
    } else {
      Vnew   = PRIMME::eigs_sym(A-rhoold*B,p)$vectors

      #Vnew   = RSpectra::eigs_sym(A-rhoold*B,k=p)$vectors
    }
    num <- crossprod(Vnew, A) %*% Vnew
    denom <- crossprod(Vnew, B) %*% Vnew
    rhonew = sum(diag(num)/diag(denom))
    message("rho: ", rhonew)
    rhoinc = abs(rhonew-rhoold)
    Vold   = Vnew
    rhoold = rhonew
    message("delta: ", (rhoinc))
    if (rhoinc < eps){
      break
    }
  }
  
  V <- Vold
  
  values <- sapply(1:ncomp, function(i) {
    drop((V[,i] %*% A %*% V[,i])/ (V[,i] %*% B %*% V[,i]))
    
  })
  
  #values = diag(t(V)%*%A%*%V)/diag(t(V)%*%B%*%V)
  list(vectors=V,values=values)
}
