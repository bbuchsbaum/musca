
#' Linear Similarity Embedding using ADAM optimization
#'
#' @param X Input data matrix (n x d)
#' @param T Target similarity matrix (n x n)
#' @param M Mask matrix (n x n)
#' @param sigma_P Scaling parameter for the Gaussian kernel
#' @param ncomp Number of reduced dimensions
#' @param ap Weighting parameter for the orthogonality regularizer
#' @param maxit Maximum number of iterations for optimization
#' @param tol Tolerance for optimization
#' @return Optimized weight matrix W (d x m)
linear_sim_embed <- function(X, T, M, sigma_P=1, ncomp=2, ap=.2, maxit = 100, tol = 1e-8, batch_size=.1, use_cpp=TRUE) {
  # Perform PCA to get initial W
  preproc <- multivarious::standardize()
  proc <- prep(preproc, X)
  Xs <- init_transform(proc, X)
  
  w <- if (use_cpp) {
    ret <- linear_sim_embed_cpp(Xs, as.matrix(T), as.matrix(M), as.numeric(sigma_P), as.integer(ncomp), 
                         as.numeric(alpha_p), as.integer(maxit), as.numeric(tol), as.numeric(batch_size))
    ret$W
  } else {
    pca_result <- prcomp(Xs)
    initial_W <- pca_result$rotation[, 1:ncomp]  # Take the first m principal components
    optimize_W(Xs, T, M, initial_W, sigma_P, ap, maxit, tol)
  }
  
  sc <- Xs %*% w
  sdev <- apply(sc,2,sd)
  ret <- multivarious::bi_projector(as.matrix(w), sc, sdev, proc, target_sim=T, mask=M, sigma_P=sigma_P, ap=ap, classes="linear_sim_embedding")
}


#' Main function to perform optimization using ADAM
#'
#' @param X Input data matrix (n x d)
#' @param T Target similarity matrix (n x n)
#' @param M Mask matrix (n x n)
#' @param initial_W Initial weight matrix (d x m)
#' @param sigma_P Scaling parameter for the Gaussian kernel
#' @param ap Weighting parameter for the orthogonality regularizer
#' @param maxit Maximum number of iterations for optimization
#' @param tol Tolerance for optimization
#' @return Optimized weight matrix W (d x m)
optimize_W <- function(X, T, M, initial_W, sigma_P=1, ap=.2, maxit = 200, tol = 1e-8) {
  # Convert initial_W to a vector for optimg
  initial_W_vec <- as.vector(initial_W)
  d <- ncol(X)
  m <- ncol(initial_W)
  
  
  P <<- NULL
  Y <<- NULL
  W_mat <<- NULL
  count <- 1
  objective_function <- function(W, X, T, M, sigma_P, d, m) {
    count <<- count+1
    W_mat <<- matrix(W, nrow = ncol(X), ncol = ncol(initial_W))
    
    Y <<- X %*% W_mat
    
    P <<- compute_P(Y, sigma_P)
    
    Js <- sum((P - T)^2 * M) * 1/(2*sum(M))
    Jp <- sum((t(W_mat) %*% W_mat - diag(ncol(W_mat)))^2) * 1/(2*ncol(W_mat)^2)
    
    obj <- (2-ap) *Js + ap*Jp
    
  
    #if (count  == 1 || count %% 50 == 0) {
    #  print(W_mat)
    #  print(paste("objective=", obj, "Js=", Js, "Jp=", Jp))
    #}
    
    return(obj)
  }
  
  gradient_function <- function(W_vec, X, T, M, sigma_P, d, m) {
    if (is.null(W_mat)) {
      W_mat <<- matrix(W_vec, nrow = d, ncol = m)  
    }
    
    # Reshape vector back to matrix
    grad_Js <- compute_gradient_Js(W_mat, X, T, M, sigma_P)
    
    grad_Jp <- compute_grad_Jp(W_mat)
    total_grad <- (2-ap)*grad_Js + ap*grad_Jp
    return(as.vector(total_grad))  # Convert gradient matrix back to vector for optimg
  }
  
  # Call to optimg for optimization using ADAM
   # result <- optimg(par = initial_W_vec,
   #                  fn = function(W_vec) objective_function(W_vec, X, T, M, sigma_P, d, m),
   #                  gr = function(fn, startvalue, Interval) gradient_function(startvalue, X, T, M, sigma_P, d, m),
   #                  method = "ADAM",
   #                  maxit = maxit,
   #                  tol = tol,
   #                  Interval=1e-3,
   #                  verbose = TRUE)

  result <- optim(par = initial_W_vec,
                   fn = function(W_vec) objective_function(W_vec, X, T, M, sigma_P, d, m),
                   gr = function(W_vec) gradient_function(W_vec, X, T, M, sigma_P, d, m),
                   method="L-BFGS-B",
                   control = list(maxit = maxit, trace = 1))

  print(result)
  # Reshape the optimized parameters back to matrix form
  optimized_W <- matrix(result$par, nrow = d, ncol = m)
  return(optimized_W)
}

#' Compute the similarity matrix P
#'
#' @param W Matrix of weights (d x m)
#' @param X Matrix of input data (n x d)
#' @param sigma_P Scaling parameter for the Gaussian kernel
#' @return Matrix P of similarity scores (n x n)
compute_P <- function(Y, sigma_P) {
  dists <- as.matrix(dist(Y)^2)  # Squared Euclidean distances (n x n)
  P <- exp(-dists / sigma_P)
  zapsmall(P)# Gaussian kernel similarity (n x n)
  return(P)
}


compute_gradient_Js <- function(W, X, T, M, sigma_P) {

  n <- nrow(X)
  d <- ncol(X)
  m <- ncol(W)
  
  #Y <- X %*% W  # Projected data (n x m)
  if (is.null(Y)) {
    Y <<- X %*% W
  }
  
  P <<- compute_P(Y, sigma_P)  # Similarity matrix (n x n)
  
  # Matrix to hold the gradient of Js
  grad_Js <- matrix(0, nrow = d, ncol = m)
  
  for (i in 1:n) {
    for (j in 1:n) {
      Y_diff <- Y[i,] - Y[j,]  # (1 x m)
      X_diff <- X[i,] - X[j,]  # (1 x d)
      grad_component <- -2 / sigma_P * P[i,j] * outer(X_diff, Y_diff)  # (d x m)
      grad_Js <- grad_Js + M[i,j] * (P[i,j] - T[i,j]) * grad_component
    }
  }
  
  grad_Js <- grad_Js / sum(M)  # Normalize by sum of M
  return(grad_Js)
}

# Compute the gradient using the fast version
compute_gradient_Js_fast <- function(W, X, T, M, sigma_P) {
  N <- nrow(X)
 
  #Y <- X %*% W  # Projected data (n x m)
  if (is.null(Y)) {
    Y <<- X %*% W
  }
  
  P <<- compute_P(Y, sigma_P)
  i <- rep(1:N, each=N)
  j <- rep(1:N, N)
  
  
  X_diff <- X[i,] - X[j,]
  Y_diff <- Y[i,] - Y[j,]
  P_diff <- P[cbind(i,j)] - T[cbind(i,j)]
  
  dJp_dW <- 2/sigma_P * t(X_diff) %*% (Y_diff * P[cbind(i,j)] * M[cbind(i,j)] * P_diff)
  -dJp_dW / sum(M)
}



#' Compute the gradient of the orthogonality regularizer Jp
#'
#' @param W Matrix of weights (d x m)
#' @param m Number of reduced dimensions
#' @return Gradient of Jp with respect to W (d x m)
#compute_grad_Jp <- function(W, m) {
#  WT_W <- t(W) %*% W  # W'W (m x m)
#  WT_W_minus_I <- WT_W - diag(m)  # Subtract the identity matrix (m x m)
#  grad_Jp <- 2 / m^2 * W %*% WT_W_minus_I  # Gradient computation (d x m)
#  return(grad_Jp)
#}

compute_grad_Jp <- function(W) {
  #dJp_dW <- 2/m^2 * W %*% (t(W) %*% W - diag(m))
  m <- ncol(W)
  WtW <- t(W) %*% W
  diag(WtW) <- diag(WtW) - 1
  grad_Jp <- 2 / m^2 * W %*% WtW
  return(grad_Jp)
}






