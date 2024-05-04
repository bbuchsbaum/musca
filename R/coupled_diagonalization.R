#' Construct Laplacian Matrix
#'
#' This function constructs the normalized Laplacian matrix from a given data matrix.
#'
#' @param X A numeric matrix representing the data samples.
#' @param sigma The bandwidth parameter for the Gaussian kernel.
#' @param knn The number of nearest neighbors to consider.
#' @return The normalized Laplacian matrix.
#' @export
construct_laplacian <- function(X, sigma, knn=5) {
  # Check input dimensions
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix.")
  }
  
  g <- neighborweights::graph_weights(X, k=knn, sigma=sigma, neighbor_mode="knn")
  W <- neighborweights::adjacency(g)
  # Compute pairwise similarities using Gaussian kernel
  #W <- exp(-as.matrix(dist(X))^2 / (2 * sigma^2))
  
  # Compute degree matrix
  D <- Matrix::Diagonal(x=rowSums(W))
 
  # Construct normalized Laplacian matrix
  L <- D - W
  Disq_half <- Matrix::Diagonal(x = 1/sqrt(diag(D)))
  L <- Disq_half %*% L %*% Disq_half
  
  return(L)
}

#' Compute Eigenvectors and Eigenvalues
#'
#' This function computes the first k' eigenvectors and eigenvalues of a given Laplacian matrix.
#'
#' @param L The Laplacian matrix.
#' @param k_prime The number of eigenvectors and eigenvalues to compute.
#' @return A list containing the eigenvectors (U_bar) and eigenvalues (Lambda_bar).
#' @export
compute_eigenvectors <- function(L, k_prime) {
  # Check input dimensions
  if ( (!is.matrix(L) && !inherits(L, "Matrix")) || !isSymmetric(L)) {
    stop("L must be a symmetric matrix.")
  }
  
 
  # Compute first k' eigenvectors and eigenvalues
  eig <- RSpectra::eigs_sym(L, k = k_prime, which = "SM")
  
  U_bar <- eig$vectors
  Lambda_bar <- diag(eig$values)
  
  return(list(U_bar = U_bar, Lambda_bar = Lambda_bar))
}

#' Coupled Diagonalization Optimization
#'
#' This function performs the optimization step of the coupled diagonalization method.
#'
#' @param U_bars A list of matrices containing the eigenvectors of the Laplacian matrices.
#' @param Lambda_bars A list of diagonal matrices containing the eigenvalues of the Laplacian matrices.
#' @param F_mats A list of matrices representing the coupling constraints.
#' @param G_mats A list of matrices representing the decoupling constraints.
#' @param k The desired reduced dimension.
#' @param mu_c The coupling strength parameter.
#' @param mu_d The decoupling strength parameter.
#' @return A list of optimal transformation matrices A_i.
optimization <- function(U_bars, Lambda_bars, F_mats, G_mats, k, k_prime, mu_c, mu_d) {
  # Check input dimensions
  if (length(U_bars) != length(Lambda_bars) || length(U_bars) != length(F_mats) || length(U_bars) != length(G_mats)) {
    stop("Input lists must have the same length.")
  }
  
  # Initialize A_list
  A_list <- lapply(U_bars, function(U_i) matrix(rnorm(k_prime*k), nrow = k_prime, ncol = k))
  
  cost_function <- function(A_i_vec, U_i, Lambda_i_bar, F_i, G_i, A_j_list, F_j_list, G_j_list) {
    A_i <- matrix(A_i_vec, nrow = k_prime, ncol = k)
    
    # Compute the first term of the cost function
    term1 <- sum((t(A_i) %*% Lambda_i_bar %*% A_i - diag(diag(t(A_i) %*% Lambda_i_bar %*% A_i)))^2)
    
    # Compute the second term of the cost function
    term2 <- mu_c * sum(sapply(seq_along(A_j_list), function(j) {
      norm(t(F_i) %*% U_i %*% A_i - t(F_j_list[[j]]) %*% U_bars[[j]] %*% A_j_list[[j]], type = "F")^2
    }))
    
    # Compute the third term of the cost function
    term3 <- -mu_d * sum(sapply(seq_along(A_j_list), function(j) {
      norm(t(G_i) %*% U_i %*% A_i - t(G_j_list[[j]]) %*% U_bars[[j]] %*% A_j_list[[j]], type = "F")^2
    }))
    
    return(term1 + term2 + term3)
  }
  
  gradient_function <- function(A_i_vec, U_i, Lambda_i_bar, F_i, G_i, A_j_list, F_j_list, G_j_list) {
    A_i <- matrix(A_i_vec, nrow = k_prime, ncol = k)
    
    # Compute the gradient of the first term
    grad1 <- 4 * (t(A_i) %*% Lambda_i_bar %*% A_i - diag(diag(t(A_i) %*% Lambda_i_bar %*% A_i))) %*% t(A_i) %*% Lambda_i_bar
    #grad1 <- 4 * (Lambda_i_bar %*% A_i %*% t(A_i) %*% Lambda_i_bar %*% A_i - Lambda_i_bar %*% A_i %*% Lambda_i_bar)
    
    
    # Compute the gradient of the second term
    grad2 <- 2 * mu_c * t(U_i) %*% F_i %*% (t(F_i) %*% U_i %*% A_i - rowSums(sapply(seq_along(A_j_list), function(j) {
      t(F_j_list[[j]]) %*% U_bars[[j]] %*% A_j_list[[j]]
    })))
    
    # Compute the gradient of the third term
    grad3 <- -2 * mu_d * t(U_i) %*% G_i %*% (t(G_i) %*% U_i %*% A_i - rowSums(sapply(seq_along(A_j_list), function(j) {
      t(G_j_list[[j]]) %*% U_bars[[j]] %*% A_j_list[[j]]
    })))

    
    # Combine the gradients and reshape to a vector
    return(as.vector(t(grad1) + grad2 + grad3))
  }
  
  # Perform block coordinate optimization
  for (i in seq_along(U_bars)) {
    U_i <- U_bars[[i]]
    Lambda_i_bar <- Lambda_bars[[i]][1:k_prime, 1:k_prime]
    F_i <- F_mats[[i]]
    G_i <- G_mats[[i]]
    
    A_j_list <- A_list[-i]
    F_j_list <- F_mats[-i]
    G_j_list <- G_mats[-i]
    
    opt_result <- optim(
      par = as.vector(A_list[[i]]),
      fn = cost_function,
      gr = gradient_function,
      method = "BFGS",
      U_i = U_i,
      Lambda_i_bar = Lambda_i_bar,
      F_i = F_i,
      G_i = G_i,
      A_j_list = A_j_list,
      F_j_list = F_j_list,
      G_j_list = G_j_list
    )
    
    A_list[[i]] <- matrix(opt_result$par, nrow = k_prime, ncol = k)
  }
  
  return(A_list)
}

#' Coupled Diagonalization
#'
#' This function performs coupled diagonalization of Laplacian matrices for multimodal manifold analysis.
#'
#' @param X_list A list of numeric matrices representing the data samples from each modality.
#' @param F_mats A list of matrices representing the coupling constraints.
#' @param G_mats A list of matrices representing the decoupling constraints.
#' @param k The desired reduced dimension.
#' @param k_prime The number of eigenvectors and eigenvalues to compute.
#' @param mu_c The coupling strength parameter.
#' @param mu_d The decoupling strength parameter.
#' @param sigma The bandwidth parameter for the Gaussian kernel.
#' @return A list of coupled basis vectors V_i for each modality.
#' @export
coupled_diagonalization <- function(X_list, F_mats, G_mats, k, k_prime, mu_c, mu_d, sigma=1, knn=.1*nrow(X_list[[1]]), lambda=1e-4) {
  # Check input dimensions
  if (length(X_list) != length(F_mats) || length(X_list) != length(G_mats)) {
    stop("Input lists must have the same length.")
  }
  
  X_list <- multidesign::multiblock(X_list)
  bind <- lapply(1:length(X_list), function(i) block_indices(X_list, i))

  # Construct Laplacian matrices
  L_list <- lapply(X_list, construct_laplacian, sigma = sigma, knn=knn)
  
  # Compute eigenvectors and eigenvalues
  eig_list <- lapply(L_list, compute_eigenvectors, k_prime = k_prime)
  U_bars <- lapply(eig_list, function(x) x$U_bar)
  Lambda_bars <- lapply(eig_list, function(x) x$Lambda_bar)
  
  # Perform optimization
  A_list <- optimization(U_bars, Lambda_bars, F_mats, G_mats, k, k_prime, mu_c, mu_d)
  
  # Construct coupled basis vectors
  V_list <- mapply(function(U_i, A_i) U_i %*% A_i, U_bars, A_list, SIMPLIFY = FALSE)
  #Z <- Matrix::bdiag(V_list)
  
  vectors <- lapply(seq_along(X_list), function(i) {
    Z <- V_list[[i]]
    rfit <- glmnet(Z, X_list[[i]], family = "mgaussian", alpha=0, lambda=lambda)
    cfs <- coef(rfit)
    cfs <- do.call(cbind, cfs)[-1,,drop=FALSE]
  })
    
  vectors <- do.call(cbind, vectors)

  multivarious::multiblock_projector(
    v=vectors,
    V_list=V_list,
    preproc=multivarious::pass(),
    block_indices=bind,
    classes="coupled_diagonalization"
  )
}