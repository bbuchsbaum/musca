#' Align two graphs using the GRASP algorithm
#'
#' @param G1 An igraph object representing the first graph.
#' @param G2 An igraph object representing the second graph.
#' @param k The number of eigenvectors to use for alignment (default: 20).
#' @param q The number of corresponding functions to use (default: 100).
#' @param t A vector of time steps for the heat kernel (default: seq(0.1, 50, length.out = q)).
#' @param mu The regularization parameter for the optimization problem (default: 0.132).
#'
#' @return A list containing the following elements:
#'   - alignment: A vector of length |V1| indicating the aligned node indices in G2 for each node in G1.
#'   - P: The permutation matrix representing the node alignment.
#'   - M: The alignment matrix for the eigenvectors.
#'   - C: The diagonal mapping matrix for the corresponding functions.
#' 
GRASP_align <- function(G1, G2, k = 20, q = 100, t = seq(0.1, 50, length.out = q), mu = 0.132) {
  n <- igraph::vcount(G1)  # Number of nodes in the graphs
  
  # Step 1: Compute eigenvectors and eigenvalues
  L1 <- igraph::laplacian_matrix(G1)
  L2 <- igraph::laplacian_matrix(G2)
  
  ei1 <- eigen(L1)
  ei2 <- eigen(L2)
  
  phi <- ei1$vectors
  lambda1 <- diag(ei1$values)
  
  psi <- ei2$vectors
  lambda2 <- diag(ei2$values)
  
  # Step 2: Compute corresponding functions
  F <- matrix(0, nrow = n, ncol = q)
  G <- matrix(0, nrow = n, ncol = q)
  
  for (i in 1:q) {
    F[, i] <- colSums(phi * exp(-t[i] * diag(lambda1)) * phi)
    G[, i] <- colSums(psi * exp(-t[i] * diag(lambda2)) * psi)
  }
  
  # Step 3: Base alignment
  phi_k <- phi[, 1:k]
  psi_k <- psi[, 1:k]
  lambda2_k <- lambda2[1:k, 1:k]
  
  M_init <- diag(rep(1, k))
  
  objective_fn <- function(M) {
    off_diag_sum <- sum(diag(M %*% lambda2_k %*% t(M)))
    alignment_term <- norm(t(F) %*% phi_k - t(G) %*% psi_k %*% M, type = "F")^2
    return(off_diag_sum + mu * alignment_term)
  }
  
  gradient_fn <- function(M) {
    off_diag_grad <- 2 * (lambda2_k %*% M - diag(diag(M %*% lambda2_k %*% t(M))))
    alignment_grad <- 2 * mu * t(psi_k) %*% G %*% (t(G) %*% psi_k %*% M - t(F) %*% phi_k)
    return(off_diag_grad + alignment_grad)
  }
  
  opt_result <- optim(par = M_init, fn = objective_fn, gr = gradient_fn, method = "BFGS")
  
  M <- matrix(opt_result$par, nrow = k, ncol = k)
  psi_hat <- psi_k %*% M
  
  # Step 4: Calculate mapping matrix
  A <- matrix(0, nrow = q, ncol = k)
  b <- matrix(0, nrow = q, ncol = 1)
  
  for (i in 1:q) {
    A[i, ] <- diag(G[, i] %*% psi_hat)
    b[i, 1] <- sum(phi_k * matrix(F[, i], ncol = 1))
  }
  
  c <- solve(t(A) %*% A) %*% t(A) %*% b
  C <- diag(as.vector(c))
  
  # Step 5: Node alignment
  alignment <- linear_assignment(phi_k, C %*% t(psi_hat))
  
  P <- matrix(0, nrow = n, ncol = n)
  P[cbind(1:n, alignment)] <- 1
  
  return(list(alignment = alignment, P = P, M = M, C = C))
}

#' Compute the normalized Laplacian matrix of a graph
#'
#' @param G An igraph object representing the graph.
#'
#' @return The normalized Laplacian matrix of the graph.
normalized_laplacian <- function(G) {
  A <- igraph::as_adjacency_matrix(G)
  D <- diag(rowSums(A))
  
  L <- D - A
  D_inv_sqrt <- diag(1 / sqrt(diag(D)))
  
  L_norm <- D_inv_sqrt %*% L %*% D_inv_sqrt
  
  return(L_norm)
}

#' Perform linear assignment using the Jonker-Volgenant algorithm
#'
#' @param C1 The first matrix of size (n x k) for assignment.
#' @param C2 The second matrix of size (n x k) for assignment.
#'
#' @return A vector of length n indicating the assigned indices.
linear_assignment <- function(C1, C2) {
  cost_matrix <- matrix(0, nrow = nrow(C1), ncol = nrow(C2))
  
  for (i in 1:nrow(C1)) {
    for (j in 1:nrow(C2)) {
      cost_matrix[i, j] <- sqrt(sum((C1[i, ] - C2[j, ])^2))
    }
  }
  
  assignment <- clue::solve_LSAP(cost_matrix, maximum = FALSE)
  
  return(assignment)
}