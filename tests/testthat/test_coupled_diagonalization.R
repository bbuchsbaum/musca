library(testthat)

# Test case 1: Identical Laplacian matrices
test_that("Coupled diagonalization works for identical Laplacian matrices", {
  n <- 100
  L1 <- L2 <- matrix(rnorm(n^2), nrow = n, ncol = n)
  L1 <- L1 + t(L1)
  L2 <- L2 + t(L2)
  X_list <- list(X1 = L1, X2 = L2)
  F_mats <- list(F1 = diag(n), F2 = diag(n))
  G_mats <- list(G1 = matrix(0, nrow = n, ncol = n), G2 = matrix(0, nrow = n, ncol = n))
  k <- 10
  k_prime <- 20
  mu_c <- 1
  mu_d <- 0
  sigma <- 1
  
  V_list <- coupled_diagonalization(X_list, F_mats, G_mats, k, k_prime, mu_c, mu_d, sigma, knn=13)
  
  expect_equal(length(V_list), 2)
  expect_equal(ncol(V_list[[1]]), k)
  expect_equal(ncol(V_list[[2]]), k)
  expect_true(all(abs(V_list[[1]] - V_list[[2]]) < 1e-6))
})

# Test case 2: Commuting Laplacian matrices
test_that("Coupled diagonalization works for commuting Laplacian matrices", {
  n <- 100
  D1 <- diag(sort(runif(n)))
  D2 <- diag(sort(runif(n)))
  Q <- qr.Q(qr(matrix(rnorm(n^2), nrow = n, ncol = n)))
  L1 <- Q %*% D1 %*% t(Q)
  L2 <- Q %*% D2 %*% t(Q)
  X_list <- list(X1 = L1, X2 = L2)
  F_mats <- list(F1 = diag(n), F2 = diag(n))
  G_mats <- list(G1 = matrix(0, nrow = n, ncol = n), G2 = matrix(0, nrow = n, ncol = n))
  k <- 10
  k_prime <- 20
  mu_c <- 1
  mu_d <- 0
  sigma <- 1
  
  V_list <- coupled_diagonalization(X_list, F_mats, G_mats, k, k_prime, mu_c, mu_d, sigma=.1, knn=50)
  
  testthat::expect_equal(length(V_list), 2)
  testthat::expect_equal(ncol(V_list[[1]]), k)
  testthat::expect_equal(ncol(V_list[[2]]), k)

})



