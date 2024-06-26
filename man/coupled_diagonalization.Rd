% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coupled_diagonalization.R
\name{coupled_diagonalization}
\alias{coupled_diagonalization}
\title{Coupled Diagonalization}
\usage{
coupled_diagonalization(
  X_list,
  F_mats,
  G_mats,
  k,
  k_prime,
  mu_c,
  mu_d,
  sigma = 1,
  knn = 0.1 * nrow(X_list[[1]]),
  lambda = 1e-04
)
}
\arguments{
\item{X_list}{A list of numeric matrices representing the data samples from each modality.}

\item{F_mats}{A list of matrices representing the coupling constraints.}

\item{G_mats}{A list of matrices representing the decoupling constraints.}

\item{k}{The desired reduced dimension.}

\item{k_prime}{The number of eigenvectors and eigenvalues to compute.}

\item{mu_c}{The coupling strength parameter.}

\item{mu_d}{The decoupling strength parameter.}

\item{sigma}{The bandwidth parameter for the Gaussian kernel.}
}
\value{
A list of coupled basis vectors V_i for each modality.
}
\description{
This function performs coupled diagonalization of Laplacian matrices for multimodal manifold analysis.
}
