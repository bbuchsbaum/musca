
#' Barycentric Discriminant Analysis
#' 
#' A component technique that maximizes the between group variance over a set of variables. 
#' 
#' @details 
#' 
#' The \code{subject} argument is used to model multilevel structure. When `rescomp` > 0, then the group barycenters
#' will be subtracted from each subject's data and a residual pca-lda will be computed to estimate components of variation
#' that deviate from the group barycenters. This can be used to capture individual level variation that can be useful
#' for modeling idiosyncratic structure in the multilevel data.
#' 
#' @param data the object (e.g. `multidesign`) containing the X and Y data.
#' @param y the categorical response variable
#' @param subject the subject variable
#' @param ncomp number of components to estimate
#' @param preproc pre-processing function, defaults to \code{center}
#' 
#' @references
#' Abdi, H., Williams, L. J., & Bera, M. (2017). Barycentric discriminant analysis. \emph{Encyclopedia of Social Network Analysis and Mining}, 1-20.
#' 
#' @examples 
#' 
#' X <- matrix(rnorm(50*100), 50, 100)
#' Y <- tibble(condition=rep(letters[1:5], 10), subject=rep(1:5, each=10))
#' 
#' des <- multivarious::multidesign(X,Y)
#' bres1 <- bada(des, condition, subject=subject, rescomp=0)
#' bres2 <- bada(des, condition, subject=subject, resdim=20, rescomp=5)
bada <- function(data, y, subject, preproc, ncomp) UseMethod("bada")


#' Statis for covariance matrices
#' 
#' @param data a set of covariance matrices
#' @param the number of componentss to compute
#' @param normalize whether to normalize so that the sums-of-squares of each matrix is 1.
#' @param dcenter whether to double center the covariance matrix.
#' @export
covstatis <- function(data, ncomp, normalize, dcenter, ...) UseMethod("covstatis")


#' Multiple Factor Analysis
#' 
#' principal component analysis for multiple blocks of data collected over th same instances
#' 
#' 
#' @inheritParams bada
#' @param normalization the normalization method: MFA, RV, RV-MFA, or None (see details).
#' @param A custom weight matrix for the columns
#' @param M custom weight matrix for the rows
#' @export
#' 
#' @references
#' Abdi, H., Williams, L. J., & Valentin, D. (2013). Multiple factor analysis: principal component analysis for multitable and multiblock data sets. 
#' \emph{Wiley Interdisciplinary reviews: computational statistics}, 5(2), 149-179.
#' @export
mfa <- function(data, preproc, ncomp, normalization, A, M) UseMethod("mfa")


#' @inheritParams bada
#' @param knn the number of nearest neighbors to use in construction of local similarity graph
#' @param sigma the bandwidth of the heat kernel used during local similarity graph construction
#' @param u the relative weight to place on local topology preservation vs cross-block alignment (0 = all cross-block, 1 = all local)
#' @param kernel the kernel function (e.g. `rbfdot` or any other kernel function from `kernlab` package)
#' @param sample_frac the fraction of samples to use during estimation 
#' @param dweight the weight to place on the "repulsion" graph
#' @param 
kema <- function(data, y, preproc, ncomp, knn, sigma, u, kernel, sample_frac, dweight, ...) UseMethod("kema")

#' project covariance matrix
#' 
#' project a new covariance matrix onto a component fit
#' 
#' @param x the model fit
#' @param new_data the covariance matrix
#' @param ... extra args
#' @export
project_cov <- function(x, new_data, ...) UseMethod("project_cov")


#' @importFrom multivarious project
#' @export
#' @rdname project
multivarious::project


