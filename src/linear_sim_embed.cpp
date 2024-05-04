#include <RcppArmadillo.h>
#include <roptim.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
class LinearSimEmbed : public roptim::Functor {
public:
  LinearSimEmbed(const arma::mat& X, const arma::mat& T, const arma::mat& M, const double sigma_P, const int ncomp, const double alpha_p)
    : X_(X), T_(T), M_(M), sigma_P_(sigma_P), ncomp_(ncomp), alpha_p_(alpha_p) {}
  
  double operator()(const arma::vec& w) override {
    arma::mat W = arma::reshape(w, X_.n_cols, ncomp_);
    arma::mat Y = X_ * W;
    arma::mat P = pairwise_distances(Y);
    double Js = arma::accu(M_ % square(P - T_)) / (2 * arma::accu(M_));
    double Jp = arma::accu(square(W.t() * W - arma::eye(W.n_cols, W.n_cols))) / (2 * std::pow(W.n_cols, 2));
    //Rcpp::Rcout << "Objective function value: " << (2 - alpha_p_) * Js + alpha_p_ * Jp << std::endl;
    //Rcpp::Rcout << "W:\n" << W << std::endl;
    return (2 - alpha_p_) * Js + alpha_p_ * Jp;
   
  }
  
  void Gradient(const arma::vec& w, arma::vec& grad) override {
    //Rcpp::Rcout << "Gradient called: " << std::endl;
    arma::mat W = arma::reshape(w, X_.n_cols, ncomp_);
    arma::mat Y = X_ * W;
    arma::mat P = pairwise_distances(Y);
    
    arma::mat grad_Js = arma::zeros(X_.n_cols, ncomp_);
    
    for (size_t i = 0; i < X_.n_rows; ++i) {
      for (size_t j = 0; j < X_.n_rows; ++j) {
        arma::rowvec Y_diff = Y.row(i) - Y.row(j);
        arma::rowvec X_diff = X_.row(i) - X_.row(j);
        arma::mat grad_component = -2 / sigma_P_ * P(i, j) * X_diff.t() * Y_diff;
        grad_Js += M_(i, j) * (P(i, j) - T_(i, j)) * grad_component;
      }
    }
    
    grad_Js /= arma::accu(M_);
    
    arma::mat WtW = W.t() * W;
    arma::mat grad_Jp = (2 / std::pow(W.n_cols, 2)) * W * (WtW - arma::eye(W.n_cols, W.n_cols));
    
    grad = arma::vectorise((2 - alpha_p_) * grad_Js + alpha_p_ * grad_Jp);
    
    //Rcpp::Rcout << "Gradient norm: " << arma::norm(grad) << std::endl;
  }
  
  // void Gradient(const arma::vec& w, arma::vec& grad) override {
  //   Rcpp::Rcout << "Gradientcalled: " << std::endl;
  //   arma::mat W = arma::reshape(w, X_.n_cols, m_);
  //   arma::mat Y = X_ * W;
  //   arma::mat P = pairwise_distances(Y);
  //   
  //   arma::mat dJs_dW = (2 / sigma_P_) * X_.t() * (P % M_ % (P - T_)) * Y;
  //   dJs_dW /= arma::accu(M_);
  //   
  //   arma::mat WtW = W.t() * W;
  //   arma::mat dJp_dW = (2 / std::pow(W.n_cols, 2)) * W * (WtW - arma::eye(W.n_cols, W.n_cols));
  //   Rcpp::Rcout << "Gradient norm: " << arma::norm(grad);
  //   grad = arma::vectorise((2 - alpha_p_) * dJs_dW + alpha_p_ * dJp_dW);
  //   
  // }
  
  
  
private:
  arma::mat X_;
  arma::mat T_;
  arma::mat M_;
  int ncomp_;
  double sigma_P_;
  double alpha_p_;
  
  arma::mat pairwise_distances(const arma::mat& Y) {
    arma::mat D(Y.n_rows, Y.n_rows);
    for (size_t i = 0; i < Y.n_rows; ++i) {
      for (size_t j = 0; j < Y.n_rows; ++j) {
        D(i, j) = std::exp(-std::pow(arma::norm(Y.row(i) - Y.row(j)), 2) / sigma_P_);
      }
    }
    return D;
  }
};

// [[Rcpp::export]]
Rcpp::List linear_sim_embed_cpp(const arma::mat& X, const arma::mat& T, const arma::mat& M,
                                double sigma_P, int ncomp, double alpha_p, int maxit=200, double tol=1e-8, double batch_size=1) {
  // Perform PCA to get initial W
  arma::mat U, V;
  arma::vec s;
  arma::svd(U, s, V, arma::cov(X));
  arma::mat initial_W = V.cols(0, ncomp - 1);
  //Rcpp::Rcout << "Initial W:\n" << initial_W << std::endl;
  
  

  // Create an instance of LinearSimEmbed
  LinearSimEmbed task(X, T, M, sigma_P, ncomp, alpha_p);
    
  // Create an instance of Roptim with BFGS method
  //roptim::Roptim<LinearSimEmbed> opt("L-BFGS-B");
  roptim::Roptim<LinearSimEmbed> opt("L-BFGS-B");
    
  // Set optimization control parameters
  opt.control.maxit = maxit;
  opt.control.trace = 1;
  opt.control.reltol = tol;
  // Perform optimization
  arma::vec initial_W_vec = arma::vectorise(initial_W);
  opt.minimize(task, initial_W_vec);
    
  // Extract optimized W
  arma::mat optimized_W = arma::reshape(opt.par(), X.n_cols, ncomp);
    
  return Rcpp::List::create(
    Rcpp::Named("W") = optimized_W,
    Rcpp::Named("convergence") = opt.convergence(),
    Rcpp::Named("message") = opt.message()
  );
  
}