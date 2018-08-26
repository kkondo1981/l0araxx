#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List l0araxxC(arma::mat X, arma::vec y, arma::vec weights, arma::vec offset,
              String family, double lambda, int maxit, double eps)
{
  // initialization
  int n = X.n_rows;
  int m = X.n_cols;
  arma::mat beta = zeros(m, 1);
  arma::mat old_beta = zeros(m, 1);
  arma::mat beta_hist = zeros(m, maxit + 1);
  arma::mat DXt(m, n);
  arma::mat yy(n, 1);
  arma::mat Xbeta(n, 1);
  arma::mat mu(n, 1);
  arma::mat V(n, n);
  arma::mat z(n, 1);
  arma::mat Pm(m, m);

  if (family == "gaussian") {
    beta(0) = mean(y);
  } else if (family == "poisson") {
    beta(0) = log(mean(y/exp(offset)));
  } else if (family == "gamma") {
    beta(0) = 1 / mean(y);
  } else if (family == "gamma(log)") {
    beta(0) = log(mean(y));
  }
  beta_hist.col(0) = beta;

  Pm = eye(m,m);
  Pm(0, 0) = 0;
  DXt = trans(X);

  int iter;
  Rcout << "Iteration... " << std::flush;
  for (iter = 1;; ++iter) {
    Rcout << '+' << std::flush;

    old_beta = beta;
    Xbeta = X * beta;

    if (family == "gaussian") {
      mu = Xbeta;
      yy = weights % (y - mu);
      V = diagmat(vec(weights));
    } else if (family == "poisson") {
      mu = exp(Xbeta + offset);
      yy = weights % (y - mu);
      V = diagmat(vec(weights % mu));
    } else if (family == "gamma") {
      mu = 1 / Xbeta;
      yy = -weights % (y - mu);
      V = diagmat(vec(weights % mu % mu));
    } else if (family == "gamma(log)") {
      mu = exp(Xbeta);
      yy = weights % (y / mu - 1);
      V = diagmat(vec(weights % (y / mu)));
    }

    z = V * Xbeta + yy;
    beta = solve(DXt * V * X + lambda * Pm, DXt * z);
    beta_hist.col(iter) = beta;

    DXt = trans(repmat(trans(beta % beta), n, 1) % X);

    // check for convergence
    if (norm(beta - old_beta, 2) < eps) {
      break;
    } else if (iter >= maxit) {
      warning("Did not converge. Increase maxit.");
      break;
    }
  }
  Rcout << std::endl;

  for (int i = 0; i < m; ++i) {
    if (std::abs(beta(i, 0)) < 1e-3) {
      beta(i,0) = 0;
    }
  }

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j <= iter; ++j) {
      if (std::abs(beta_hist(i, j)) < 1e-3) {
        beta_hist(i, j) = 0;
      }
    }
  }

  return List::create(Named("beta") = beta,
                      Named("iter") = iter,
                      Named("beta_hist") = beta_hist);
}
