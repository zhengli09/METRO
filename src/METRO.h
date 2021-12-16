#ifndef _METRO_H_
#define _METRO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List METROSummaryStats(
  const Rcpp::NumericMatrix betaeQTLin, // p x M matrix of eQTL marginal effect size estimates in M populations
  const Rcpp::NumericVector betaGWASin, // p x 1 vector of GWAS marginal effect size estimates
  const Rcpp::List Dzin, // length M list of p x p SNP correlation matrix (LD) in each population
  const Rcpp::NumericMatrix Din, // p x p matrix of SNP correlation matrix (LD) in GWAS data
  const Rcpp::NumericVector nzin, // M x 1 vector of sample size in M populations
  const double n, // sample size in GWAS data
  const double nu, // shrinkage parameter for SNP correlation matrix
  const double hthre = 1e-3, // threshold of expression heritability
  const int maxIter = 1000, // maximum iterations for EM algorithm
  const double tol = 1e-3, // convergence threshold for EM algorithm
  const bool verbose = false // whether to print intermediate information
  );

// Extension to adjust for horizontal pleiotropic effect
Rcpp::List METROSummaryStatsPleio(
  const Rcpp::NumericMatrix betaeQTLin,
  const Rcpp::NumericVector betaGWASin,
  const Rcpp::List Dzin,
  const Rcpp::NumericMatrix Din,
  const Rcpp::NumericVector nzin,
  const double n,
  const double nu,
  const double hthre = 1e-3,
  const int maxIter = 1000,
  const double tol = 1e-3,
  const bool verbose = false
  );

// Extension to adjust for additional covariates (e.g., ancestry proportions)
Rcpp::List METROCovars(
  const Rcpp::NumericMatrix betaeQTLin,
  const Rcpp::NumericVector betaGWASin,
  const Rcpp::List Dzin,
  const Rcpp::NumericMatrix Din,
  const Rcpp::NumericVector nzin,
  const double n,
  const Rcpp::NumericMatrix covGCovarsin, // p x p0 covariance matrix between GWAS genoypte G
                                          // and covariates matrix X0 (i.e., G^T * X0 / (n - 1))
  const Rcpp::NumericMatrix covYCovarsin, // p0 x 1 covariance matrix between GWAS phenotype y
                                          // and covariates matrix X0 (i.e., X0^T * y / (n - 1))
  const Rcpp::NumericMatrix covCovarsin, // p0 x p0 covariance matrix of covariates matrix X0
  const double nu,
  const double hthre = 1e-3,
  const int maxIter = 1000,
  const double tol = 1e-3,
  const bool verbose = false
  );

Rcpp::List METRO_EM_Algorithm(
  // input data
  const arma::vec &betaeQTL,
  const arma::vec &betaGWAS,
  const arma::cube &Dz,
  const arma::mat &D,
  const arma::vec nz,
  const double n,
  // outputs
  arma::vec &sigma2m, // residual variance in the m'th expression data
  arma::vec &sigma2beta, // variance of SNP effects
  // hyperparameters
  const std::string mode, // either "null" or "alternative"
  const std::string constrain,  // either "positive" or "negative" to constrain 
                                // the sign of alpha under the alternative
  const int maxIter, // maximum iterations for EM algorithm
  const double tol, // convergence threshold for EM algorithm
  const bool verbose // whether to print out intermediate information
  );

// Extension to adjust for horizontal pleiotropic effect
Rcpp::List METRO_Pleio_EM_Algorithm(
  // input data
  const arma::vec &betaeQTL,
  const arma::vec &betaGWAS,
  const arma::cube &Dz,
  const arma::mat &D,
  const arma::vec nz,
  const double n,
  arma::vec &sigma2m,
  arma::vec &sigma2beta,
  const std::string mode,
  const std::string constrain,
  const int maxIter,
  const double tol,
  const bool verbose
  );

// Extension to adjust for additional covariates (e.g., ancestry proportions)
Rcpp::List METRO_Covars_EM_Algorithm(
  const arma::vec &betaeQTL,
  const arma::vec &betaGWAS,
  const arma::cube &Dz,
  const arma::mat &D,
  const arma::vec nz,
  const double n,
  const arma::mat &covGCovars,
  const arma::mat &covYCovars,
  const arma::mat &covCovars,
  arma::vec &sigma2m,
  arma::vec &sigma2beta,
  const std::string mode,
  const std::string constrain,
  const int maxIter,
  const double tol,
  const bool verbose
  );

arma::mat blockDiagonal(
  const arma::cube X // m x n x l cube and each slice is a block matrix
  );

arma::mat blockRowCol(
  const arma::cube X, // m x n x l cube and each slice is a block matrix
  std::string mode // either "row" or "column"
  );

void constrainVec(
  arma::vec &alpha,
  std::string sign // either "positive" or "negative"
  );
#endif