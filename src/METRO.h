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