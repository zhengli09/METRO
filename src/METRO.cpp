// Author: Zheng Li
// Date: 2021-01-05
// Multiple ancEstry group TRanscriptOme-wide analysis (METRO)
// Note:
// 1. Please refer to supplementary file for detailed derivations
// and parameter notations.
// 2. Please check header file for function specifications.

#include <RcppArmadillo.h>
#include <string>
#include <chrono>
#include "METRO.h"
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List METROSummaryStats(
  const Rcpp::NumericMatrix betaeQTLin,
  const Rcpp::NumericVector betaGWASin,
  const Rcpp::List Dzin,
  const Rcpp::NumericMatrix Din,
  const Rcpp::NumericVector nzin,
  const double n, 
  const double nu,
  const double hthre,
  const int maxIter,
  const double tol,
  const bool verbose
  )
{
  // starting check
  cout << "Starting METRO..." << endl;
  auto start = chrono::steady_clock::now();
  int M = betaeQTLin.ncol(); // number of populations
  int p = betaeQTLin.nrow(); // number of SNPs
  if(verbose)
  {
    cout << "***** info *****" << endl;
    cout << "  - Handling data with " << p << " SNPs" << endl;
    cout << "  - Handling data with " << M << " expression studies " << endl;
  }

  // input processing
  arma::vec nz = Rcpp::as<arma::vec>(nzin);
  arma::vec betaeQTL = Rcpp::as<arma::vec>(betaeQTLin);
  arma::vec betaGWAS = Rcpp::as<arma::vec>(betaGWASin);
  arma::cube Dz (p, p, M);
  for(int m = 0; m < M; m++)
  {
    Dz.slice(m) = Rcpp::as<arma::mat>(Dzin(m));
  }
  arma::mat D = Rcpp::as<arma::mat>(Din);
  D = nu * D + (1.0 - nu) * arma::mat(p, p, fill::eye);

  // outputs
  // define sigma2m and sigma2beta here to provide the estimates under null
  // as initial values in the alternative scenario.
  arma::vec sigma2m(M, fill::ones);
  arma::vec sigma2beta(M);
  sigma2beta.fill(0.01 / p);
  Rcpp::List nullRes;
  Rcpp::List altResPos;
  Rcpp::List altResNeg;
  Rcpp::List altRes;
  Rcpp::List output;
  Rcpp::Function pchisq("pchisq");
  double LRTStat;
  double pvalueLRT;
  arma::vec h2m(M); // heritability of gene expression
  double h2y; // heritability of GWAS outcome
  double alpha; // gene effect on GWAS outcome
  arma::vec w(M); // population weights
  arma::vec beta(p, fill::zeros); // SNP effects on gene expression in GWAS
  double zscore; // TWAS zscore from two stage model

  // METRO algorithm
  // 1. Estimation under the null
  nullRes = METRO_EM_Algorithm(
    betaeQTL, betaGWAS, Dz, D, nz, n, sigma2m, sigma2beta,
    "null", "no", maxIter, tol, verbose);

  // 2.Estimation under the alternative
  // positive effects
  altResPos = METRO_EM_Algorithm(
    betaeQTL, betaGWAS, Dz, D, nz, n, sigma2m, sigma2beta,
    "alternative", "positive", maxIter, tol, verbose);
  // negative effects
  altResNeg = METRO_EM_Algorithm(
    betaeQTL, betaGWAS, Dz, D, nz, n, sigma2m, sigma2beta,
    "alternative", "negative", maxIter, tol, verbose);
  Rcpp::NumericVector altlogliksPos = altResPos["logliks"];
  Rcpp::NumericVector altlogliksNeg = altResNeg["logliks"];
  altRes = (Rcpp::max(altlogliksPos) > Rcpp::max(altlogliksNeg)) 
    ? altResPos : altResNeg;

  // 3. Heritability estimation
  arma::vec mubeta = altRes["mubeta"];
  arma::mat Sigmabeta = altRes["Sigmabeta"];
  arma::vec alpham = altRes["alpha"];
  int idxStart;
  int idxEnd;
  for(int m = 0; m < M; m++)
  {
    idxStart = m * p;
    idxEnd = (m + 1) * p - 1;
    h2m(m) = 
      arma::as_scalar(
        mubeta.subvec(idxStart, idxEnd).t() * Dz.slice(m) * 
        mubeta.subvec(idxStart, idxEnd)
        ) +
      arma::trace(
        Dz.slice(m) * Sigmabeta.submat(idxStart, idxStart, idxEnd, idxEnd)
        );
  }
  h2y = 
    arma::as_scalar(
      mubeta.t() * (arma::kron(alpham * alpham.t(), D)) * mubeta 
      ) +
    arma::as_scalar(
      arma::trace(arma::kron(alpham * alpham.t(), D) * Sigmabeta)
      );

  // 4. Hypothesis testing: likelihood ratio test
  Rcpp::NumericVector nulllogliks = nullRes["logliks"];
  Rcpp::NumericVector altlogliks = altRes["logliks"];
  int df = arma::sum(h2m > hthre);
  LRTStat = -2.0 * (Rcpp::max(nulllogliks) - Rcpp::max(altlogliks));
  pvalueLRT = Rcpp::as<double>(
    pchisq(_["q"] = LRTStat, _["df"] = df, _["lower.tail"] = false)
    );
  if(df == 0)
  {
    cout << "WARNING: heritabilities from all genetic ancestries ";
    cout << "are below the threshold " << hthre << endl;
    // produce a conservative p-value
    pvalueLRT = Rcpp::as<double>(
      pchisq(_["q"] = LRTStat, _["df"] = M, _["lower.tail"] = false)
    );
  }

  // 5. weights and gene effect on outcome
  alpha = arma::sum(alpham);
  w.fill(1.0 / M);
  if(alpha != 0)
  {
    w = alpham / alpha;
  }

  // 6.SNP effects in GWAS
  for(int m = 0; m < M; m++)
  {
    beta += w(m) * mubeta.subvec(m * p, (m + 1) * p - 1);
  }

  // 7.TWAS zscores
  zscore = arma::as_scalar(beta.t() * betaGWAS * std::sqrt(n - 1)) / 
    arma::as_scalar(arma::sqrt(beta.t() * D * beta));

  auto end = chrono::steady_clock::now();
  output = Rcpp::List::create(
    _["alpha"] = alpha,
    _["weights"] = w,
    _["beta"] = beta,
    _["zscore"] = zscore,
    _["LRTStat"] = LRTStat,
    _["df"] = df,
    _["pvalueLRT"] = pvalueLRT,
    _["sigma2beta"] = altRes["sigma2beta"],
    _["sigma2m"] = altRes["sigma2m"],
    _["sigma2"] = altRes["sigma2"],
    _["h2m"] = h2m,
    _["h2y"] = h2y,
    _["p"] = p,
    _["M"] = M,
    _["nullLoglik"] = double(Rcpp::max(nulllogliks)),
    _["altLoglik"] = double(Rcpp::max(altlogliks)),
    _["elapsedTime"] = chrono::duration_cast<chrono::seconds>(end - start).count()
    );
  cout << "***** done *****" << endl;
  return output;
}


Rcpp::List METRO_EM_Algorithm(
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
  )
{
  cout << "***** Starting EM algorithm unter the " << mode << " (" << constrain << " effects)";
  cout << " *****" << endl;
  int p = betaGWAS.n_elem;
  int M = Dz.n_slices;

  // intializations
  double sigma2 = 1.0;
  arma::vec alpha(M, fill::zeros);
  if(constrain == "positive")
  {
    alpha.fill(0.01);
  }
  else if(constrain == "negative")
  {
    alpha.fill(-0.01);
  }

  // intermediate matrices and quantities
  arma::mat Sigmabeta (M * p, M * p);
  arma::mat Sigmabetainv (M * p, M * p);
  arma::vec mubeta(M * p);
  arma::mat A (M * p, M * p); // component of Sigmabetainv
  arma::cube Acompo (p, p, M); // component of A matrix
  arma::mat Ip (p, p, fill::eye); // p x p identity matrix
  int idxStart; // start index to extract sub-vectors/matrices
  int idxEnd; // end index to extract sub-vectors/matrices
  // components in updating alpha
  arma::mat alphacompo1 (M, M);
  arma::mat alphacompo2 (M, M);
  arma::vec alphacompo3 (M);
  // components in updating loglik 
  Rcpp::NumericVector logliks;
  double logdet;
  double sign;
  // components in updating lambda
  arma::vec lambda(M, fill::ones); // expanded parameter in PX-EM algorithm
  double lambdacompo1;
  double lambdacompo2;
  int iter = 0;


  while(iter < 2 || (iter < maxIter && std::abs(logliks(iter - 1) - logliks(iter - 2)) >= tol))
  {
    // E-step
    // E-step 1: update Sigma_beta
    // (1) update component of A matrix
    for(int m = 0; m < M; m++)
    {
      Acompo.slice(m) = ((nz(m) - 1.0) / sigma2m(m)) * Dz.slice(m) + 
      // Ip / sigma2beta;
      Ip / sigma2beta(m);
    }
    // (2) update A inverse matrix
    A = blockDiagonal(Acompo);

    // (3) update Sigmabeta
    Sigmabetainv = A + ((n - 1.0) / sigma2) * arma::kron(alpha * alpha.t(), D);
    Sigmabeta = arma::inv(Sigmabetainv);

    // E-step 2: update mu_beta
    mubeta = Sigmabeta * 
      (
        arma::kron((nz - 1.0) % (1.0 / sigma2m), arma::ones<arma::vec>(p)) % betaeQTL + 
        ((n - 1.0) / sigma2) * arma::kron(alpha, betaGWAS)
      );

    // Compute log-likelihood
    arma::log_det(logdet, sign, Sigmabetainv);
    logliks.push_back(
      -0.5 * arma::sum(nz % arma::log(sigma2m)) -
      0.5 * n * std::log(sigma2) -
      0.5 * p * arma::sum(arma::log(sigma2beta)) -
      0.5 * logdet - 
      0.5 * arma::sum((nz - 1.0) / sigma2m) -
      0.5 * (n - 1) / sigma2 +
      0.5 * arma::as_scalar(mubeta.t() * Sigmabetainv * mubeta)
      );

    // M-step
    // 1. update alpha
    if(mode == "alternative")
    {
      for(int i = 0; i < M; i++)
      {
        for(int j = 0; j < M; j++)
        {
          alphacompo1(i, j) = 
            arma::as_scalar(
              mubeta.subvec(i * p, (i + 1) * p - 1).t() * D * 
              mubeta.subvec(j * p, (j + 1) * p - 1)
              );
          alphacompo2(i, j) =
            arma::trace(
              D * Sigmabeta.submat(i * p, j * p, (i + 1) * p - 1, (j + 1) * p - 1)
              );
        }
        alphacompo3(i) = arma::as_scalar(
          mubeta.subvec(i * p, (i + 1) * p - 1).t() * betaGWAS);
      }
      // handle the situation when optimal parameters are on the boundary
      if(arma::cond(alphacompo1 + alphacompo2) > 1e15) break;
      alpha = arma::solve(alphacompo1 + alphacompo2, alphacompo3);
    }

    // 2. update sigma2
    sigma2 = (1.0 + 
      arma::as_scalar(mubeta.t() * arma::kron(alpha * alpha.t(), D) * mubeta) -
      2.0 * arma::as_scalar(mubeta.t() * arma::kron(alpha, betaGWAS)) +
      arma::trace(arma::kron(alpha * alpha.t(), D) * Sigmabeta)
      );

    // 3. update sigma2beta
    for(int m = 0; m < M; m++)
    {
      idxStart = m * p;
      idxEnd = (m + 1) * p - 1;
      sigma2beta(m) =
        arma::as_scalar(
          mubeta.subvec(idxStart, idxEnd).t() * mubeta.subvec(idxStart, idxEnd)
          ) +
        arma::trace(
          Sigmabeta.submat(idxStart, idxStart, idxEnd, idxEnd)
          );
      sigma2beta(m) = sigma2beta(m) / p;
    }

    // 4. update sigma2m
    for(int m = 0; m < M; m++)
    {
      idxStart = m * p;
      idxEnd = (m + 1) * p - 1;
      sigma2m(m) = (1.0 +
        arma::as_scalar(
          mubeta.subvec(idxStart, idxEnd).t() * Dz.slice(m) * 
          mubeta.subvec(idxStart, idxEnd)
          ) -
        arma::as_scalar(
          2.0 * mubeta.subvec(idxStart, idxEnd).t() * 
          betaeQTL.subvec(idxStart, idxEnd)
          ) +
        arma::trace(
          Dz.slice(m) * Sigmabeta.submat(idxStart, idxStart, idxEnd, idxEnd)
          ));
    }

    // 5. update lambda
    for(int m = 0; m < M; m++)
    {
      idxStart = m * p;
      idxEnd = (m + 1) * p - 1;
      lambdacompo1 = 
        arma::as_scalar(
          mubeta.subvec(idxStart, idxEnd).t() * 
          betaeQTL.subvec(idxStart, idxEnd)
        );
      lambdacompo2 = 
        arma::as_scalar(
          mubeta.subvec(idxStart, idxEnd).t() * Dz.slice(m) * 
          mubeta.subvec(idxStart, idxEnd)
          ) + 
        arma::trace(
          Dz.slice(m) * Sigmabeta.submat(idxStart, idxStart, idxEnd, idxEnd)
          );
      lambda(m) = lambdacompo1 / lambdacompo2;
    }

    // reduction step
    alpha = alpha / lambda;
    constrainVec(alpha, constrain);
    sigma2beta = sigma2beta % lambda % lambda;
    lambda.fill(1.0);

    iter += 1;
  }


  if(verbose)
  {
    cout << "    log-likelihood: " << Rcpp::max(logliks) << endl;
    cout << "    sigma2: " << sigma2 << endl;
    cout << "    sigma2beta: ";
    sigma2beta.t().raw_print();
    cout << "    sigma2m: ";
    sigma2m.t().raw_print();
    cout << "    alpha: ";
    alpha.t().raw_print(); 
  }


  Rcpp::List res = Rcpp::List::create(
    _["sigma2"] = sigma2,
    _["sigma2beta"] = sigma2beta,
    _["sigma2m"] = sigma2m,
    _["alpha"] = alpha,
    _["logliks"] = logliks,
    _["mubeta"] = mubeta,
    _["Sigmabeta"] = Sigmabeta
    );
  return res;
}


arma::mat blockDiagonal(
  const arma::cube X
  )
{
  int m = X.n_rows;
  int n = X.n_cols;
  int l = X.n_slices;
  arma::mat blockDiagMat (l * m, l * n, fill::zeros);
  for(int i = 0; i < l; i++)
  {
    blockDiagMat.submat(i * m, i * n, (i + 1) * m - 1, (i + 1) * n - 1) = X.slice(i);
  }
  return blockDiagMat;
}


arma::mat blockRowCol(
  const arma::cube X,
  std::string mode
  )
{
  int m = X.n_rows;
  int n = X.n_cols;
  int l = X.n_slices;
  arma::mat blockMat (m, n * l, fill::zeros);
  for(int i = 0; i < l; i++)
  {
    blockMat.submat(0, i * n, m - 1, (i + 1) * n - 1) = X.slice(i);
  }
  if(mode == "column")
  {
    return blockMat.t();
  }
  return blockMat;
}


void constrainVec(
  arma::vec &alpha,
  std::string sign
  )
{
  if(sign == "positive")
  {
    alpha.transform( [](double val) { return (val < 0.0) ? 0.0 : val; } );
  }
  else if(sign == "negative")
  {
    alpha.transform( [](double val) { return (val > 0.0) ? 0.0 : val; } );
  }
}