# Author: Zheng Li
# Date: 2021-06-22
# METRO with individual statistics

#' METRO with individual level expression data and GWAS data
#'
#' @param eQTLGeno A list of Genotype matrices from multiple genetic ancestries.
#'    Rows are samples and columns are cis-SNPs shared across all genotype
#'    matrices in expression studies and GWAS. 
#' @param eQTLExpression A list of numeric vectors. Gene expression values from
#'    multiple genetic ancestries.
#' @param GWASGeno A numeric matrix of genotype data in GWAS. Rows are samples
#'    and columns are the same set of cis-SNPs as in the expression studies.
#' @param GWASPheno A numeric vector of GWAS outcome variable.
#' @param covars A numeric matrix of covariates in the GWAS data. 
#'    Rows are samples and columns are different covariates. 
#' @param hthre A numeric scalar. Heritability threshold to determine whether
#'    cis-SNPs are predictive of gene expression.
#' @param maxIter An integer. Maximum number of iterations for PX-EM algorithm.
#' @param tol A numeric scalar. Convergence tolerance for PX-EM algorithm.
#' @param verbose A logical scalar. Whether to show verbose information.
#'
#' @return \code{METROIndStat} returns a list of estimated parameters and 
#'    statistics from METRO:
#' \item{alpha}{A numeric scalar of expression-on-outcome effect.}
#' \item{delta}{A numeric vector of effect size estimates for covariates if 
#'    covars are provided.}
#' \item{weights}{A numeric vector of weights from multiple genetic ancestries.}
#' \item{beta}{A numeric vector of genetic effects on gene expression among GWAS
#'    individuals.}
#' \item{LRTStat}{Likelihood ratio test statistic.}
#' \item{df}{Degree of freedoms of the reference chi2 distribution for LRTStat.}
#' \item{pvalueLRT}{P-value of the likelihood ratio test.}
#' \item{sigma2beta}{A numeric vector of variance estimates of cis-SNP effects
#'    in different genetic ancestries.}
#' \item{sigma2m}{A numeric vector of residual variance estimates in gene 
#'    expression models.}
#' \item{sigma2}{The residual variance estiamte in the GWAS-expression 
#'    association model.}
#' \item{h2m}{A numeric vector of cis-SNP heritabilities of gene expression in
#'    different genetic ancestries.}
#' \item{h2y}{The proportion of variance in GWAS outcome explained by gene
#'    expression.}
#' \item{p}{The number of cis-SNPs.}
#' \item{M}{The number of genetic ancesties.}
#' \item{nullLoglik}{The log-likelihood under the null.}
#' \item{altLoglik}{The log-likelihood under the alternative.}
#' \item{elapsedTime}{Time (seconds) used to run METRO.}
#'
#' @export
#'
METROIndStat <- function(
  eQTLGeno,
  eQTLExpression,
  GWASGeno,
  GWASPheno,
  covars = NULL,
  hthre = 2e-3,
  maxIter = 1000,
  tol = 1e-3,
  verbose = FALSE
  )
{
  M <- length(eQTLGeno)
  nz <- sapply(eQTLGeno, nrow)
  n <- nrow(GWASGeno)
  p <- ncol(GWASGeno)

  # error checking
  if(p == 0){
    stop("Number of SNPs is zero")
  }
  if(length(eQTLExpression) != M){
    stop("Number of datasets does not match in expression studies")
  }
  if(any(sapply(eQTLExpression, length) != nz)){
    stop("Number of individuals does not match in expression studies")
  }
  if(length(GWASPheno) != n){
    stop("Number of individuals does not match in GWAS data")
  }
  if(any(sapply(eQTLGeno, ncol) != p)){
    stop("SNPs do not match across expression studies and GWAS data")
  }

  # preprocessing:
  # 1.center and standardize each SNPs
  # 2.quantile normalize expression data
  # 3.center and standardize GWAS phenotype
  eQTLGeno <- lapply(eQTLGeno, scale, center = T, scale = T)
  GWASGeno <- scale(GWASGeno, center = T, scale = T)
  eQTLExpression <- lapply(eQTLExpression, function(X){
    qqnorm(X, plot = FALSE)$x
    })
  GWASPheno <- scale(GWASPheno)

  # generate summary statistics
  # 1.eQTL:
  eQTLEst <- sapply(1:M, function(m){
    X <- eQTLExpression[[m]]
    sapply(1:p, function(j){
      snpj <- eQTLGeno[[m]][, j]
      summary(lm(X ~ snpj))$coefficients["snpj", "Estimate"]
      })
    })
  # 2.GWAS
  GWASEst <- sapply(1:p, function(j){
    snpj <- GWASGeno[, j]
    summary(lm(GWASPheno ~ snpj))$coefficients["snpj", "Estimate"]
    })
  if(p == 1) eQTLEst <- t(eQTLEst)
  # 3.LD matrices
  eQTLLD <- lapply(eQTLGeno, cor)
  GWASLD <- cor(GWASGeno)

  # Handle covariates data
  if(!is.null(covars)){
    if(!is.matrix(covars)){
      covars <- as.matrix(covars)
    }
    if(nrow(covars) != n){
      stop("Number of individuals does not match in covariates data")
    }
    covars <- scale(covars, center = T, scale = F)
    covGCovars <- cov(GWASGeno, covars)
    covYCovars <- cov(covars, GWASPheno)
    covCovars <- cov(covars)
    METRORes <- METROCovars(eQTLEst, GWASEst, eQTLLD, GWASLD, nz, n,
      covGCovars, covYCovars, covCovars, nu = 1, hthre = hthre,
      maxIter = maxIter, tol = tol, verbose = verbose)
  } else{
    # run METRO
    METRORes <- METROSummaryStats(eQTLEst, GWASEst, eQTLLD, GWASLD, nz, n,
      nu = 1, hthre = hthre, maxIter = maxIter, tol = tol, verbose = verbose)
  }
  METRORes
}