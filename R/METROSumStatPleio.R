# Author: Zheng Li
# Date: 2021-11-30
# METRO with GWAS summary statistics and is extended
# to adjust for horizontal pleiotropic effects

#' METRO with individual level expression data and summary level GWAS data
#' and is extended to adjust for horizontal pleiotropic effects
#'
#' @param eQTLGeno A list of Genotype matrices from multiple genetic ancestries.
#'    Rows are samples and columns are cis-SNPs shared across all genotype
#'    matrices in expression studies and GWAS. 
#' @param eQTLExpression A list of numeric vectors. Gene expression values from
#'    multiple genetic ancestries.
#' @param GWASzscores A numeric vector of marginal zscores from GWAS.
#' @param LDMatrix A SNP-SNP correlation (LD) matrix estimated from a reference
#'    panel with individuals in the same ancestry as GWAS data.
#' @param n The number of individuals in GWAS data.
#' @param nu The shrinkage intensity parameter ranging from 0 to 1 for the LD
#'    matrix. A smaller value indicates higher intensity.
#' @param hthre A numeric scalar. Heritability threshold to determine whether
#'    cis-SNPs are predictive of gene expression.
#' @param maxIter An integer. Maximum number of iterations for PX-EM algorithm.
#' @param tol A numeric scalar. Convergence tolerance for PX-EM algorithm.
#' @param verbose A logical scalar. Whether to show verbose information.
#'
#' @return \code{METROSumStatPleio} returns a list of estimated parameters and 
#'    statistics from METRO:
#' \item{alpha}{A numeric scalar of expression-on-outcome effect.}
#' \item{gamma}{A numeric scalar of horizontal pleiotropic effect.}
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
METROSumStatPleio <- function(
  eQTLGeno,
  eQTLExpression,
  GWASzscores,
  LDMatrix,
  n,
  nu = 0.8,
  hthre = 2e-3,
  maxIter = 1000,
  tol = 1e-3,
  verbose = FALSE
  )
{
  M <- length(eQTLGeno)
  nz <- sapply(eQTLGeno, nrow)
  p <- length(GWASzscores)

  # error checking
  if(length(eQTLExpression) != M){
    stop("Number of datasets does not match in expression studies")
  }
  if(any(sapply(eQTLExpression, length) != nz)){
    stop("Number of individuals does not match in expression studies!")
  }
  if(any(sapply(eQTLGeno, ncol) != p) | any(dim(LDMatrix) != p)){
    stop("SNPs do not match across expression studies and GWAS data")
  }

  # preprocessing:
  # 1.center and standardize each SNP
  # 2.quantile normalize expression data
  # 3.transform GWAS summary statistics as if the data are centered
  #   and standardized
  eQTLGeno <- lapply(eQTLGeno, scale, center = T, scale = T)
  eQTLExpression <- lapply(eQTLExpression, function(X){
    qqnorm(X, plot = FALSE)$x
    })
  GWASData <- GWASzscores / sqrt(n - 1)

  # generate summary statistics for eQTL data
  eQTLEst <- sapply(1:M, function(m){
    X <- eQTLExpression[[m]]
    sapply(1:p, function(j){
      snpj <- eQTLGeno[[m]][, j]
      summary(lm(X ~ snpj))$coefficients["snpj", "Estimate"]
      })
    })
  if(p == 1) eQTLEst <- t(eQTLEst)
  eQTLLD <- lapply(eQTLGeno, cor)

  # run METRO
  METRORes <- METROSummaryStatsPleio(eQTLEst, GWASData, eQTLLD, LDMatrix, nz, n,
    nu = nu, hthre = hthre, maxIter = maxIter, tol = tol, verbose = verbose)
}

