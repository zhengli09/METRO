# Author: Zheng Li
# Date: 2021-11-21
# METRO with GWAS summary statistics and eQTL summary statistics
# and is extended to adjust for horizontal pleiotropic effects

#' METRO with summary level expression data and summary level GWAS data
#' and is extended to adjust for horizontal pleiotropic effects
#'
#' @param eQTLzscores A pxM numeric matrix of marginal zscores for p SNPs
#'    from M eQTL studies conducted in different ancestries.
#' @param eQTLLDs M-list of SNP-SNP correlation (LD) matrices estimated from
#'    reference panels with individuals in the same ancestries as expression
#'    studies. The order has to match across eQTL zscores, eQTL LD matrices, 
#'    and the number of individuals.
#' @param GWASzscores A numeric vector of marginal zscores for p SNPs from GWAS.
#' @param GWASLD A SNP-SNP correlation (LD) matrix estimated from a reference
#'    panel with individuals in the same ancestry as GWAS data.
#' @param ns A numeric vector. The number of individuals in each expression study.
#'    The order has to match across eQTL zscores, eQTL LD matrices, and the
#'    number of individuals.
#' @param n The number of individuals in GWAS data.
#' @param nu The shrinkage intensity parameter ranging from 0 to 1 for the LD
#'    matrix. A smaller value indicates higher intensity.
#' @param hthre A numeric scalar. Heritability threshold to determine whether
#'    cis-SNPs are predictive of gene expression.
#' @param maxIter An integer. Maximum number of iterations for PX-EM algorithm.
#' @param tol A numeric scalar. Convergence tolerance for PX-EM algorithm.
#' @param verbose A logical scalar. Whether to show verbose information.
#'
#' @return \code{METROEgger2SumStat} returns a list of estimated parameters and 
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
METROEgger2SumStat <- function(
  eQTLzscores,
  eQTLLDs,
  GWASzscores,
  GWASLD,
  ns,
  n,
  nu = 0.9,
  hthre = 2e-3,
  maxIter = 1000,
  tol = 1e-3,
  verbose = FALSE
  )
{
  if(!is.matrix(eQTLzscores)){
    eQTLzscores <- as.matrix(eQTLzscores)
  }
  M <- ncol(eQTLzscores)
  p <- length(GWASzscores)

  # error checking
  if(nrow(eQTLzscores) != p){
    stop("SNPs do not match in eQTL zscores")
  }
  if(any(sapply(eQTLLDs, nrow)!= p) | any(sapply(eQTLLDs, ncol)!= p)){
    stop("SNPs do not match in reference LD matrices for expression data")
  }
  if(any(dim(GWASLD) != p)){
    stop("SNPs do not match in reference LD matrix for GWAS data")
  }

  # preprocessing:
  if(M == 1){
    eQTLData <- eQTLzscores / sqrt(ns - 1)
  } else if(M > 1){
    eQTLData <- eQTLzscores %*% diag(1 / sqrt(ns - 1))
  }
  GWASData <- GWASzscores / sqrt(n - 1)

  # run METRO
  METRORes <- METROSummaryStatsPleio(eQTLData, GWASData, eQTLLDs, GWASLD, ns, n,
    nu = nu, hthre = hthre, maxIter = maxIter, tol = tol, verbose = verbose)
}

