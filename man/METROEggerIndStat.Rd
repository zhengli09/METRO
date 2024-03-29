% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/METROEggerIndStat.R
\name{METROEggerIndStat}
\alias{METROEggerIndStat}
\title{METRO with individual level expression data and GWAS data
and is extended to adjust for horizontal pleiotropic effects}
\usage{
METROEggerIndStat(
  eQTLGeno,
  eQTLExpression,
  GWASGeno,
  GWASPheno,
  hthre = 0.002,
  maxIter = 1000,
  tol = 0.001,
  verbose = FALSE
)
}
\arguments{
\item{eQTLGeno}{A list of Genotype matrices from multiple genetic ancestries.
Rows are samples and columns are cis-SNPs shared across all genotype
matrices in expression studies and GWAS.}

\item{eQTLExpression}{A list of numeric vectors. Gene expression values from
multiple genetic ancestries.}

\item{GWASGeno}{A numeric matrix of genotype data in GWAS. Rows are samples
and columns are the same set of cis-SNPs as in the expression studies.}

\item{GWASPheno}{A numeric vector of GWAS outcome variable.}

\item{hthre}{A numeric scalar. Heritability threshold to determine whether
cis-SNPs are predictive of gene expression.}

\item{maxIter}{An integer. Maximum number of iterations for PX-EM algorithm.}

\item{tol}{A numeric scalar. Convergence tolerance for PX-EM algorithm.}

\item{verbose}{A logical scalar. Whether to show verbose information.}
}
\value{
\code{METROEggerIndStat} returns a list of estimated parameters and 
   statistics from METRO:
\item{alpha}{A numeric scalar of expression-on-outcome effect.}
\item{gamma}{A numeric scalar of horizontal pleiotropic effect.}
\item{weights}{A numeric vector of weights from multiple genetic ancestries.}
\item{beta}{A numeric vector of genetic effects on gene expression among GWAS
   individuals.}
\item{LRTStat}{Likelihood ratio test statistic.}
\item{df}{Degree of freedoms of the reference chi2 distribution for LRTStat.}
\item{pvalueLRT}{P-value of the likelihood ratio test.}
\item{sigma2beta}{A numeric vector of variance estimates of cis-SNP effects
   in different genetic ancestries.}
\item{sigma2m}{A numeric vector of residual variance estimates in gene 
   expression models.}
\item{sigma2}{The residual variance estiamte in the GWAS-expression 
   association model.}
\item{h2m}{A numeric vector of cis-SNP heritabilities of gene expression in
   different genetic ancestries.}
\item{h2y}{The proportion of variance in GWAS outcome explained by gene
   expression.}
\item{p}{The number of cis-SNPs.}
\item{M}{The number of genetic ancesties.}
\item{nullLoglik}{The log-likelihood under the null.}
\item{altLoglik}{The log-likelihood under the alternative.}
\item{elapsedTime}{Time (seconds) used to run METRO.}
}
\description{
METRO with individual level expression data and GWAS data
and is extended to adjust for horizontal pleiotropic effects
}
