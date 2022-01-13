\name{METRO-package}
\alias{METRO-package}
\alias{METRO}
\docType{package}
\title{
  METRO:Multi-ancEstry TRanscriptOme-wide analysis
}
\description{
  METRO is a new computational method that leverages expression data collected 
  from multiple genetic ancestries to enhance Transcriptome-wide association 
  studies (TWAS). METRO incorporates expression prediction models constructed 
  in multiple genetic ancestries through a likelihood-based inference framework, 
  producing calibrated test statistics with substantially improved TWAS power.
}
\details{
  Refer to:\cr
  (1) "METROIndStat" that performs multi-ancestry TWAS analysis with individual level 
      gene expression data and individual level GWAS data. Additional covariates can
      be included (e.g genetic PCs) for variable adjustment.\cr
  (2) "METROSumStat" that performs multi-ancestry TWAS analysis with individual level 
      gene expression data and GWAS summary statistics.\cr
  (3) "METRO2SumStat" that performs multi-ancestry TWAS analysis with summary level 
      gene expression data and GWAS summary statistics.\cr
  (4) "METROEggerIndStat" that performs multi-ancestry TWAS analysis with individual 
      level gene expression data and individual level GWAS data and controls for
      horizontal pleiotropic effects.\cr
  (5) "METROEggerSumStat" that performs multi-ancestry TWAS analysis with individual 
      level gene expression data and GWAS summary statistics and controls for
      horizontal pleiotropic effects.\cr
  (6) "METROEgger2SumStat" that performs multi-ancestry TWAS analysis with summary 
      level gene expression data and GWAS summary statistics and controls for
      horizontal pleiotropic effects.
}
\author{
Zheng Li, Xiang Zhou.

Maintainer: Zheng Li <zlisph@umich.edu>
}
\references{
  METRO: Multi-ancestry transcriptome-wide association studies for
  powerful gene-trait association detection
}
\keyword{ package }
\examples{
  \dontrun{
    data(PLTP_GEUVADIS)
    METRORes <- METROSumStat(eQTLGeno, eQTLExpression, GWASzscores, 
      LDMatrix, n, verbose = T)
    METRORes$alpha
    METRORes$weights
    METRORes$pvalueLRT
  }
}
