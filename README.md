# METRO (Multi-ancEstry TRanscriptOme-wide analysis)
<p align="justify"> METRO is a new computational method that leverages expression data collected from multiple genetic ancestries to enhance Transcriptome-wide association studies (TWAS). METRO incorporates expression prediction models constructed in multiple genetic ancestries through a likelihood-based inference framework, producing calibrated test statistics with substantially improved TWAS power.</p>

## Installation
Install METRO R package maintained in github through the "devtools" package.
```r
if(!require(devtools))
  install.packages(devtools)
devtools::install_github("zhengli09/METRO")
library(METRO)
```

## Usage
There are two main functions in METRO package. One is **METROIndStat** that performs multi-ancestry TWAS analysis with individual level gene expression data and individual level GWAS data. The other one is **METROSumStat** that performs multi-ancestry TWAS analysis with individual level gene expression data and GWAS summary statistics in the form of marginal z-scores. Please refer to the paper and documentations for the detailed instructions.
```r
library(METRO)
?METROIndStat
?METROSumStat
```

## Example
We exemplify the usage of METRO by analyzing the gene *ABCA1* with GEUVADIS expression data and GLGC HDL GWAS summary statistics. The two datasets are publically available. Please refer to the paper for the detailed processing of the two datasets.
```r
data(example)
METRORes <- METROSumStat(eQTLGeno, eQTLExpression, GWASzscores, LDMatrix, n,
  nu = 0.8, verbose = T)
Starting METRO...
***** info *****
  - Handling data with 315 SNPs
  - Handling data with 2 expression studies 
***** Starting EM algorithm unter the null (no effects) *****
    log-likelihood: -47373.7
    sigma2: 1
    sigma2beta: 0.00105908 1.07259e-05
    sigma2m: 0.698478 0.996663
    alpha: 0 0
***** Starting EM algorithm unter the alternative (positive effects) *****
    log-likelihood: -46824.8
    sigma2: 0.982954
    sigma2beta: 0.00105659 7.28043e-06
    sigma2m: 0.69898 0.999371
    alpha: 4.62758e-05 4.20893
***** Starting EM algorithm unter the alternative (negative effects) *****
    log-likelihood: -46827.6
    sigma2: 0.982952
    sigma2beta: 9.90301e-06 2.73241e-06
    sigma2m: 0.998661 0.999145
    alpha: -3.60259 -0.441739
***** done *****
# Key parameter estiamtes and statistics:
METRORes$alpha # 4.2
METRORes$weights # c(0, 1)
METRORes$pvalueLRT # 9.6e-241
```
