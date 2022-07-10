# METRO (Multi-ancEstry TRanscriptOme-wide analysis)
<p align="justify"> METRO is a new computational method that leverages expression data collected from multiple genetic ancestries to enhance Transcriptome-wide association studies (TWAS). METRO incorporates expression prediction models constructed in multiple genetic ancestries through a likelihood-based inference framework, producing calibrated test statistics with substantially improved TWAS power.</p>

## Installation
Install METRO R package maintained in github through the "devtools" package.
```r
if(!require(devtools))
  install.packages(devtools)
devtools::install_github("zhengli09/METRO")
library(METRO)
?METRO
```

## Usage
There are three main functions in METRO package:
1. **METROIndStat**: multi-ancestry TWAS analysis with individual level gene expression data and individual level GWAS data.
2. **METROSumStat**: multi-ancestry TWAS analysis with individual level gene expression data and GWAS summary statistics.
3. **METRO2SumStat**: multi-ancestry TWAS analysis with summary level gene expression data and GWAS summary statistics.

We have also extended METRO to adjust for SNP horizontal pleiotropy with the Egger assumption. The relevant functions are
**METROEggerIndStat**, **METROEggerSumStat**, and **METROEgger2SumStat**.

[Please refer to the paper and documentations for the detailed instructions]


## Example
We exemplify the usage of METRO by analyzing the gene *PLTP* with GEUVADIS expression data and GLGC HDL GWAS summary statistics. The two datasets are publically available. Please refer to the paper for the detailed processing of the two datasets.
```r
data(PLTP_GEUVADIS)
METRORes <- METROSumStat(eQTLGeno, eQTLExpression, GWASzscores, 
  LDMatrix, n, verbose = T)
Starting METRO...
***** info *****
  - Handling data with 124 SNPs
  - Handling data with 2 expression studies 
***** Starting EM algorithm unter the null (no effects) *****
    log-likelihood: -47352
    sigma2: 1
    sigma2beta: 0.000788926 0.0015376
    sigma2m: 0.883093 0.836214
    alpha: 0 0
***** Starting EM algorithm unter the alternative (positive effects) *****
    log-likelihood: -47287.6
    sigma2: 0.999843
    sigma2beta: 0.000782837 0.00157566
    sigma2m: 0.88352 0.835781
    alpha: 0 0
***** Starting EM algorithm unter the alternative (negative effects) *****
    log-likelihood: -47118.6
    sigma2: 0.992376
    sigma2beta: 0.000132166 0.00161972
    sigma2m: 0.993455 0.835366
    alpha: -1.67004 -0.0288259
***** done *****
# Key parameter estiamtes and statistics:
METRORes$alpha # -1.70
METRORes$weights # c(0.98, 0.02)
METRORes$pvalueLRT # 4.2e-102
```

## Citing the work
If you find the `METRO` package or any of the source code in this repository 
useful for your work, please cite:

> Zheng Li, Wei Zhao, Lulu Shang, Thomas H. Mosley, Sharon L.R. Kardia, 
> Jennifer A. Smith, Xiang Zhou# (2022). METRO: Multi-ancestry transcriptome-wide 
> association studies for powerful gene-trait association detection. American 
> Journal of Human Genetics. https://doi.org/10.1016/j.ajhg.2022.03.003.

Visit our [group website](http://www.xzlab.org) for more statistical tools on 
analyzing genetics and genomics data.