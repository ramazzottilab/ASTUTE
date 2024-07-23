# ASTUTE: Association of SomaTic mUtaTions to gene Expression profiles

| Branch | Status |
| --- | --- |
| master | [![R-CMD-check-bioc](https://github.com/ramazzottilab/ASTUTE/actions/workflows/check-bioc.yml/badge.svg?branch=master)](https://github.com/ramazzottilab/ASTUTE/actions/workflows/check-bioc.yml) |
| development | [![R-CMD-check-bioc](https://github.com/ramazzottilab/ASTUTE/actions/workflows/check-bioc.yml/badge.svg?branch=development)](https://github.com/ramazzottilab/ASTUTE/actions/workflows/check-bioc.yml) |

ASTUTE is an R package designed to integrate cancer genomic and transcriptomic data in order to perform genotype-phenotype mapping. It leverages regularized regression with LASSO penalty to uncover associations between somatic mutations and gene expression profiles.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Contact](#contact)

## Features

- **Genotype-Phenotype Mapping**: Perform sophisticated analyses to uncover the relationship between somatic mutations and gene expression profiles.
- **Regularized Regression**: Employ LASSO penalty to prevent overfitting and select significant features.
- **Bootstrapping**: Calculate p-values and confidence intervals for expression changes.
- **Patient Stratification**: Identify and utilize expression signatures for prognosis and therapeutic decision-making.

## Installation

You can install the development version of ASTUTE from GitHub using the `devtools` package:

```r
# Install devtools if you haven't already
install.packages("devtools")
library("devtools")

# Install ASTUTE from GitHub
devtools::install_github("ramazzottilab/ASTUTE")
```

## Contact

Please feel free to contact us if you have any questions regarding our tool at daniele.ramazzotti@unimib.it.
