# SlopeHunter: An R package for collider bias correction in conditional GWAS

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5617862.svg)](https://doi.org/10.5281/zenodo.5617862)

`SlopeHunter` is an R package containing tools for correcting for collider bias in conditional Genome-Wide Association Studies (GWAS). For bug reports, feature requests, and questions on technical issues of using the software, please open an [Issue](/../../issues).

-   [Overview](#overview)
-   [Documentation](#documentation)
-   [System Requirements](#system-requirements)
-   [Installation Guide](#installation-guide)
-   [Illustrative Example](#illustrative-example)
-   [Citation](#citation)
-   [License](#license)
-   [Issues](https://github.com/Osmahmoud/SlopeHunter/issues)

# Overview

Studying genetic associations conditioned on another phenotype, e.g. a study of blood pressure conditional on body mass index, may be influenced by selection bias. An important example - of a growing interest - of these conditional studies is the study of genetic associations with prognosis (e.g. survival, subsequent events). Studies of prognosis, of necessity, can be conducted only in those who have the disease. Selection on disease status can induce associations between causes of incidence with prognosis, potentially leading to selection bias - also termed index event bias or collider bias. The `SlopeHunter` R package implements the Slope-Hunter method described in [Mahmoud et al. 2020](https://www.biorxiv.org/content/10.1101/2020.01.31.928077v1) that is developed for correcting collider bias in such conditional Genome-Wide Association Studies.

# Documenation

Detailed documentation with usage is at the [Slope-Hunter website](http://osmahmoud.com/SlopeHunter/tutorial.html).

# System Requirements

## Hardware requirements

The `SlopeHunter` R package requires only a standard computer for most applications.

## Software requirements

The `SlopeHunter` R package is supported for *Windows*, *Linux* and *macOS*. The package has been tested in R under the following systems:

-   Linux: Ubuntu 16.04 (R 4.0.4)
-   macOS: Mojave 10.14.6 (R 4.0.4)
-   Windows: 10 (R 4.0.4)

# R Dependencies

This software requires `R (>= 3.5.0)`. If you do have an older version of R installed on your machine, you may need to install the latest R version [from here](https://cloud.r-project.org/).

The `SlopeHunter` R package dependens on the following R packages:

    ggplot2 (>= 2.1.0)
    mclust
    plotly
    stats
    ieugwasr
    data.table
    dplyr

# Installation Guide:

To install the package within R, run:

```{r}
# required only once per machine!
devtools::install_github("Osmahmoud/SlopeHunter")
```

# Illustrative Example

Get started with an example analysis by running the following:

```{r}
# load the package into your R session
require(SlopeHunter)

# load data of a toy example
data(data_example, package = "SlopeHunter")

# view the first few rows of the data example
head(data_example)

# display description of the data example
help(data_example)
```

This will show you a detailed description on:

-   The instructions on how to run the software.
-   The expected outcomes.
-   How to explore your outputs.

For more details on usage of the software, run:

```{r}
help(slopehunter)
```

# Citation

> Mahmoud, O., Dudbridge, F., Davey Smith, G. et al. A robust method for collider bias correction in conditional genome-wide association studies. Nature Communications 13, 619 (2022). (doi: <https://doi.org/10.1038/s41467-022-28119-9>)
 
# License

This project is covered under the **GNU General Public License, version 3.0 (GPL-3.0)**.
