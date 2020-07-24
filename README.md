# SlopeHunter: An R package for collider bias correction in GWAS of subsequent traits

This R package implements the Slope-Hunter method described in [Mahmoud et al. 2020](https://www.biorxiv.org/content/10.1101/2020.01.31.928077v1). Tutorials on the software usage can be viewed from the [Slope-Hunter website](http://osmahmoud.com/SlopeHunter/).

For bug reports, feature requests, and questions on technical issues of using the software, please open an [Issue](/../../issues).

## Installation instructions:

To install the package within R, run:

```{r}
# required only once per machine!
devtools::install_github("Osmahmoud/SlopeHunter")
```
This software requires `R (>= 3.5.0)`. If you do have an older version of R installed on your machine, you may need to install the latest R version [from here](https://cloud.r-project.org/).


## Illustative example

Get started with an example analysis by running the following:

```{r}
# load the package into your R session
library(SlopeHunter)

# display description of a toy example
help(data_example)
```

For more details on usage of the software, run:
```{r}
help(slopehunter)
```

## Citation

> Mahmoud O. et al, Slope-Hunter: A robust method for index-event bias correction in genome-wide association studies of subsequent traits.  Submitted for publication.
