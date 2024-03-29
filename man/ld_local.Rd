% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LD_prune.R
\name{ld_local}
\alias{ld_local}
\title{clump function using local plink binary and ld reference dataset
This function is modified from: https://github.com/MRCIEU/ieugwasr/blob/master/R/ld_clump.R}
\usage{
ld_local(dat, clump_kb = 250, clump_r2 = 0.1, clump_p1 = 1, bfile)
}
\arguments{
\item{dat}{Dataframe. Must have a variant name column (\code{rsid}) and pval column called (\code{pval}).}

\item{clump_kb}{Clumping window, default is \code{250}.}

\item{clump_r2}{Clumping r-squared threshold, default is \code{0.1}.}

\item{clump_p1}{Clumping sig level for index SNPs, default is \code{1}.}

\item{bfile}{Path to the bed/bim/fam LD reference (e.g. "1kg.v3/EUR" for local 1000 EUR ref. population file).}
}
\description{
clump function using local plink binary and ld reference dataset
This function is modified from: https://github.com/MRCIEU/ieugwasr/blob/master/R/ld_clump.R
}
