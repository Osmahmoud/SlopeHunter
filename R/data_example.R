#' Simulated effects on quantitative incidence and prognosis traits
#'
#' A simulated dataset for 10,000 independent variables (e.g. SNPs) consisting of regression coefficients on
#' incidence and prognosis, with their standard errors. Among all the SNPs, 5% (500 variables) have effects on
#' incidence only, 5% (500 variables) on prognosis only, and 5% have correlated effects on both with a correlation
#' coeficient of '-0.5'. The estimates are obtained from linear regression in a simulated dataset of 20,000 individuals.
#'
#' @format A data frame with 10,000 rows and 5 variables:
#' \describe{
#'   \item{xbeta}{Regression coefficient on incidence}
#'   \item{xse}{Standard error of \code{xbeta}}
#'   \item{ybeta}{Regression coefficient on prognosis}
#'   \item{yse}{Standard error of \code{ybeta}}
#'   \item{yp}{P-value of the association with prognosis}
#' }
#'
#' @examples
#' # Load the \code{SlopeHunter} package
#' require(SlopeHunter)
#'
#' # Load the input data set
#' data(data_example, package = "SlopeHunter")
#' head(data_example)
#'
#' # Implement the Slope-Hunter method
#' Sh.Model <- hunt(dat = data_example, xbeta_col="xbeta", xse_col="xse",
#'                         ybeta_col="ybeta", yse_col="yse", yp_col="yp",
#'                         xp_thresh = 0.001, Bootstrapping = TRUE, show_adjustments = TRUE, seed=2021)
#'
#' # [1] "Estimated slope: -0.274120383700514"
#' # [1] "SE of the slope: 0.0229566376478153"
#' # [1] "95% CI: -0.319115393490232, -0.229125373910796"
#'
#' # Display the estimated slope (adjustment factor)
#' Sh.Model$b
#' # [1] -0.2741204
#'
#' # Extract information about cluster memberships of SNPs included in the analysis
#' Adj <- Sh.Model$Fit
#'
#' # Show the first 6 values of the unadjusted estimated effects on prognosis
#' head(data_example$ybeta)
#' # [1] -0.0092889266  0.0005575032  0.0112203795 -0.0095533069  0.0082635203  0.0026550045
#'
#'
#' # Show results of the first 6 corrected variants:
#' head(Sh.Model$est)
#'
#' #   xbeta  xse   ybeta  yse   yp    xp    SNP  ybeta_adj   yse_adj yp_adj
#' # 1 -0.007 0.007 -0.009 0.006 0.136 0.300 snp1 -0.011      0.006   0.083
#' # 2  0.014 0.007  0.000 0.006 0.928 0.042 snp2  0.004      0.006   0.492
#' # 3 -0.011 0.007  0.011 0.006 0.072 0.097 snp3  0.008      0.006   0.220
#' # 4  0.004 0.007 -0.009 0.006 0.125 0.493 snp4 -0.008      0.006   0.208
#' # 5 -0.025 0.007  0.008 0.006 0.185 0.000 snp5  0.001      0.006   0.851
#' # 6  0.013 0.007  0.002 0.006 0.670 0.054 snp6  0.006      0.006   0.329
#'
#' # Generate an interactive plot for the estimated clusters (hover on the data points to view info)
#' require(ggplot2)
#' require(plotly)
#' ggplotly(Sh.Model$plot)
"data_example"
