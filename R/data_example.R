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
#' Sh.Model <- slopehunter(dat = data_example, xbeta_col="xbeta", xse_col="xse",
#'                         ybeta_col="ybeta", yse_col="yse", yp_col="yp",
#'                         comp.size = seq(0.03, 0.10, 0.01), xp.thresh = 0.1, coef.diff = 1,
#'                         correct.reg.dill = TRUE, show_adjustments = TRUE, seed=2019)
#'
#' # [1] "Estimated slope: -0.262779432124327"
#' # [1] "SE of the slope: 0.0147636543987382"
#' # [1] "95% CI: -0.291716194745854, -0.2338426695028"
#'
#' # Display the estimated slope (adjustment factor)
#' Sh.Model$Sh.b
#' # [1] -0.2627794
#'
#' # Extract resulted data with adjusted estimates and information of SNPs included in the analysis
#' Adj <- Sh.Model$Estimates
#'
#' # Show the first 6 values of the unadjusted estimated effects on prognosis
#' head(Adj$ybeta)
#' # [1] -0.0092889266  0.0005575032  0.0112203795 -0.0095533069  0.0082635203  0.0026550045
#'
#' # Show the first 6 values of the adjusted estimated effects on prognosis
#' head(Adj$ybeta.Adj)
#' # [1] -0.011214090  0.004324046  0.008137912 -0.008280179  0.001521822  0.006230996
#'
#' # Generate an interactive plot for the estimated clusters (hover on the data points to view info)
#' plot(Sh.Model, what = "clusters")
#'
#' # Generate plot for classification uncertainty
#' plot(Sh.Model, what = "uncertainty")
"data_example"
