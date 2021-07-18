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
#' #    xbeta         xse         ybeta         yse         yp        xp          SNP    ybeta_adj     yse_adj     yp_adj
#' # 1 -0.007326158 0.007071232 -0.0092889266 0.006237317 0.13643722 0.3001782664 snp1 -0.011297176 0.006535750 0.08389502
#' # 2  0.014333475 0.007070695  0.0005575032 0.006238135 0.92878862 0.0426454048 snp2  0.004486601 0.006542603 0.49286972
#' # 3 -0.011730245 0.007070935  0.0112203795 0.006237420 0.07205249 0.0971282073 snp3  0.008004880 0.006539207 0.22090085
#' # 4  0.004844855 0.007071338 -0.0095533069 0.006237203 0.12562087 0.4932556974 snp4 -0.008225233 0.006534432 0.20811973
#' # 5 -0.025655351 0.007069094  0.0082635203 0.006239275 0.18537346 0.0002842704 snp5  0.001230866 0.006561766 0.85120477
#' # 6  0.013608340 0.007070767  0.0026550045 0.006238045 0.67039309 0.0542804363 snp6  0.006385328 0.006541707 0.32901732
#'
#' # Generate an interactive plot for the estimated clusters (hover on the data points to view info)
#' require(ggplot2)
#' require(plotly)
#' ggplotly(Sh.Model$plot)
"data_example"
