#' Correct index event bias for new data
#' @param x an pbject of the class SH obtained from the \code{slopehunter} function.
#' @param dat A data.frame with harmonised effects and alleles, formatted using the \code{harmonise_effects} function.
#' @param snp_col Name of column with SNP IDs.
#' @param xbeta_col Required name of column with effects on the incidence trait.
#' @param xse_col Required name of column with standard errors of \code{xbeta}.
#' @param ybeta_col Required name of column with unadjusted effects on the prognosis trait.
#' @param yse_col Required name of column with standard errors of \code{ybeta}.
#'
#' @importFrom stats pchisq
#' @return data.frame with adjusted estimates
#' @export

SHadj = function(x, dat, snp_col="SNP", xbeta_col="BETA.incidence", xse_col="SE.incidence",
                 ybeta_col="BETA.prognosis", yse_col="SE.prognosis")
{
  if(!inherits(x, "SH")) stop("x is not of class \"SH\"!")

  # Check if columns required for adjustments are present
  cols_req <- c(xbeta_col, xse_col, ybeta_col, yse_col)
  if (!all(cols_req %in% names(dat)))
  {
    stop("The following columns are not present and are required for the Slope-Hunter adjustment:\n", paste(cols_req[!cols_req %in% names(dat)]), collapse="\n")
  }

  if(!snp_col %in% names(dat)) dat[[snp_col]] <- paste0("snp", 1:nrow(dat))

  names(dat)[names(dat) == snp_col] <- "SNP"
  names(dat)[names(dat) == xbeta_col] <- "xbeta"
  names(dat)[names(dat) == xse_col] <- "xse"
  names(dat)[names(dat) == ybeta_col] <- "ybeta"
  names(dat)[names(dat) == yse_col] <- "yse"

  Sh.b <- x$b
  bse  <- x$bse

  dat$ybeta.adj <- dat$ybeta - Sh.b * dat$xbeta
  dat$yse.adj   <- sqrt(dat$yse^2 + (Sh.b^2 * dat$xse^2) + (bse^2 * dat$xbeta^2) + (bse^2 * dat$xse^2))
  dat$yp.adj    <- pchisq((dat$ybeta.adj/dat$yse.adj)^2, 1, lower.tail = FALSE)

  return(dat)
}
