#' Estimate collider bias
#' @param dat Data frame. Must have header with at least the `xbeta`, `xse`, `ybeta` and `yse` columns present.
#' @param snp_col Name of column with SNP IDs.
#' @param xbeta_col Required name of column with effects on the incidence trait.
#' @param xse_col Required name of column with standard errors of \code{xbeta}.
#' @param xp_col Name of column with p-value for \code{xbeta} (optional). If not given, It will be inferred from \code{xbeta} and \code{xse}.
#' @param ybeta_col Required name of column with unadjusted effects on the prognosis trait.
#' @param yse_col Required name of column with standard errors of \code{ybeta}.
#' @param yp_col Name of column with p-value for \code{ybeta} (optional). If not given, It will be inferred from \code{ybeta} and \code{yse}.
#' @param xp_thresh p-value threshold for SNP-incidence associations. Effects with p-values larger than \code{xp.thresh} will be excluded prior to fitting the main model-based clustering.
#' @param init_pi initial value for the weight of the mixture component that represents the cluster of SNPs affecting x only.
#' @param init_sigmaIP initial value for the covariance between x and y.
#' @param Bootstrapping Logical, if TRUE estimate the standard error of the adjustment factor using the Bootstrap method.
#' @param M Number of bootstrap samples drawn to estimate the standard error of the adjustment factor.
#' @param seed Random number seed used for drawing the bootstrap samples.
#' @param Plot Logical, if TRUE (the default), calling the function should plot the final clusters.
#' @param show_adjustments Logical indicating whether to show adjusted effects of the given SNPs in the outputs.

#' @importFrom mclust Mclust mclustBIC
#' @importFrom plotly ggplotly
#' @import dplyr
#' @import ggplot2
#' @importFrom stats dist pchisq lm sd cov na.omit pt
#'
#' @export
#'
#' @return List of the following:
#' \itemize{
#' \item est: estimated adjusted associations, their standard errors and p-values (only if \code{show_adjustments} is TRUE).
#' \item b: The estimated slope (adjustment factor).
#' \item bse: Standard error of the estimated slope.
#' \item b_CI: 95\% confidence interval of the adjustment factor.
#' \item pi: Estimated probability of the mixture component of SNPs affecting only incidence.
#' \item entropy: The entropy of the estimated clusters.
#' \item plot: Generated plot of the SlopeHunter fitted model.
#' \item Fit: a Data frame  summarising the fitted model-based clustering with the following columns:
#' \itemize{
#' \item cluster: cluster of the variants defined as follows:
#' \itemize{
#' \item `Hunted` = assigned to the cluster of SNPs affecting only incidence.
#' \item `Pleiotropic` = assigned to the cluster affecting both incidence and prognosis - i.e. variants that affect incidence and have direct effect on prognosis.
#' }
#' \item pt and p0: membership probabilities of the variants for the hunted and pleiotropic clusters respectively.
#' \item associations of variants with x and y, their standard errors and p-values.
#' }
#' \item iter: Number of the EM algorithm's iterations.
#' \item Bts.est: Details on the bootstrap estimate of the standard error of the adjustment factor, if `Bootstrapping` is TRUE.
#' }

hunt = function(dat, snp_col="SNP", xbeta_col="BETA.incidence", xse_col="SE.incidence", xp_col="Pval.incidence",
                ybeta_col="BETA.prognosis", yse_col="SE.prognosis", yp_col="Pval.prognosis",
                xp_thresh=0.001, init_pi = 0.6, init_sigmaIP = 1e-5, Bootstrapping = TRUE, M = 100, seed = 777,
                Plot=TRUE, show_adjustments = FALSE){

  # binding variable locally to the function:
  ## To avoid Notes: e.g. "hunt: no visible binding for global variable ‘xp’"; "hunt: no visible binding for global variable ‘xbeta’"
  xp <- xbeta <- ybeta <- clusters <- yse <- xse <- ybeta_adj <- yse_adj <- NULL

  all_cols <- c(snp_col, xbeta_col, xse_col, xp_col, ybeta_col, yse_col, yp_col)
  i <-  names(dat) %in% all_cols
  if (sum(i) == 0)
  {
    stop("None of the specified columns present!")
  }
  dat <- dat[, i]

  # Check if columns required for SH are present
  cols_req <- c(xbeta_col, xse_col, ybeta_col, yse_col)
  if (!all(cols_req %in% names(dat)))
  {
    stop("The following columns are not present and are required for the Slope-Hunter analysis:\n", paste(cols_req[!cols_req %in% names(dat)]), collapse="\n")
  }

  # generate p-values and SNP IDs if are not given
  cols_desired <- c(xp_col, yp_col, snp_col)
  i <- cols_desired %in% names(dat)
  if (!all(i))
  {
    message("The following column(s) is not present and will be generated:\n", paste(cols_desired[!i]))
    if(!i[1]){dat[[xp_col]] = pchisq((dat[[xbeta_col]]/dat[[xse_col]])^2, 1, lower.tail = FALSE)}
    if(!i[2]){dat[[yp_col]] = pchisq((dat[[ybeta_col]]/dat[[yse_col]])^2, 1, lower.tail = FALSE)}
    if(!i[3]){dat[[snp_col]] = paste0("snp", 1:nrow(dat))}
  }

  names(dat)[names(dat) == snp_col] <- "SNP"
  names(dat)[names(dat) == xbeta_col] <- "xbeta"
  names(dat)[names(dat) == xse_col] <- "xse"
  names(dat)[names(dat) == xp_col] <- "xp"
  names(dat)[names(dat) == ybeta_col] <- "ybeta"
  names(dat)[names(dat) == yse_col] <- "yse"
  names(dat)[names(dat) == yp_col] <- "yp"

  # filter in the variants associated with x
  gwas <- dat[, c("SNP", "xbeta", "xse", "xp", "ybeta", "yse", "yp")] %>%
    filter(xp <= xp_thresh)

  # fit slope-hunter model
  Model = shclust(gwas, pi0=init_pi, sxy1=init_sigmaIP)
  b = Model$b

  # estimate bse
  if (Bootstrapping){
    b.bts = vector("numeric", M)
    for(i in 1:M){
      print(paste("Bootstrap sample", i, "of", M, "samples ..."))   ### replace it by a progress bar ...
      set.seed(seed+i)
      Bts = sample(1:nrow(gwas), size = nrow(gwas), replace = TRUE)
      DT = gwas[Bts,]
      Model.DT = shclust(DT, pi0=init_pi, sxy1=init_sigmaIP)
      b.bts[i] = Model.DT$b
    }

    # If any of the model fits for the bootstrap samples generated NA
    if(any(is.na(b.bts))){
      b.bts = na.omit(b.bts)
      warning(paste("Only", length(b.bts), "bootstrap samples - out of", M, "- produced converged models used for estimating the standard error." ))
    }

    # calculate bse
    bse = sd(abs(b.bts))
    b_CI = c(b - 1.96*bse, b + 1.96*bse)
  } else{
    bse = Model$bse
    b_CI = Model$b_CI
  }

  # model fit
  Fit = Model$Fit

  g = ggplot2::ggplot(Fit, aes(xbeta, ybeta, label = pt, color=clusters)) +
    geom_point() +
    geom_line(aes(y=xbeta * b), color = 'black') +
    theme_bw() +
    scale_color_manual(values = c("black", "red"))

  # plot
  if (Plot) print(plotly::ggplotly(g))

  if(show_adjustments)
  {
    ##### Adjust
    est = dat %>% mutate(ybeta_adj = ybeta - (b * xbeta),
                         yse_adj = sqrt(yse^2 + (b^2 * xse^2) + (xbeta^2 * bse^2) + (xse^2 * bse^2)),
                         yp_adj = pchisq((ybeta_adj/yse_adj)^2, 1, lower.tail = FALSE))
  }

  # Print main results
  print(paste0("Estimated slope: ", b))
  print(paste0("SE of the slope: ", bse))
  print(paste0("95% CI: ", b_CI[1], ", ", b_CI[2]))

  # Return
  if(Bootstrapping & show_adjustments){
    SH = list(est=est, b=b, bse=bse, b_CI=b_CI, pi=Model$pi0, entropy=Model$entropy, plot=g, Fit=Model$Fit, iter = Model$iter, Bts.est = b.bts)
  }
  if(Bootstrapping & !show_adjustments){
    SH = list(b=b, bse=bse, b_CI=b_CI, pi=Model$pi0, entropy=Model$entropy, plot=g, Fit=Model$Fit, iter = Model$iter, Bts.est = b.bts)
  }
  if(!Bootstrapping & show_adjustments){
    SH = list(est=est, b=b, bse=bse, b_CI=b_CI, pi=Model$pi0, entropy=Model$entropy, plot=g, Fit=Model$Fit, iter = Model$iter)
  }
  if(!Bootstrapping & !show_adjustments){
    SH = list(b=b, bse=bse, b_CI=b_CI, pi=Model$pi0, entropy=Model$entropy, plot=g, Fit=Model$Fit, iter = Model$iter)
  }
  class(SH) <- "SH"
  return(SH)
}
