#' Estimate index even bias
#' @param dat Data frame. Must have header with at least the `xbeta`, `xse`, `ybeta` and `yse` columns present.
#' @param snp_col Name of column with SNP IDs.
#' @param xbeta_col Required name of column with effects on the incidence trait.
#' @param xse_col Required name of column with standard errors of \code{xbeta}.
#' @param xp_col Name of column with p-value for \code{xbeta} (optional). If not given, It will be inferred from \code{xbeta} and \code{xse}.
#' @param ybeta_col Required name of column with unadjusted effects on the prognosis trait.
#' @param yse_col Required name of column with standard errors of \code{ybeta}.
#' @param yp_col Name of column with p-value for \code{ybeta} (optional). If not given, It will be inferred from \code{ybeta} and \code{yse}.
#' @param comp.size Vector of complementary set sizes as a proportion of the noisy SNPs (neither affecting incidence nor subsequent traits).
#' The function will iteratively fit the `slopehunter` model using each of these sizes.
#' @param xp.thresh p-value threshold for SNP-incidence associations. Effects with p-values higher than th\code{xp.thresh} will be dropped from the analysis.
#' @param coef.diff A multiplier (coefficient) of the minimum standard deviation of the fitted clusters (say 'sd.min'). The product of
#' \code{coef.diff} with 'sd.min' giving a ceiling for the Euclidean distance between the two cluster means of an accepted model fit.
#' It is used for trimming the set of candidate models corresponding to each size in \code{comp.size}. The default value is 1. Larger values result in
#' larger set of possible model solutions.
#' @param seed Random number seed for the cluster-based models and weighted selections from the noisy SNPs
#' @param correct.reg.dill Logical indicating whether to correct for regression dillution using the Hedges Olkin method (the default is TRUE).
#' @param show_adjustments Logical indicating whether to show adjusted effects of the given SNPs in the outputs.
#' @param Plot Logical, if TRUE (the default), calling the function should plot the final clusters
#' @importFrom mclust Mclust mclustBIC
#' @importFrom plotly ggplotly
#' @import ggplot2
#' @importFrom stats dist pchisq lm
#'
#' @export
#'
#' @return List of the following:
#' \itemize{
#' \item Estimates: Data frame with:
#' \itemize{
#' \item status of the variants defined as follows:
#' \itemize{
#' \item `G4-Drop`: assigned to the noisy cluster (excluded from the analysis).
#' \item `G3-Drop` = assigned to the cluster affecting prognosis with no effect on incidence (excluded from the analysis).
#' \item `G1` = the target cluster.
#' \item `G2` = assigned to the cluster affecting both incidence and prognosis - i.e. variants with pleiotropic effect.
#' \item `G1-Retained` = Retained from the excluded G4 cluster and assigned to the target cluster.
#' \item `G2-Retained` = Retained from the excluded G4 cluster and assigned to the cluster affecting both incidence and prognosis traits).
#' }
#' \item membership probabilities of the vaiants based on which they are assigned to the target cluster (obviously these are calculated for only the variants included in the analysis)
#' \item adjusted effect sizes, their standard errors and p-values (only if \code{show_adjustments} is TRUE).
#' }
#' \item Iterations: Details of the iterations.
#' \item Model: the fitted model.
#' \item b.raw: The estimated slope (adjustment factor).
#' \item Sh.b: The estimated slope adjusted for regression dilution, if `correct.reg.dill` is TRUE.
#' \item bse: Standard error of the estimated slope.
#' }

slopehunter = function(dat, snp_col="SNP", xbeta_col="BETA.incidence", xse_col="SE.incidence", xp_col="Pval.incidence",
                       ybeta_col="BETA.prognosis", yse_col="SE.prognosis", yp_col="Pval.prognosis",
                       comp.size = seq(0.03, 0.10, 0.01), xp.thresh=0.10, coef.diff = 1, seed=2019,
                       correct.reg.dill = TRUE, show_adjustments = FALSE, Plot = TRUE){

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
    if(!i[1]){dat[[xp_col]] = pchisq((dat[[xbeta_col]]/dat[[xse_col]])^2, 1, lower = F)}
    if(!i[2]){dat[[yp_col]] = pchisq((dat[[ybeta_col]]/dat[[yse_col]])^2, 1, lower = F)}
    if(!i[3]){dat[[snp_col]] = paste0("snp", 1:nrow(dat))}
  }

  names(dat)[names(dat) == snp_col] <- "SNP"
  names(dat)[names(dat) == xbeta_col] <- "xbeta"
  names(dat)[names(dat) == xse_col] <- "xse"
  names(dat)[names(dat) == xp_col] <- "xp"
  names(dat)[names(dat) == ybeta_col] <- "ybeta"
  names(dat)[names(dat) == yse_col] <- "yse"
  names(dat)[names(dat) == yp_col] <- "yp"

  InData.raw <- dat[, c("SNP", "xbeta", "xse", "xp", "ybeta", "yse", "yp")]
  InData <- InData.raw

  # standardise effects
  x.Stand <- std(beta = InData$xbeta, se = InData$xse)
  xSigma <- x.Stand$Sigma
  InData$xbeta <- x.Stand$beta.std
  InData$xse   <- x.Stand$se.std
  InData$xp    <- pchisq((InData$xbeta/InData$xse)^2, 1, lower = F)

  y.Stand <- std(beta = InData$ybeta, se = InData$yse)
  ySigma <- y.Stand$Sigma
  InData$ybeta <- y.Stand$beta.std
  InData$yse   <- y.Stand$se.std
  InData$yp    <- pchisq((InData$ybeta/InData$yse)^2, 1, lower = F)

  ### Phase 1: Identify SNPs neither affect x nor y (noisy) ###
  set.seed(seed)
  Fit = mclust::Mclust(InData[c("xbeta", "ybeta")], 2, modelNames = "VVV")
  InData$class = Fit$classification

  # Euclidean distance from the origin
  InData$dist = Euclid.dist(InData$xbeta, InData$ybeta)
  KeepClass = which.max(tapply(InData$dist, InData$class, min))

  # Weighted selection from the noisy SNPs
  InData$preprocess = ifelse(InData$class == KeepClass, "Included", "Excluded")
  InData$Weights = NA
  # generate weights for noisy SNPs using the Euclidean dist calculated with higher influence for the y-dimension)))
  Excluded <- InData$preprocess == "Excluded"
  InData$Weights[Excluded] = Euclid.dist(InData$xbeta[Excluded], InData$ybeta[Excluded], ydominant = TRUE)

  # Fit the model iteratively using sizes from comp.size
  Valid0 = FALSE    # initial status of cluster solution
  Store = list()    # initial list to store details of iterations

  for (Size in comp.size){
    cat(paste0("\n Fitting Slope-Hunter by retaining ", Size*100, "% of the SNPs with neither effects on X nor Y: "))
    set.seed(seed+1)
    Sel = sample(which(Excluded), sum(Excluded)*Size, replace = FALSE, prob = InData$Weights[Excluded]) # perform weighted selections
    DT = InData
    DT$preprocess[Sel] = "Retained"
    DT$preprocess[Excluded & DT$preprocess != "Retained"] = "Dropped"

    # prepare data for model-based clustering
    DT$analysis = ifelse((DT$preprocess == "Included" & DT$xp < xp.thresh) | DT$preprocess == "Retained", "In", "Out")

    ### Phase 2: Estimating slope of the cluster of SNPs affecting x only: Clusters with different orientations for both cluster components are only allowed.
    # modelNames are formatted as: "Volume Shape Orientation", e.g. EVV = equal volume, variant shape and orientation.
    set.seed(seed+2)
    Fit2 = mclust::Mclust(DT[DT$analysis == "In", c("xbeta", "ybeta")], 2, modelNames = c("EEV", "VEV", "EVV", "VVV"))

    # Check model convergence
    VarCov = Fit2$parameters$variance$sigma
    MinVar = min(VarCov[1,1,1], VarCov[2,2,1], VarCov[1,1,2], VarCov[2,2,2])
    # valid candidate solution if Euclidean dist between cluster means <= 'coef.diff' * min cluster SDs
    Valid  = (dist(t(Fit2$parameters$mean)) <= coef.diff * sqrt(MinVar))

    # save valid solutions into the 'Store' list
    if(Valid){
      DT$Fitclass = -9
      DT$Fitclass[DT$analysis == "In"] = Fit2$classification
      Store[[paste0(Size)]] = list(DT = DT, Fit2 = Fit2)
    }

    # stop the itrations if hits an invalid fit after at least one vaild fit
    if(Valid - Valid0 == -1)
    {
      message(paste0("The Slope-Hunter model couldn't converge using ", Size*100, "% of the SNPs with neither effects on X nor Y\n
                     The models converged using ", paste0(as.numeric(names(Store))*100, collapse = ", "), "% will be used as possible solutions."))
      break
    }
    # update the solution status
    Valid0 = Valid
  }

  if(length(Store) == 0)
  {
    stop("The Slope-Hunter model couldn't converge! \nTry higher values for the 'coef.diff', but the results would need to be checked carfully! and/Or try
    different range(s) for the 'comp.size'. \nIf there still no convergence, then the Slope-Hunter might be inadequate for analysing these data")
  }

  # Identify the optimal solution among valid cluster solutions as the one with the highest precision (lowest SE) of the estimated slopes
  # calc. SE(slope) for each class
  SEs = lapply(Store, function(x)
  {
    c(summary(lm(ybeta ~ xbeta - 1, data = x$DT[x$DT$Fitclass == 1, ]))$coefficients[,c("Estimate", "Std. Error")],
      summary(lm(ybeta ~ xbeta - 1, data = x$DT[x$DT$Fitclass == 2, ]))$coefficients[,c("Estimate", "Std. Error")])
  })

  SEs <- as.data.frame(do.call(rbind, SEs))
  SEs$Targ <- apply(SEs[,c(2,4)], 1, which.min)

  Iter <- data.frame(t(apply(SEs, 1, function(x){
    if(x["Targ"] == 1){x[1:2]} else {x[3:4]}
  })))

  Iter$Targ <- SEs$Targ
  #row.names(Iter) <- names(Store)
  names(Iter) <- c("Raw.Slope", "SE", "Targ.class")
  Bst.size <- row.names(Iter)[which.min(Iter$SE)]
  G1.class <- Iter[Bst.size, "Targ.class"]
  Iter$Optimal <- ifelse(row.names(Iter) %in% Bst.size, TRUE, FALSE)

  Res = Store[[Bst.size]]$DT
  Fit = Store[[Bst.size]]$Fit2

  Res$Status <- NA
  Res$Status[Res$preprocess == "Dropped"] <- "G4-Drop"
  Res$Status[Res$preprocess == "Included" & Res$analysis == "Out"] <- "G3-Drop"
  Res$Status[Res$preprocess == "Included" & Res$analysis == "In" & Res$Fitclass == G1.class] <- "G1"
  Res$Status[Res$preprocess == "Included" & Res$analysis == "In" & Res$Fitclass != G1.class] <- "G2"
  Res$Status[Res$preprocess == "Retained" & Res$Fitclass == G1.class] <- "G1-Retained"
  Res$Status[Res$preprocess == "Retained" & Res$Fitclass != G1.class] <- "G2-Retained"

  Res$Clusters <- NA
  Res$Clusters[Res$Status == "G1" | Res$Status == "G1-Retained"] <- "G1"
  Res$Clusters[Res$Status == "G2" | Res$Status == "G2-Retained"] <- "G2"

  b.std   <- Iter[Iter$Optimal, "Raw.Slope"] # slope for standardised data
  bse.std <- Iter[Iter$Optimal, "SE"]

  b.raw <- b.std * (ySigma/xSigma)  # calc. equivelent slope on the unstandardised scale
  bse   <- bse.std * (ySigma/xSigma)

  ## extract prbabilities
  Res$Pr.G1 <- NA
  Res$Pr.G1[Res$analysis == "In"] <- round(Fit$z[, G1.class], digits = 4)

  ## produce results
  Res <- Res[, c("SNP", "xbeta", "xse", "xp", "ybeta", "yse", "yp", "Clusters", "Status", "Pr.G1", "Fitclass")]
  Res[, c("xbeta", "xse", "xp", "ybeta", "yse", "yp")] <- InData.raw[,c("xbeta", "xse", "xp", "ybeta", "yse", "yp")]  # data on the original scale

  if(correct.reg.dill){
    Sh.b <- hOlkin.adj(Data=Res, b.raw=b.raw, xbeta="xbeta", xse="xse", G1.class = G1.class)
  } else {
    Sh.b <- b.raw
  }

  Sh.b_CI <- c(Sh.b - 1.96*bse, Sh.b + 1.96*bse)

  G = ggplot2::ggplot(Res[!is.na(Res$Clusters),], aes(x = xbeta, ybeta, label=Pr.G1, color=Clusters)) + geom_point() +
    geom_line(aes(y=xbeta * Sh.b), color = 'black')
  G <- plotly::ggplotly(G)

  if(show_adjustments)
  {
    ##### Adjust
    Res$ybeta.Adj <- Res$ybeta - (Sh.b * Res$xbeta)
    Res$yse.Adj = sqrt(Res$yse^2 + Sh.b^2 * Res$xse^2)
    Res$yp.Adj = pchisq((Res$ybeta.Adj/Res$yse.Adj)^2, 1, lower = F)
  }

  Res <- Res[,-which(names(Res) %in% "Fitclass")]
  row.names(Res) <- NULL
  print(paste0("Estimated slope: ", Sh.b))
  print(paste0("SE of the slope: ", bse))
  print(paste0("95% CI: ", Sh.b_CI[1], ", ", Sh.b_CI[2]))

  if(Plot){
    print(G)
  }

  SH <- list(Estimates=Res, Iterations=Iter, Model=Fit, b.raw=b.raw, Sh.b=Sh.b, bse=bse, Sh.b_CI=Sh.b_CI, plot.clusters=G)
  class(SH) <- "SH"
  return(SH)
}
