#' Calculate Euclidean distance from origin
#' @param x Vector of x-axis values
#' @param y Vector of y-axis values
#' @param ydominant Logical, if TRUE, larger influence will be given to the y-axis values
#' @return Vector of distances
Euclid.dist = function(x, y, ydominant = FALSE){
  if(ydominant){
    xw = max(abs(x))    # derive scale for x-values
    yw = max(abs(y))    # derive scale for y-values
    D = (x/xw)^2 + abs(y/yw)
    D / max(D)
  }
  else{
    sqrt(x^2 + y^2)
  }
}

#' @importFrom stats var
hOlkin.adj <- function(Data, b.raw, xbeta, xse, G1.class)
{
  Dat <- subset(Data, Fitclass == G1.class)
  Vx  <- var(Dat[[xbeta]])
  cVx <- Vx - (mean(Dat[[xse]])^2)
  b.raw * Vx / cVx
}

#' @importFrom stats sd
std <- function(beta, se)
{
  Mu <- mean(beta, na.rm = TRUE)
  Sigma <- sd(beta, na.rm = TRUE)
  beta.std <- (beta-Mu)/Sigma
  se.std <- se/Sigma
  return(list(beta.std=beta.std, se.std=se.std, Mu=Mu, Sigma=Sigma))
}

#' Plotting model for Slope-Hunter clustering
#'
#' @param x Output from \code{slopehunter}.
#' @param what A string specifying the type of graph requested. Available choices are:
#' "clusters": showing clusters. The plot can display membership probabilities of each variable (e.g. SNP) to the target cluster (G1) by hovering over the points.
#' "classification": A plot showing point assigned to each cluster (class).
#' "uncertainty": A plot of classification uncertainty.
#' "density": A plot of estimated density.
#' @param xlab Optional label for the x-axis in case of "classification", "uncertainty", or "density" plots.
#' @param ylab Optional label for the y-axis in case of "classification", "uncertainty", or "density" plots.
#' @param addEllipses A logical indicating whether or not to add ellipses with axes corresponding to the within-cluster covariances in case of "classification" or "uncertainty" plots.
#' @param main A logical or NULL indicating whether or not to add a title to the plot identifying the type of plot drawn in case of "classification", "uncertainty", or "density" plots.
#' @param ... Other graphics parameters.
#' @export
#' @importFrom utils menu
plot.SH <- function(x, what= c("clusters", "classification", "uncertainty", "density"),
                    xlab = NULL, ylab = NULL, addEllipses = TRUE, main = FALSE, ...)
{
  object <- x
  if(!inherits(object, "SH"))
    stop("x is not of class \"SH\"!")

  if(interactive() & length(what) > 1){
    title <- "Slope-Hunter fitting plots:"
    # present menu waiting user choice
    choice <- menu(what, graphics = FALSE, title = title)
    while(choice != 0){
      if(what[choice] == "clusters")       print(object$plot.clusters)
      if(what[choice] == "classification") mclust::plot.Mclust(object$Model, what = "classification", xlab=xlab, ylab=ylab, addEllipses=addEllipses, main=main)
      if(what[choice] == "uncertainty")    mclust::plot.Mclust(object$Model, what = "uncertainty", xlab=xlab, ylab=ylab, addEllipses=addEllipses, main=main)
      if(what[choice] == "density")        mclust::plot.Mclust(object$Model, what = "density", xlab=xlab, ylab=ylab, main=main)
      # re-present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
    }
  } else {
    if(any(what == "clusters"))       print(object$plot.clusters)
    if(any(what == "classification")) mclust::plot.Mclust(object$Model, what = "classification", xlab=xlab, ylab=ylab, addEllipses=addEllipses, main=main)
    if(any(what == "uncertainty"))    mclust::plot.Mclust(object$Model, what = "uncertainty", xlab=xlab, ylab=ylab, addEllipses=addEllipses, main=main)
    if(any(what == "density"))        mclust::plot.Mclust(object$Model, what = "density", xlab=xlab, ylab=ylab, main=main)
  }
  invisible()
}

#' Implement the EM algorithm for the SlopeHunter model-based clustering
#'
#' @param gwas a data frame with columns: xbeta; xse; ybeta; yse.
#' @param pi0 initial value for the weight of the mixture component that represents the cluster of SNPs affecting x only.
#' @param sxy1 initial value for the covariance between x and y.
#' @return EM fit for SlopeHunter estimator
#' @importFrom mclust dmvnorm
#' @import dplyr
#' @importFrom stats sd cov na.omit pt
shclust <- function(gwas, pi0, sxy1){
  sx0 = sx1 = gwas %>% summarise(sd(xbeta)) %>% pull
  sy0 = sy1 = gwas %>% summarise(sd(ybeta)) %>% pull
  dir0 = gwas %>% summarise(cov(xbeta, ybeta)) %>% pull %>% sign()
  if (dir0==0) stop("All associations with at least either x or y are constant")

  # convergence criterion
  loglkl_ck = 0

  ### EM algorithm
  for(iter in 1:50000){
    #### The E step:
    # covariance matrix for the target component (f0)
    sxy0 = sx0*sy0*dir0*0.95       # the x & y perfectly correlated under 1st component  #===========
    sigma0 = matrix(c(sx0^2,sxy0,sxy0,sy0^2), 2, 2)

    # 1st component
    f0 = gwas %>%
      dplyr::select(xbeta, ybeta) %>%
      mclust::dmvnorm(mean=c(0,0), sigma=sigma0)
    f0[f0<1e-300] = 1e-300

    # covariance matrix for the component (f1)
    sigma1 = matrix(c(sx1^2,sxy1,sxy1,sy1^2), 2, 2)

    # 2nd component
    f1 = gwas %>%
      dplyr::select(xbeta, ybeta) %>%
      mclust::dmvnorm(mean=c(0,0), sigma=sigma1)
    f1[f1<1e-300] = 1e-300

    # loglik of the mixture model: pi0 * f0 + (1-p0) * f1
    loglkl = sum(log(pi0*f0+(1-pi0)*f1))

    ## proportional contribution of density of f0 (for every point) to the total mixture
    pt = pi0*f0/(pi0*f0+(1-pi0)*f1)
    pt[pt>0.9999999] = 0.9999999

    #### The M step:
    # update pi0
    pi0 = mean(pt)
    if (pi0<0.0001) pi0 = 0.0001
    if (pi0>0.9999) pi0 = 0.9999

    # update sx0 & sy0
    sx0 = gwas %>% summarise(sqrt(sum(pt*(xbeta^2))/sum(pt))) %>% pull
    sy0 = gwas %>% summarise(sqrt(sum(pt*(ybeta^2))/sum(pt))) %>% pull
    dir0 = gwas %>% summarise(sum(pt*xbeta*ybeta)/sum(pt)) %>% pull %>% sign()
    if (dir0==0) dir0=sample(c(1,-1), 1)   # avoid slope = 0 (horizontal line)

    # update sx1, sy1 & sxy1
    sx1 = gwas %>% summarise(sqrt(sum((1-pt)*(xbeta^2))/(length(xbeta)-sum(pt)))) %>% pull
    sy1 = gwas %>% summarise(sqrt(sum((1-pt)*(ybeta^2))/(length(ybeta)-sum(pt)))) %>% pull
    sxy1 = gwas %>% summarise(sum((1-pt)*xbeta*ybeta)/(length(xbeta)-sum(pt))) %>% pull
    if (abs(sxy1) > 0.75*sx1*sy1)  sxy1 = sign(sxy1)*0.75*sx1*sy1     #===========

    ## Check convergence
    if (iter%%10==0){
      if ((loglkl - loglkl_ck)/loglkl < 1e-10){
        break
      } else {
        loglkl_ck = loglkl
      }
    }
  }

  # Diagnosis
  if (iter == 50000) warning("The algorithm may not have converged.\n")

  ### Results
  Fit = gwas %>% mutate(pt = pt) %>%
    mutate(po = 1-pt, clusters = factor(ifelse(pt >= 0.5, "Hunted", "Pleiotropic")))

  # slope of the eigenvector
  b = dir0*sy0/sx0
  bse = 0
  b_CI = c(b - 1.96*bse, b + 1.96*bse)
  entropy = Fit %>% filter(clusters == "Hunted") %>% summarise(mean(pt)) %>% pull

  return(list(b=b, bse=bse, b_CI=b_CI, iter=iter, pi0=pi0, entropy=entropy, Fit=Fit))
}
