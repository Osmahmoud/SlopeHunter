#' Perform LD pruning on SNP data
#'
#' Uses PLINK clumping method ('--clump' command), where a greedy search algorithm is implemented to randomly select a variant
#' (or the variant with the lowest p-value, if a user wish to), referred to as the index SNP, and remove all variants within a
#' certain kb distance in linkage disequilibrium with the index SNP, based on an r-squared threshold from the 1000 Genomes
#' reference panel phase 3 data. Then repeats until no variants are left.
#'
#' @md
#' @param dat Output from \code{harmonise_effects}. Must have a SNP name column (SNP).
#' @param clump_kb Clumping window, default is `250`.
#' @param clump_r2 Clumping r-squared threshold, default is `0.1`.
#' @param Random Logical, if `TRUE` (the default), SNPs will be randomly pruned. Otherwise, based on p-values.
#' @param clump_p1 Clumping sig level for index SNPs, default is `1`.
#' @param local Logical, if `FALSE` (the default), the MRC-IEU API 'http://gwas-api.mrcieu.ac.uk/'
#' will be used for clumping. Otherwise, your local machine will be used for clumping given that you
#' provide a bed/bim/fam LD reference dataset.
#' @param ref_pop Super-population to use as reference panel at the API (when \code{local} is FALSE). Default = `"EUR"`.
#' @param seed Random number seed for random pruning
#'
#' @export
#' @return Data frame
#' @importFrom ieugwasr ld_clump

LD_prune <- function(dat, clump_kb=250, clump_r2=0.1, Random = TRUE, clump_p1=1, local=FALSE, ref_pop="EUR", seed = 77777)
{

  if(!is.data.frame(dat)) stop("Expecting data.frame or data.table returned from 'harmonise_effects'")
  # Not Random but P-val not present!
  if(!Random & !"PVAL.incidence" %in% names(dat)) {
    stop("Pruning based on P-values, but the 'PVAL.incidence' column is not present!\n Set 'Random = TRUE' to prune randomly!")
  }

  # Not Random Pruning
  if(!Random & "PVAL.incidence" %in% names(dat)) {
    message("Pruning based on P-values ...\n Using the 'PVAL.incidence' column!")
    pval_col <- "PVAL.incidence"
  }

  # Random Pruning
  if(Random) {
    message("Pruning for LD randomly ...")
    if(!is.numeric(seed)) stop("seed should be only numeric")
    l <- length(dat$SNP)
    set.seed(seed)
    dat$RanGen <- sample(seq(l), size = l, replace = FALSE)/l
    pval_col <- "RanGen"
  }

  variants <- data.frame(rsid=dat$SNP, pval=dat[[pval_col]], id="is performed for")
  message(paste0("Information for ", nrow(variants), " SNPs have been supplied.\n"))

  if(local) {
    # Pruning using local machine
    message("pruning on the local machine using PLINK")
    out <- NA  # TEMP value - to be replaced by correct coding
  } else {
    # Pruning using API
    out <- ieugwasr::ld_clump(variants, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p1, pop = ref_pop)
  }

  Keep <- dat$SNP %in% out$rsid
  return(dat[Keep, ])
}
