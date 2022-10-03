#' Check and download PLINK 1.90 executable suitable for the operating system, and return its path
#' Inspired by https://github.com/MRCIEU/genetics.binaRies
#' @importFrom utils download.file
download_plink <- function() {
  # Identify operating system
  os <- Sys.info()[["sysname"]]
  # Identify the target executable based on os
  exename <- ifelse(os == "Windows", "plink.exe", "plink")
  # Initiate a destination directory in which the executable will be stored
  dest <- file.path(system.file(package = "SlopeHunter"), "bin")
  # Create bin folder in the package directory and download specified executable into it (if not exist)
  if (!dir.exists(dest)) dir.create(dest)
  # Full path of the executable file after downloaded
  destfile <- file.path(dest, exename)

  if(!file.exists(destfile)) {
    # Get the full URL to download the executable file
    urlpref <- "https://data-science.essex.ac.uk/fileserve/ld/plink_binaries/"
    urlfull <- paste0(urlpref, os, "/", exename)
    err <- try(download.file(url = urlfull, destfile = destfile))
    if (!inherits(err, "try-error")) {
      message("PLINK 1.90 executable has been downloaded successifully.")

      # If Unix/Linux, Set mode to allow execute plink
      if(os != "Windows") {
        Chmod_command <- paste0("chmod 755 ", shQuote(destfile, type="sh"))
        system(Chmod_command)
      }

    } else {
      stop("PLINK 1.90 executable couldn't be downlaoded!")
    }

  } else {
    message("PLINK executable is already available on your local system and will be used for clumping.")
  }

  return(destfile)
}


#' clump function using local plink binary and ld reference dataset
#' This function is modified from: https://github.com/MRCIEU/ieugwasr/blob/master/R/ld_clump.R
#' @param dat Dataframe. Must have a variant name column (`rsid`) and pval column called (`pval`).
#' @param clump_kb Clumping window, default is `250`.
#' @param clump_r2 Clumping r-squared threshold, default is `0.1`.
#' @param clump_p1 Clumping sig level for index SNPs, default is `1`.
#' @param bfile Path to the bed/bim/fam LD reference (e.g. "1kg.v3/EUR" for local 1000 EUR ref. population file).
#' @importFrom data.table fwrite data.table fread
#'
ld_local <- function(dat, clump_kb=250, clump_r2=0.1, clump_p1=1, bfile) {
  # Identify operating system
  os <- Sys.info()[["sysname"]]
  shell <- ifelse(os == "Windows", "cmd", "sh")

  # create input text file for PLINK from dat
  Input <- tempfile()
  data.table::fwrite(data.table::data.table(SNP=dat$rsid, P=dat$pval), file=Input, quote=FALSE, sep = " ")

  # get path to the executable PLINK or download it and return its path
  plink_bin = download_plink()

  # make PLINK command
  PLINK_command <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --clump ", shQuote(Input, type=shell),
    " --clump-p1 ", clump_p1,
    " --clump-r2 ", clump_r2,
    " --clump-kb ", clump_kb,
    " --out ", shQuote(Input, type=shell)
  )

  # run PLINK
  system(PLINK_command)

  # Extract clumped SNPs
  Output <- data.table::fread(file = paste0(Input, ".clumped"), header = TRUE)

  # clean
  message("cleaning ...")
  unlink(paste0(Input, "*"))

  # Get pruned SNPs
  Pruned <- dat[dat$rsid %in% Output$SNP, ]

  RM <- nrow(dat) - nrow(Pruned)
  if(nrow(Pruned) > 0) {
    message("Removing ", RM, " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel ...")
  }
  return(Pruned)
}


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
#' @param ref_bfile Path to the bed/bim/fam LD reference (e.g. "1kg.v3/EUR" for local 1000 EUR ref. population file). If \code{local}=TRUE, then this should be provided.
#' @param seed Random number seed for random pruning
#'
#' @export
#' @return Data frame
#' @importFrom ieugwasr ld_clump

LD_prune <- function(dat, clump_kb=250, clump_r2=0.1, Random = TRUE, clump_p1=1, local=FALSE, ref_pop="EUR", ref_bfile, seed = 77777)
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
    if(missing(ref_bfile)) {stop("ref_bfile should be given if you set local = TRUE. You can use the MRC-IEU GWAS API for pruning by setting local = FALSE")}
    out <- ld_local(variants, clump_kb=clump_kb, clump_r2=clump_r2, clump_p1=clump_p1, bfile=ref_bfile)
  } else {
    # Pruning using API
    out <- ieugwasr::ld_clump(variants, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p1, pop = ref_pop)
  }

  Keep <- dat[dat$SNP %in% out$rsid, ]
  return(as.data.frame(Keep))
}
