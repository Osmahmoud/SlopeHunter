#' Format input data
#'
#' Reads in and format input data. It checks and organises columns for use with Slope-Hunter analyses.
#' Infers p-values when possible from beta and se.
#'
#' @md
#' @param dat Data frame. Must have header with at least the `SNP`, `beta`, `se` and `EA` columns present.
#' @param type Is this the incidence or the prognosis data that is being read in? The default is `"incidence"`.
#' @param snps SNPs to extract. If NULL, then it keeps all. The default is `NULL`.
#' @param snp_col Required name of column with SNP rs IDs. The default is `"SNP"`.
#' @param beta_col Required name of column with effect sizes. The default is `"BETA"`.
#' @param se_col Required name of column with standard errors. The default is `"SE"`.
#' @param pval_col Name of column with p-value (optional). The default is `"PVAL"`. It will be Inferred when possible from beta and se.
#' @param eaf_col Name of column with effect allele frequency (optional). The default is `"EAF"`.
#' @param effect_allele_col Required for harmonisation. Name of column with effect allele. Must be "A", "C", "T" or "G". The default is `"EA"`.
#' @param other_allele_col Required for harmonisation. Name of column with non-effect allele. Must be "A", "C", "T" or "G". The default is `"OA"`.
#' @param gene_col Optional column for gene name. The default is `"GENE"`.
#' @param chr_col Optional column for chromosome number. The default is `"CHR"`.
#' @param pos_col Optional column for SNP position. The default is `"POS"`.
#' @param min_pval Minimum allowed p-value. The default is `1e-200`.
#' @param log_pval The p-value is -log10(P). The default is `FALSE`.
#'
#' @importFrom data.table data.table setnames setcolorder
#' @importFrom stats pnorm
#' @export
#' @return data frame

format_data = function (dat, type = "incidence", snps = NULL,
                        snp_col="SNP", beta_col="BETA", se_col="SE",
                        pval_col="PVAL", eaf_col="EAF", effect_allele_col="EA",
                        other_allele_col="OA", gene_col="GENE", chr_col = "CHR", pos_col="POS",
                        min_pval=1e-200, log_pval=FALSE)
{
  all_cols <- c(snp_col, beta_col, se_col, pval_col, eaf_col, effect_allele_col, other_allele_col, gene_col, chr_col, pos_col)

  # check given names ####
  i <- names(dat) %in% all_cols
  ## Error if non of the given names are present in the data
  if (sum(i) == 0)
  {
    stop("None of the specified columns present!")
  }
  ## data.table (FAQ 1.5): select the columns whose names are given by the user (stored in i)
  Sel = names(dat)[i]
  dat <- dat[, Sel, with = FALSE]
  ## Check if columns required for SH are present
  SH_cols_req = c(snp_col, beta_col, se_col, effect_allele_col, other_allele_col)
  ## Error if not ALL of the essential columns for SH present
  if (!all(SH_cols_req %in% names(dat))){
    stop("The following columns are not present and are required for the Slope-Hunter analysis\n", paste(SH_cols_req[!SH_cols_req %in% names(dat)]), collapse="\n")
  }

  # check SNP ####
  ## Format SNP IDs to the standard: lower-case letters, remove spaces & exclude NAs ####
  data.table::setnames(dat, snp_col, "SNP")
  ## Extract only a list of SNPs if <snps> is given
  if (!is.null(snps)) {
    dat <- dat[SNP %in% snps]
  }
  dat <- dat[, SNP:= tolower(SNP)]
  dat <- dat[, SNP:= gsub("[[:space:]]", "", SNP)]
  ## Remove NAs
  if (any(is.na(dat$SNP))) {
    warning("Excluding SNPs with missing rsID ...")
    dat <- dat[!is.na(SNP)]
  }
  ## Remove duplicated SNPs
  dup <- duplicated(dat$SNP)
  if (any(dup)) {
    warning(sum(dup), " duplicated SNP rsIDs present: Just keeping the first instance ...")
    dat <- dat[!dup]
  }


  # Check beta ####
  data.table::setnames(dat, beta_col, "BETA.outcome")
  ## Coercing to numeric if not
  if(!is.numeric(dat$BETA.outcome)) {
    warning("beta column is not numeric. Coercing...")
    dat <- dat[, BETA.outcome:= as.numeric(BETA.outcome)]
  }
  ## replace infinite values (NA, NaN, Inf or -Inf), if any, with NA
  dat <- dat[, BETA.outcome:= ifelse(is.finite(BETA.outcome), BETA.outcome, NA)]


  # Check se ####
  data.table::setnames(dat, se_col, "SE.outcome")
  ## Coercing to numeric if not
  if(!is.numeric(dat$SE.outcome)) {
    warning("se column is not numeric. Coercing...")
    dat <- dat[, SE.outcome:= as.numeric(SE.outcome)]
  }
  ## replace infinite values (NA, NaN, Inf or -Inf), if any, with NA
  dat <- dat[, SE.outcome:= ifelse(is.finite(SE.outcome), SE.outcome, NA)]


  # Check pval ####
  ## pval column is given
  if(pval_col %in% names(dat)) {
    data.table::setnames(dat, pval_col, "PVAL.outcome")
    ## Coercing to numeric if not
    if(!is.numeric(dat$PVAL.outcome)) {
      warning("pval column is not numeric. Coercing...")
      dat <- dat[, PVAL.outcome:= as.numeric(PVAL.outcome)]
    }
    ## Transforming to the p-value if it was given as `-log10(p-value)`
    if (log_pval) {
      message("The p-value column is given as -log10(p-value). Transforming to the p-values ...")
      dat <- dat[, PVAL.outcome:= 10^-PVAL.outcome]
    }
    ## replace infinite values (NA, NaN, Inf or -Inf) and implausible values, if any, with NA
    dat <- dat[, PVAL.outcome:= ifelse(is.finite(PVAL.outcome) & PVAL.outcome >= 0 & PVAL.outcome <= 1, PVAL.outcome, NA)]
    ## replace extremely small values (smaller than `min_pval`), if any, with `min_pval`
    dat <- dat[, PVAL.outcome:= ifelse(PVAL.outcome < min_pval, min_pval, PVAL.outcome)]
    ## Label p-values as reported/inferred
    dat <- dat[, PVAL_origin.outcome:= ifelse(!is.na(PVAL.outcome), "reported", "inferred")]
    ## inferring p-value from `beta` and `se` if there is any missing p-values
    dat <- dat[, PVAL.outcome:= ifelse(!is.na(PVAL.outcome), PVAL.outcome, pchisq((BETA.outcome^2)/(SE.outcome^2), 1, lower.tail = FALSE))]

  ## No pval column is given - inferring it
  } else {
    message("Inferring p-values ...")
    dat <- dat[, PVAL.outcome:= pchisq((BETA.outcome^2)/(SE.outcome^2), 1, lower.tail = FALSE)]
    dat <- dat[, PVAL_origin.outcome:= "inferred"]
  }


  # Check eaf ####
  if(eaf_col %in% names(dat)) {
    data.table::setnames(dat, eaf_col, "EAF.outcome")
    ## Coercing to numeric if not
    if(!is.numeric(dat$EAF.outcome)) {
      warning("eaf column is not numeric. Coercing...")
      dat <- dat[, EAF.outcome:= as.numeric(EAF.outcome)]
    }
    ## replace infinite values (NA, NaN, Inf or -Inf), if any, with NA
    dat <- dat[, EAF.outcome:= ifelse(is.finite(EAF.outcome), EAF.outcome, NA)]
  } else {
    # Add in EAF col (as NA) if missing
    dat <- dat[, EAF.outcome:= NA]
  }


  # Check effect_allele ####
  data.table::setnames(dat, effect_allele_col, "EA.outcome")
  ## Coercing to character if not
  if(!is.character(dat$EA.outcome)) {
    warning("effect_allele column is not character data. Coercing...")
    dat <- dat[, EA.outcome:= as.character(EA.outcome)]
  }
  dat <- dat[, EA.outcome:= toupper(EA.outcome)]
  ## replace effect_allele inputs which are not A/C/T/G or an indel into NAs
  dat <- dat[, EA.outcome:= ifelse((grepl("^[ACTG]+$", EA.outcome) | EA.outcome %in% c("D", "I")), EA.outcome, NA)]
  if(any(is.na(dat$EA.outcome))) {
    warning("effect_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. Excluding ...")
  }


  # Check other_allele ####
  data.table::setnames(dat, other_allele_col, "OA.outcome")
  ## Coercing to character if not
  if(!is.character(dat$OA.outcome)) {
    warning("other_allele column is not character data. Coercing...")
    dat <- dat[, OA.outcome:= as.character(OA.outcome)]
  }
  dat <- dat[, OA.outcome:= toupper(OA.outcome)]
  ## replace other_allele inputs which are not A/C/T/G or an indel into NAs
  dat <- dat[, OA.outcome:= ifelse((grepl("^[ACTG]+$", OA.outcome) | OA.outcome %in% c("D", "I")), OA.outcome, NA)]
  if(any(is.na(dat$OA.outcome))) {
    warning("other_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. Excluding ...")
  }


  # Report gene ####
  if(gene_col %in% names(dat)) {
    data.table::setnames(dat, gene_col, "GENE.outcome")
  }

  # Report chr ####
  if(chr_col %in% names(dat)) {
    data.table::setnames(dat, chr_col, "CHR.outcome")
  }

  # Report position ####
  if(pos_col %in% names(dat)) {
    data.table::setnames(dat, pos_col, "POS.outcome")
  }


  # Indicator for SNPs with complete valid data for the SH analysis ####
  dat <- dat[, SH_keep.outcome:= ifelse(!(is.na(BETA.outcome) | is.na(SE.outcome) | is.na(EA.outcome) | is.na(OA.outcome)), TRUE, FALSE)]
  ## check if there's any valid SNPs for the SH
  if(any(dat$SH_keep.outcome)) {
    message("Identified ", sum(dat$SH_keep.outcome), " valid SNPs for the Slope-Hunter analyses!")
  } else {
    stop("None of the provided SNPs can be used for SH analysis. They are missing essentially required information!")
  }

  # Set order of the columns
  dat <- data.table::setcolorder(dat, c("SNP", "EA.outcome", "OA.outcome", "BETA.outcome", "SE.outcome",
                                        "PVAL.outcome", "PVAL_origin.outcome", "EAF.outcome"))

  # Rename columns
  data.table::setnames(dat, names(dat), gsub("outcome", type, names(dat)))
  rownames(dat) <- NULL
  return(dat)
}
