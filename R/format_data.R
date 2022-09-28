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
#' @importFrom data.table data.table
#' @importFrom base tolower gsub subset duplicated
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

  i <- names(dat) %in% all_cols
  if (sum(i) == 0) {
    stop("None of the specified columns present!")
  }
  dat <- dat[, i, with = FALSE]    # data.table (FAQ 1.5): select column names stored in i

  # Check if columns required for SH are present
  SH_cols_req = c(snp_col, beta_col, se_col, effect_allele_col, other_allele_col)
  if (!all(SH_cols_req %in% names(dat))){
    stop("The following columns are not present and are required for the Slope-Hunter analysis\n", paste(SH_cols_req[!SH_cols_req %in% names(dat)]), collapse="\n")
  }

  # Format SNP IDs to the standard: lower-case letters, remove spaces & exclude NAs
  names(dat)[names(dat) == snp_col] <- "SNP"
  # snp_col <- "SNP"
  dat$SNP <- base::tolower(dat$SNP)
  dat$SNP <- base::gsub("[[:space:]]", "", dat$SNP)
  dat <- base::subset(dat, !is.na(SNP))

  if (!is.null(snps)){
    dat <- base::subset(dat, SNP %in% snps)
  }

  if (log_pval){
    dat$PVAL <- 10^-dat$PVAL
  }

  # Remove duplicated SNPs
  dup <- base::duplicated(dat$SNP)
  if (any(dup)){
    warning(sum(dup), " duplicated SNP rsIDs present: Just keeping the first instance ...")
    dat <- dat[!dup, ]
  }

  # initiate indicator for valid SNPs for the SH analysis
  dat$SH_keep.outcome <- TRUE

  # Check beta
  i <- which(names(dat) == beta_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "BETA.outcome"
    if(!is.numeric(dat$BETA.outcome))
    {
      warning("beta column is not numeric. Coercing...")
      dat$BETA.outcome <- as.numeric(dat$BETA.outcome)
    }
    index <- !is.finite(dat$BETA.outcome)
    index[is.na(index)] <- TRUE
    dat$BETA.outcome[index] <- NA
  }

  # Check se
  i <- which(names(dat) == se_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "SE.outcome"
    if(!is.numeric(dat$SE.outcome))
    {
      warning("se column is not numeric. Coercing...")
      dat$SE.outcome <- as.numeric(dat$SE.outcome)
    }
    index <- !is.finite(dat$SE.outcome) | dat$SE.outcome <= 0
    index[is.na(index)] <- TRUE
    dat$SE.outcome[index] <- NA
  }

  # Check pval
  i <- which(names(dat) == pval_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "PVAL.outcome"
    if(!is.numeric(dat$PVAL.outcome))
    {
      warning("pval column is not numeric. Coercing...")
      dat$PVAL.outcome <- as.numeric(dat$PVAL.outcome)
    }
    index <- !is.finite(dat$PVAL.outcome) | dat$PVAL.outcome < 0 | dat$PVAL.outcome > 1
    index[is.na(index)] <- TRUE
    dat$PVAL.outcome[index] <- NA
    index <- dat$PVAL.outcome < min_pval
    index[is.na(index)] <- FALSE
    dat$PVAL.outcome[index] <- min_pval
    dat$PVAL_origin.outcome <- "reported"
    if(any(is.na(dat$PVAL.outcome)))
    {
      index <- is.na(dat$PVAL.outcome)
      dat$PVAL.outcome[index] <- 2*pnorm(abs(dat$BETA.outcome[index])/dat$SE.outcome[index], lower.tail=FALSE)
      dat$PVAL_origin.outcome[index] <- "inferred"
    }
  }

  # If no pval column then create it from beta and se
  if(!"PVAL.outcome" %in% names(dat))
  {
    message("Inferring p-values ...")
    dat$PVAL.outcome <- 2*pnorm(abs(dat$beta.outcome)/dat$se.outcome, lower.tail=FALSE)
    dat$PVAL_origin.outcome <- "inferred"
  }

  # Check eaf
  i <- which(names(dat) == eaf_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "EAF.outcome"
    if(!is.numeric(dat$EAF.outcome))
    {
      warning("eaf column is not numeric. Coercing...")
      dat$EAF.outcome <- as.numeric(dat$EAF.outcome)
    }
    index <- !is.finite(dat$EAF.outcome) | dat$EAF.outcome <= 0 | dat$EAF.outcome >= 1
    index[is.na(index)] <- TRUE
    dat$EAF.outcome[index] <- NA
  }

  # Check effect_allele
  i <- which(names(dat) == effect_allele_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "EA.outcome"
    if(is.logical(dat$EA.outcome))
    {
      dat$EA.outcome <- substr(as.character(dat$EA.outcome), 1, 1)
    }
    if(!is.character(dat$EA.outcome))
    {
      warning("effect_allele column is not character data. Coercing...")
      dat$EA.outcome <- as.character(dat$EA.outcome)
    }

    dat$EA.outcome <- toupper(dat$EA.outcome)
    index <- !(grepl("^[ACTG]+$", dat$EA.outcome) | dat$EA.outcome %in% c("D", "I"))
    index[is.na(index)] <- TRUE
    if(any(index))
    {
      warning("effect_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. These SNPs will be excluded.")
      dat$EA.outcome[index] <- NA
      dat$SH_keep.outcome[index] <- FALSE
    }
  }

  # Check other_allele
  i <- which(names(dat) == other_allele_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "OA.outcome"
    if(is.logical(dat$OA.outcome))
    {
      dat$OA.outcome <- substr(as.character(dat$OA.outcome), 1, 1)
    }
    if(!is.character(dat$OA.outcome))
    {
      warning("other_allele column is not character data. Coercing...")
      dat$OA.outcome <- as.character(dat$OA.outcome)
    }

    dat$OA.outcome <- toupper(dat$OA.outcome)
    index <- ! (grepl("^[ACTG]+$", dat$OA.outcome) | dat$OA.outcome %in% c("D", "I"))
    index[is.na(index)] <- TRUE
    if(any(index))
    {
      warning("other_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. These SNPs will be excluded")
      dat$OA.outcome[index] <- NA
      dat$SH_keep.outcome[index] <- FALSE
    }
  }

  # Report gene
  if(gene_col %in% names(dat))
  {
    names(dat)[which(names(dat) == gene_col)[1]] <- "GENE.outcome"
  }

  # Report chr
  if(chr_col %in% names(dat))
  {
    names(dat)[which(names(dat) == chr_col)[1]] <- "CHR.outcome"
  }

  # Report position
  if(pos_col %in% names(dat))
  {
    names(dat)[which(names(dat) == pos_col)[1]] <- "POS.outcome"
  }

  # Identify SNPs with missing required info for SH analysis
  if(any(dat$SH_keep.outcome))
  {
    SHcols <- c("SNP", "BETA.outcome", "SE.outcome", "EA.outcome")
    SHcols_present <- SHcols[SHcols %in% names(dat)]
    dat$SH_keep.outcome <- dat$SH_keep.outcome & apply(dat[, SHcols_present, with=FALSE], 1, function(x) !any(is.na(x)))
    if(any(!dat$SH_keep.outcome))
    {
      warning("The following SNP(s) are missing required information for the SH analysis and will be excluded\n", paste(subset(dat, !SH_keep.outcome)$SNP, collapse="\n"))
    }
  }

  if(all(!dat$SH_keep.outcome))
  {
    warning("None of the provided SNPs can be used for SH analysis, they are missing required information.")
  }

  # Add in missing SH cols
  for(col in c("SNP", "BETA.outcome", "SE.outcome", "EA.outcome", "OA.outcome", "EAF.outcome"))
  {
    if(! col %in% names(dat))
    {
      dat[[col]] <- NA
    }
  }

  names(dat) <- gsub("outcome", type, names(dat))
  rownames(dat) <- NULL
  return(dat)
}
