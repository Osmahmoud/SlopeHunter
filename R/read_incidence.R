#' Read incidence data
#'
#' Reads in incidence data. Checks and organises columns for use with the Slope-Hunter analyses. Infers p-values when possible from beta and se.
#'
#' @md
#' @param filename Filename. Must have header with at least the `SNP`, `beta`, `se` and `EA`columns present.
#' @param gz whether the given data file is compressed with the gzip (`.gz` file). The default is `TRUE`.
#' @param sep Specify delimeter in file if `gz` is FALSE. The default is a space, i.e. `" "`.
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
#' @importFrom data.table fread
#' @importFrom utils read.table
#' @export
#' @return data frame
read_incidence = function(filename, gz = TRUE, sep= " ", snp_col="SNP", beta_col="BETA", se_col="SE",
                          pval_col="PVAL", eaf_col="EAF", effect_allele_col="EA",
                          other_allele_col="OA", gene_col="GENE", chr_col = "CHR", pos_col="POS",
                          min_pval=1e-200, log_pval=FALSE){

  if(gz){incidence_dat <- utils::read.table(gzfile(filename), header = TRUE)} else {
    incidence_dat <- data.table::fread(filename, header=TRUE, sep=sep)
  }

  incidence_dat <- format_data(
    as.data.frame(incidence_dat),
    type="incidence",
    snps=NULL,
    snp_col=snp_col,
    beta_col=beta_col,
    se_col=se_col,
    pval_col=pval_col,
    eaf_col=eaf_col,
    effect_allele_col=effect_allele_col,
    other_allele_col=other_allele_col,
    gene_col=gene_col,
    chr_col=chr_col,
    pos_col=pos_col,
    min_pval=min_pval,
    log_pval=log_pval
  )

  incidence_dat$data_source.incidence <- ifelse(gz, "compressed (gzip)", "textfile")
  return(incidence_dat)
}



