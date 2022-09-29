#' Read incidence data
#'
#' Reads in incidence data. Checks and organises columns for use with the Slope-Hunter analyses. Infers p-values when possible from beta and se.
#'
#' @md
#' @param filename Filename (formatted as .gz, .csv or .txt). Must have header with at least the `SNP`, `beta`, `se` and `EA`columns present.
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
#' @importFrom tools file_ext
#' @export
#' @return data frame
read_incidence = function(filename, snp_col="SNP", beta_col="BETA", se_col="SE",
                          pval_col="PVAL", eaf_col="EAF", effect_allele_col="EA",
                          other_allele_col="OA", gene_col="GENE", chr_col = "CHR", pos_col="POS",
                          min_pval=1e-200, log_pval=FALSE){

  incidence_dat <- data.table::fread(filename)

  incidence_dat <- format_data(
    incidence_dat,
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

  ext = tools::file_ext(filename)
  incidence_dat$data_source.incidence <- ifelse(ext == "gz", "compressed (gzip)", paste0("textfile (", ext, ")"))

  return(incidence_dat)
}
