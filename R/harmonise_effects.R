#' Harmonise and format data for Slope-Hunter
#'
#' Harmonise the alleles and effects between the incidence and prognosis
#'
#' In order to perform Slope-Hunter analysis the effect of a SNP on an incidence and prognosis traits must be harmonised to be
#' relative to the same allele.
#'
#' @md
#' @param incidence_dat data.table for incidence data. It is recommended to be an output from \code{read_incidence}. If not, it tries to format it before harmonisation.
#' @param prognosis_dat data.table for prognosis data. It is recommended to be an output from \code{read_prognosis}. If not, it tries to format it before harmonisation.
#' @param incidence_formatted Logical indicationg whether \code{incidence_dat} is formatted using \code{read_incidence}.
#' @param prognosis_formatted Logical indicationg whether \code{prognosis_dat} is formatted using \code{read_prognosis}.
#' @param by.pos Logical, if `TRUE` the harmonisation will be performed by matching the exact SNP positions between the incidence and prognosis datasets.
#' @param pos_cols A vector of length 2 specifying the name of the genetic position columns in the incidence and prognosis datasets respectively.
#' @param snp_cols A vector of length 2 specifying the name of the snp columns in the incidence and prognosis datasets respectively. This is the column on which the data will be merged if \code{by.pos} is `FASLE`.
#' @param beta_cols A vector of length 2 specifying the name of the beta columns in the incidence and prognosis datasets respectively.
#' @param se_cols A vector of length 2 specifying the name of the se columns in the incidence and prognosis datasets respectively.
#' @param EA_cols A vector of length 2 specifying the name of the effect allele columns in the incidence and prognosis datasets respectively.
#' @param OA_cols A vector of length 2 specifying the name of the non-effect allele columns in the incidence and prognosis datasets respectively.
#' @param chr_cols A vector of length 2 specifying the name of the chromosome columns in the incidence and prognosis datasets respectively.
#' @param gene_col A vector of length 2 specifying the name of the gene columns in the incidence and prognosis datasets respectively.
#' @details This function will try to harmonise the incidence and prognosis data sets on the specified columns. Where necessary, correct strand
#' for non-palindromic SNPs (i.e. flip the sign of effects so that the effect allele is the same in both datasets), and drop all palindromic SNPs from the analysis (i.e. with the allele A/T or G/C).
#' The alleles that do not match between data sets (e.g T/C in one data set and A/C in the other) will also be dropped.
#'
#' @return A data.frame with harmonised effects and alleles
#' @export

harmonise_effects <- function(incidence_dat, prognosis_dat, incidence_formatted=TRUE, prognosis_formatted=TRUE,
                              by.pos=FALSE, pos_cols=c("POS.incidence", "POS.prognosis"), snp_cols=c("SNP", "SNP"), beta_cols=c("BETA.incidence", "BETA.prognosis"),
                              se_cols=c("SE.incidence", "SE.prognosis"), EA_cols=c("EA.incidence", "EA.prognosis"), OA_cols=c("OA.incidence", "OA.prognosis"),
                              chr_cols=c("CHR.incidence", "CHR.prognosis"), gene_col=c("GENE.incidence", "GENE.prognosis")) {

  # check inputs of col names are given for both incidence and prognosis ####
  L_pos = length(pos_cols) == 2
  L_snp = length(snp_cols) == 2
  L_beta = length(beta_cols) == 2
  L_se = length(se_cols) == 2
  L_ea = length(EA_cols) == 2
  L_oa = length(OA_cols) == 2
  L_chr = length(chr_cols) == 2
  L_gene = length(gene_col) == 2
  CheCol = c(L_pos, L_snp, L_beta, L_se, L_ea, L_oa, L_chr, L_gene)
  if(!all(CheCol)) {
    stop("column names must be given as vectors of lengths 2!")
  }

  # Format incidence if not formatted ####
  if(!incidence_formatted){
    message("Formatting incidence data ...\n")
    incidence_dat <- format_data(
      as.data.table(incidence_dat),
      type="incidence",
      snps=NULL,
      snp_col=snp_cols[1],
      beta_col=beta_cols[1],
      se_col=se_cols[1],
      effect_allele_col=EA_cols[1],
      other_allele_col=OA_cols[1],
      pos_col=pos_cols[1],
      chr_col=chr_cols[1]
    )
  }

  # Format prognosis if not formatted ####
  if(!prognosis_formatted){
    message("Formatting prognosis data ...\n")
    prognosis_dat <- format_data(
      as.data.table(prognosis_dat),
      type="prognosis",
      snps=NULL,
      snp_col=snp_cols[2],
      beta_col=beta_cols[2],
      se_col=se_cols[2],
      effect_allele_col=EA_cols[2],
      other_allele_col=OA_cols[2],
      pos_col=pos_cols[2],
      chr_col=chr_cols[2]
    )
  }

  # merging by SNPs ####
  if(!by.pos){
    message("Merging incidence and prognosis datasets by SNPs ...")
    dat <- merge(incidence_dat, prognosis_dat, by="SNP")
  }

  # merging by POS ####
  if(by.pos){
    message("Merging incidence and prognosis datasets by genomic position ...")
    ## stop if pos_cols are not present in the data
    if(!((pos_cols[1] %in% names(incidence_dat)) & (pos_cols[2] %in% names(prognosis_dat)))) {
      stop("You asked to harmonise by genomic position, but pos_cols are not present in the data!")
    }
    dat <- merge(incidence_dat, prognosis_dat, by.x="POS.incidence", by.y="POS.prognosis")
  }

  message("Harmonising effects and alleles ...")

  #TEMP!!!!!!!1
  ## give unique names for common columns (e.g. SNP, CHR & POS)
  names(dat)[names(dat) == "CHR.incidence"] <- "CHR"
  names(dat)[names(dat) == "POS.incidence"] <- "POS"

  ## give unique names for common columns (e.g. SNP, CHR & POS)
  if(sum(grepl("^rs\\d+$", dat$SNP.x)) > sum(grepl("^rs\\d+$", dat$SNP.y))){
    names(dat)[names(dat) == "SNP.x"] <- "SNP"
  } else {
    names(dat)[names(dat) == "SNP.y"] <- "SNP"
  }
  names(dat)[names(dat) == "CHR.incidence"] <- "CHR"
  names(dat)[names(dat) == "POS.incidence"] <- "POS"
  #END TEMP!!!!!!!!!!


  dat <- harmonise_dataset(dat)
  return(dat)
}



#'
#' @importFrom data.table data.table
harmonise_dataset <- function(dat)
{
  A1=dat$EA.incidence
  A2=dat$OA.incidence
  B1=dat$EA.prognosis
  B2=dat$OA.prognosis
  betaA=dat$BETA.incidence
  betaB=dat$BETA.prognosis
  seA=dat$SE.incidence
  seB=dat$SE.prognosis
  pA=dat$PVAL.incidence
  pB=dat$PVAL.prognosis
  fA=dat$EAF.incidence
  fB=dat$EAF.prognosis

  if(nrow(dat) == 0) {
   stop("No common genetic variants present in incidence and prognosis data ...")
  }

  jinfo <- list()

  # indel recoding
  indel_index <- nchar(A1) > 1 | nchar(A2) > 1 | A1 == "D" | A1 == "I"
  temp <- recode_indels(A1[indel_index], A2[indel_index], B1[indel_index], B2[indel_index])

  A1[indel_index] <- temp$A1
  A2[indel_index] <- temp$A2
  B1[indel_index] <- temp$B1
  B2[indel_index] <- temp$B2
  jinfo[['indel_kept']] <- sum(temp$keep)
  jinfo[['indel_removed']] <- sum(!temp$keep)

  # Find SNPs with alleles that match in A and B (OK, e.g. A/C Vs. A/C)
  OK <- (A1 == B1) & (A2 == B2)

  # Find SNPs with alleles in B need to swap (e.g. A/C Vs. C/A)
  to_swap <- (A1 == B2) & (A2 == B1)
  jinfo[['switched_alleles']] <- sum(to_swap)

  # For B's alleles that need to swap, do swap
  Btemp <- B1[to_swap]
  B1[to_swap] <- B2[to_swap]
  B2[to_swap] <- Btemp
  betaB[to_swap] <- betaB[to_swap] * -1
  fB[to_swap] <- 1 - fB[to_swap]

  # Check Again
  OK <- (A1 == B1) & (A2 == B2)
  palindromic <- is.palindromic(A1, A2)
  jinfo[['palindromic']] <- sum(palindromic)

  # For 'NON-palindromic and alleles still DON'T match' (e.g. A/C Vs. T/G), do try flipping
  i <- !palindromic & !OK
  B1[i] <- flip_alleles(B1[i])
  B2[i] <- flip_alleles(B2[i])
  OK <- (A1 == B1) & (A2 == B2)
  jinfo[['flipped_alleles']] <- sum(i)

  # If still DON'T match (e.g. A/C Vs. G/T, which - after flipping in the previous step - becomes A/C Vs. C/A) then try swapping
  i <- !palindromic & !OK
  to_swap <- (A1 == B2) & (A2 == B1)
  Btemp <- B1[i & to_swap]
  B1[i & to_swap] <- B2[i & to_swap]
  B2[i & to_swap] <- Btemp
  betaB[i & to_swap] <- betaB[i & to_swap] * -1
  fB[i & to_swap] <- 1 - fB[i & to_swap]
  jinfo[['flipped_then_switched_alleles']] <- sum(i & to_swap)

  # Any SNPs left with un-matching alleles (e.g. A/C Vs. A/G) need to be removed
  OK <- (A1 == B1) & (A2 == B2)
  remove <- !OK
  remove[indel_index][!temp$keep] <- TRUE   # remove invalid indel recoding

  # Remove palindromic SNPs
  remove[palindromic] <- TRUE

  # check which ID cols are present in data
  snp_cols <- grep("^SNP", names(dat), ignore.case = FALSE, value = TRUE)
  pos_cols <- grep("^POS", names(dat), ignore.case = FALSE, value = TRUE)
  chr_cols <- grep("^CHR", names(dat), ignore.case = FALSE, value = TRUE)
  gene_cols <- grep("^GENE", names(dat), ignore.case = FALSE, value = TRUE)
  cols <- c(snp_cols, pos_cols, chr_cols, gene_cols)
  # Stop if nothing of the ID cols present
  if(length(cols) == 0) {
    stop("No columns for either SNP, POS, CHR or GENE are present!")
  }
  # Extract the ID cols, then arrange the harmonised data
  dat <- dat[, cols, with = FALSE]
  dat <- dat[, ':='(EA.incidence=A1, OA.incidence=A2, EA.prognosis=B1, OA.prognosis=B2,
                    BETA.incidence=betaA, SE.incidence=seA, Pval.incidence=pA, EAF.incidence=fA,
                    BETA.prognosis=betaB, SE.prognosis=seB, Pval.prognosis=pB, EAF.prognosis=fB,
                    remove=remove, palindromic=palindromic)]

  dat <- as.data.frame(dat)
  attr(dat, "info") <- jinfo
  return(dat)
}

is.palindromic <- function(A1, A2)
{
  (A1 == "T" & A2 == "A") |
  (A1 == "A" & A2 == "T") |
  (A1 == "G" & A2 == "C") |
  (A1 == "C" & A2 == "G")
}

flip_alleles <- function(x)
{
  x <- toupper(x)
  x <- gsub("C", "g", x)
  x <- gsub("G", "c", x)
  x <- gsub("A", "t", x)
  x <- gsub("T", "a", x)
  return(toupper(x))
}

recode_indels <- function(A1, A2, B1, B2)
{
  ncA1 <- nchar(A1)
  ncA2 <- nchar(A2)
  ncB1 <- nchar(B1)
  ncB2 <- nchar(B2)

  i1 <- ncA1 > ncA2 & B1 == "I" & B2 == "D"
  B1[i1] <- A1[i1]
  B2[i1] <- A2[i1]

  i1 <- ncA1 < ncA2 & B1 == "I" & B2 == "D"
  B1[i1] <- A2[i1]
  B2[i1] <- A1[i1]

  i1 <- ncA1 > ncA2 & B1 == "D" & B2 == "I"
  B1[i1] <- A2[i1]
  B2[i1] <- A1[i1]

  i1 <- ncA1 < ncA2 & B1 == "D" & B2 == "I"
  B1[i1] <- A1[i1]
  B2[i1] <- A2[i1]


  i1 <- ncB1 > ncB2 & A1 == "I" & A2 == "D"
  A1[i1] <- B1[i1]
  A2[i1] <- B2[i1]

  i1 <- ncB1 < ncB2 & A1 == "I" & A2 == "D"
  A2[i1] <- B1[i1]
  A1[i1] <- B2[i1]

  i1 <- ncB1 > ncB2 & A1 == "D" & A2 == "I"
  A2[i1] <- B1[i1]
  A1[i1] <- B2[i1]

  i1 <- ncB1 < ncB2 & A1 == "D" & A2 == "I"
  A1[i1] <- B1[i1]
  A2[i1] <- B2[i1]

  keep <- rep(TRUE, length(A1))
  keep[(ncA1 > 1 & ncA1 == ncA2 & (B1 == "D" | B1 == "I"))] <- FALSE
  keep[(ncB1 > 1 & ncB1 == ncB2 & (A1 == "D" | A1 == "I"))] <- FALSE
  keep[A1 == A2] <- FALSE
  keep[B1 == B2] <- FALSE

  return(data.frame(A1=A1, A2=A2, B1=B1, B2=B2, keep=keep, stringsAsFactors=FALSE))
}
