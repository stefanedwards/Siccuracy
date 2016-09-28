# Phase-o-to-geno #####################

#' Convert phase files to genotype files.
#' 
#' Simply sums every pair of rows together.
#' 
#' A phase file has format similar to SNP files, expect the genotype on each allele are listed on two pairs of rows.
#' Each individual has therefore two rows, one for each allele, see following example:
#'
#' \code{Genotype file:}
#' \tabular{lccccc}{
#'  1003 \tab 0 \tab 1 \tab 1 \tab 2 \tab 0 \cr
#'  1004 \tab 1 \tab 1 \tab 0 \tab 1 \tab 0 \cr
#' }
#'
#' \code{Phase file:}
#' \tabular{lccccc}{
#'  1003 \tab 0 \tab 0 \tab 1 \tab 1 \tab 0 \cr
#'  1003 \tab 0 \tab 1 \tab 0 \tab 1 \tab 0 \cr
#'  1004 \tab 0 \tab 1 \tab 0 \tab 0 \tab 0 \cr
#'  1004 \tab 1 \tab 0 \tab 0 \tab 1 \tab 0 \cr
#' }
#'
#' By default assumes genotypes are whole integers (cf. \code{int} argument).
#' If numeric values are required, e.g. gene dosages from imputation, use \code{int=FALSE}. 
#' Use \code{numeric.format} to change from default 5 characters width and 2 decimals (\code{"5.2"}). 
#' NB! Function does not test for validity of this argument, so change with caution!
#'
#' @param phasefn Filename of input file, every two rows are for same animal.
#' @param genofn Filename of intended output.
#' @param ncol Number of SNPs in file, if \code{NULL} (default), number is estimated from file.
#' @param nrow Number of rows to maximally read from \code{phasefn}. If \code{NULL}, no limit is used.
#' @param naval Missing values; genotypes equal to or larger than this value are set to this value.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param numeric.format Character describing \code{<width>.<decimals>}. Default \code{'5.2'}. Has no effect when \code{int=TRUE}.
#' @return Number of rows written.
#' @export
convert_phases <- function(phasefn, genofn, ncol=NULL, nrow=NULL, int=TRUE, naval=9, numeric.format='5.2') {
  stopifnot(file.exists(phasefn))
  
  if (is.null(ncol)) ncol <- get_ncols(phasefn)-1
  if (is.null(nrow)) nrow <- 0
  
  subroutine <- ifelse(int, 'convert_phase_int', 'convert_phase_num')
  res <- .Fortran(subroutine, phasefn=as.character(phasefn), genofn=as.character(genofn), ncol=as.integer(ncol), nrow=as.integer(nrow), naval=as.integer(naval), numfmt=as.character(numeric.format))
  res$nrow
}

#' @export
#' @rdname depcreated
#' @inheritParams imputation_accuracy
phasotogeno <- function(phasefn, genofn, ncol=NULL, nrow=NULL) {
  Deprecated('imputation_accuracy', package='Siccuracy')
  convert_phases(phasefn, genofn, ncol, nrow, int=FALSE)
}
#' @export
#' @rdname depcreated
#' @inheritParams imputation_accuracy
phasotogeno_int <- function(phasefn, genofn, ncol=NULL, nrow=NULL) {
  Deprecated('imputation_accuracy', package='Siccuracy')
  convert_phases(phasefn, genofn, ncol, nrow, int=TRUE)
}



# Convert plink A format #########################