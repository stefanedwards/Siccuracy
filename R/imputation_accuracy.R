#' Imputation accuracy, aka. correlations
#' 
#' Calculation of column-wise, row-wise, and matrix-wise correlations between
#' matrix in file \code{truefn} and \code{imputefn}.
#' Assumes first column in both files is an integer ID column and thus excluded.
#' Standardization (subtract mean, divide by variance) is done column-wise based
#' on means and variances of \code{truefn}.
#' Correlations are only performed on those rows that are found in \emph{both} files,
#' based on the first column (ID column).
#'
#' This method stores the "true" matrix in memory with a low-precision real type,
#' and rows in the "imputed" matrix are read and matched by ID.
#' If there are no extra rows in either matrix and order of IDs is the same,
#' consider using \code{fast}-option that does not store either.
#'
#' @param truefn Filename to true matrix. NB! Max. 255 characters!
#' @param imputefn Filename to imputed matrix.
#' @param nSNPs Integer, number of snps (i.e. number of columns minus 1).
#'        If \code{NULL} (default), detected from \code{truefn}.
#' @param nAnimals Integer, number of rows.
#'        If \code{NULL} (default), detected from \code{truefn}.
#' @param NAval Integer, value of missing genotype.
#' @param standardized Logical, whether to center and scale genotypes by dataset in \code{true}-matrix.
#'        Currently by subtracting column mean and dividing by column variance.
#' @param fast Use method that does not check row IDs.
#' @return List with following elements:
#' \describe{
#'   \item{\code{means}}{Column means of true matrix.}
#'   \item{\code{vars}}{Column variances of true matrix.}
#'   \item{\code{rowcors}}{Row-wise (animal-wise) correlation between true and imputed matrix.}
#'   \item{\code{matcor}}{Matrix-wise correlation between true and imputed matrix.}
#'   \item{\code{colcors}}{Column-wise (locus-wise) correlation between true and imputed matrix.}
#' }
#' @export
#' @seealso \code{\link{write.snps}} for writing SNPs to a file.
imputation_accuracy <- function(truefn, imputefn, nSNPs=NULL, nAnimals=NULL, NAval=9, standardized=TRUE, fast=FALSE) {
  stopifnot(file.exists(truefn))
  stopifnot(file.exists(imputefn))
  
  standardized <- as.logical(standardized)
  if (is.null(nSNPs)) {
    nSNPs <- get_ncols(truefn)-1
  }
  if (is.null(nAnimals)) {
    nAnimals <- get_nlines(truefn)
  }
  m <- as.integer(nSNPs)
  n <- as.integer(nAnimals)
  
  subroutine <- ifelse(fast, 'imp_acc_fast','imp_acc')
  
  res <- .Fortran(subroutine,
                  truefn=as.character(truefn),
                  imputedfn=as.character(imputefn),
                  nSnps=m,
                  nAnimals=as.integer(nAnimals),
                  NAval=as.integer(NAval),
                  standardized=as.integer(standardized),
                  means=vector('numeric',m), vars=vector('numeric',m),  # Placeholders for return data.
                  rowcors=vector('numeric', n), matcor=numeric(1), colcors=vector('numeric',m),
                  rowcorID=vector('integer',m),
                  PACKAGE='Siccuracy')
  res$colcors[is.infinite(res$colcors)] <- NA
  res$colcors[is.nan(res$colcors)] <- NA
  res$rowcors[is.infinite(res$rowcors)] <- NA
  res$rowcors[is.nan(res$rowcors)] <- NA
  res[c('means','vars','rowcors','matcor','colcors','rowcorID')]
}

#' Deprecated function names
#' @export
#' @rdname depcreated_impacc
#' @inheritParams imputation_accuracy
imputation_accuracy3 <- function(truefn, imputefn, nSNPs=NULL, nAnimals=NULL, NAval=9, standardized=TRUE) {
  Deprecated('imputation_accuracy', package='Siccuracy')
  imputation_accuracy(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, fast=FALSE)
}

#' Deprecated function names
#' @export
#' @rdname depcreated_impacc
imputation_accuracy1 <- function(truefn, imputefn, nSNPs=NULL, nAnimals=NULL, NAval=9, standardized=TRUE) {
  Deprecated('imputation_accuracy', package='Siccuracy')
  imputation_accuracy(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, fast=TRUE)
}