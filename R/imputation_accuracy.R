#' Imputation accuracy, aka. correlations
#' 
#' Calculation of column-wise, row-wise, and matrix-wise correlations between
#' matrix in files \code{truefn} and \code{imputefn}.
#' Assumes first column in both files is an integer ID column and thus excluded from calculations.
#' Standardization (subtract mean, divide by standard deviation) is done column-wise based
#' on means and standard deviations of \code{truefn}.
#' Correlations are only performed on those rows that are found in \emph{both} files,
#' based on the first column (ID column).
#'
#' Genotypes equal to \code{NAval} are considered missing (i.e. \code{NA}) and are not included in the calculations.
#'
#' This method stores the "true" matrix in memory with a low-precision real type,
#' and rows in the "imputed" matrix are read and matched by ID.
#' If there are no extra rows in either matrix and order of IDs is the same,
#' consider setting \code{adaptive=FALSE}, as this has a memory usage of O(m), compared to O(nm) for the adaptive method, where 'm' is the number of SNPs and 'n' the number of animals.
#' The non-adaptive method is however, and very surprisingly, slightly slower.
#'
#' @param truefn Filename of \emph{true} matrix.
#' @param imputefn Filename of \emph{imputed} matrix.
#' @param ncol Integer, number of SNP columns in files. When \code{NULL}, automagically detected with \code{get_ncols(truefn)-1}.
#' @param nlines Integer, number of lines in \code{truefn}. When \code{NULL}, automagically detected with \code{gen_nlines(truefn)}.
#' @param na Value of missing genotypes.
#' @param standardized Logical, whether to center and scale genotypes by dataset in \code{true}-matrix.
#'        Currently by subtracting column mean and dividing by column standard deviation.
#' @param adaptive Use adaptive method (default) that stores \code{truefn} in memory and compares rows by ID in first column.
#' @return List with following elements:
#' \describe{
#'   \item{\code{means}}{Column means of true matrix.}
#'   \item{\code{vars}}{Column variances of true matrix.}
#'   \item{\code{rowcors}}{Row-wise (animal-wise) correlation between true and imputed matrix.}
#'   \item{\code{matcor}}{Matrix-wise correlation between true and imputed matrix.}
#'   \item{\code{colcors}}{Column-wise (locus-wise) correlation between true and imputed matrix.}
#'   \item{\code{rowID}}{Row IDs, corresponding to \code{rowcors}.}
#' }
#' @export
#' @seealso \code{\link{write.snps}} for writing SNPs to a file.
imputation_accuracy <- function(truefn, imputefn, ncol=NULL, nlines=NULL, na=9, standardized=TRUE, adaptive=TRUE) {
  stopifnot(file.exists(truefn))
  stopifnot(file.exists(imputefn))
  
  standardized <- as.logical(standardized)
  if (is.null(ncol)) ncol <- get_ncols(truefn)-1
  if (is.null(nlines)) nlines <- get_nlines(truefn)
  
  m <- as.integer(ncol)
  n <- as.integer(nlines)
  
  subroutine <- ifelse(adaptive, 'imp_acc', 'imp_acc_fast')
  
  res <- .Fortran(subroutine,
                  truefn=as.character(truefn),
                  imputedfn=as.character(imputefn),
                  nSnps=m,
                  nAnimals=as.integer(nlines),
                  NAval=as.integer(na),
                  standardized=as.integer(standardized),
                  means=vector('numeric',m), sds=vector('numeric',m),  # Placeholders for return data.
                  rowcors=vector('numeric', n), matcor=numeric(1), colcors=vector('numeric',m),
                  rowID=vector('integer',n),
                  PACKAGE='Siccuracy')
  res$colcors[is.infinite(res$colcors)] <- NA
  res$colcors[is.nan(res$colcors)] <- NA
  res$rowcors[is.infinite(res$rowcors)] <- NA
  res$rowcors[is.nan(res$rowcors)] <- NA
  if (standardized) res$means[(res$means- -9) < 1e-8] <- NA
  #res$sds[res$sds == 0.0] <- NA
  
  if (adaptive & any(is.na(res$rowcors))) {
    res$rowcors <- res$rowcors[!is.na(res$rowcors)]
    res$rowID <- res$rowID[!is.na(res$rowID)]
  }
  
  res[c('means','sds','rowcors','matcor','colcors','rowID')]
}

#' \code{imputation_accuracy1} and \code{imputation_accuracy3} has been replaced by \code{\link{imputation_accuracy}}.
#' The difference between the two former functions is now covered by the \code{adaptive}-argument of the latter.
#' 
#' @param nSNPs Deprecated, use \code{ncol}.
#' @param nAnimals Deprecated, use \code{nlines}.
#' @param NAval Deprecated, use \code{na}.
#' 
#' @export
#' @rdname deprecated
#' @inheritParams imputation_accuracy
imputation_accuracy3 <- function(truefn, imputefn, nSNPs=NULL, nAnimals=NULL, NAval=9, standardized=TRUE) {
  .Deprecated('imputation_accuracy', package='Siccuracy')
  imputation_accuracy(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, adaptive=TRUE)
}

#' @export
#' @rdname deprecated
imputation_accuracy1 <- function(truefn, imputefn, nSNPs=NULL, nAnimals=NULL, NAval=9, standardized=TRUE) {
  .Deprecated('imputation_accuracy', package='Siccuracy')
  imputation_accuracy(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, adaptive=FALSE)
}