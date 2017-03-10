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
#' @param center Numeric vector of \code{ncol}-length to subtract with for standardization.
#' @param scale Numeric vector of \code{ncol}-length to divide by for standardization.
#' @param p Shortcut for \code{center} and \code{scale} when using allele frequencies. \code{center=2p} and \code{scale=2p(1-p)}.
#' @param excludeIDs Integer vector, exclude these individuals from correlations. \emph{Does not affect calculation of column means and standard deviations.}
#' @param excludeSNPs Integer or logical vector, exclude these columns from correlations. \emph{Does not affect calculation of column means and standard deviations.}
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
imputation_accuracy <- function(truefn, imputefn, ncol=NULL, nlines=NULL, na=9, standardized=TRUE, adaptive=TRUE, center=NULL, scale=NULL, p=NULL, excludeIDs=NULL, excludeSNPs=NULL, tol=0.1) {
  stopifnot(file.exists(truefn))
  stopifnot(file.exists(imputefn))
  
  standardized <- as.logical(standardized)
  if (is.null(ncol)) ncol <- get_ncols(truefn)-1
  if (is.null(nlines)) nlines <- get_nlines(truefn)
  
  m <- as.integer(ncol)
  n <- as.integer(nlines)
  
  usermeans <-  (!is.null(p) | !is.null(center) | !is.null(scale))
  if (!is.null(p)) {
    stopifnot(length(p)==m)
    center <- 2*p
    scale <- 2*p*(1-p)
  }
  if (is.null(center)) center=numeric(m)
  if (is.null(scale)) {scale=numeric(m);scale[] <- 1}
  if (usermeans) {
    center[is.na(center)] <- 0
    scale[is.na(scale)] <- 0
    standardized=TRUE
  }
  
  ex_ids <- rep(0, n)
  if (!is.null(excludeIDs)) {
    imp_ids <- get_firstcolumn(truefn)
    ex_ids[imp_ids %in% excludeIDs] <- 1 
  }
  
  ex_snps <- rep(0, m)
  if (!is.null(excludeSNPs)) {
    .is.na <- function(x) {if (is.null(x)) return(logical(0)); is.na(x)}
    if (is.logical(excludeSNPs) & sum(!.is.na(excludeSNPs)) < m) stop('`excludeSNPs` as logical must be same length as SNPs in input files and without NA\'s.')
    
    if (is.logical(excludeSNPs)) {
      ex_snps <- as.integer(excludeSNPs)
    } else {
      ex_snps[excludeSNPs] <- 1
    }
  }
  
  subroutine <- ifelse(adaptive, 'imp_acc', 'imp_acc_fast')
  
  res <- .Fortran(subroutine,
                  truefn=as.character(truefn),
                  imputedfn=as.character(imputefn),
                  nSnps=m,
                  nAnimals=as.integer(nlines),
                  NAval=as.integer(na),
                  standardized=as.integer(standardized),
                  means=as.numeric(center), sds=as.numeric(scale),  # Placeholders for return data.
                  usermeans=as.integer(usermeans),
                  rowcors=vector('numeric', n), matcor=numeric(1), colcors=vector('numeric',m),
                  rowID=vector('integer',n),
                  excludeIDs=as.integer(ex_ids),
                  excludeSNPs=as.integer(ex_snps),
                  tol=as.numeric(tol),
                  colcorrect=vector('integer',m), coltruena=vector('integer',m), colimpna=vector('integer',m), colbothna=vector('integer',m), 
                  rowcorrect=vector('integer',n), rowtruena=vector('integer',n), rowimpna=vector('integer',n), rowbothna=vector('integer',n),
                  NAOK=FALSE,
                  PACKAGE='Siccuracy')
  res$colcors[is.infinite(res$colcors)] <- NA
  res$colcors[is.nan(res$colcors)] <- NA
  res$rowcors[is.infinite(res$rowcors)] <- NA
  res$rowcors[is.nan(res$rowcors)] <- NA
  if (standardized) res$means[(res$means- -9) < 1e-8] <- NA
  #res$sds[res$sds == 0.0] <- NA
  
  if (adaptive & any(is.na(res$rowID))) {
    res$rowcors <- res$rowcors[!is.na(res$rowID)]
    res$rowID <- res$rowID[!is.na(res$rowID)]
  }
  
  if (!is.null(excludeIDs)) res$rowcors[ex_ids==1] <- NA
  
  
  #res[c('means','sds','rowcors','matcor','colcors','rowID')]
  with(res, 
       list(matcor=matcor, snps=data.frame(means, sds, cors=colcors, correct=colcorrect, true.na=coltruena, imp.na=colimpna, both.na=colbothna),
            animals=data.frame(rowID, cors=rowcors, correct=rowcorrect, true.na=rowtruena, imp.na=rowimpna, both.na=rowbothna)))
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
  with(imputation_accuracy(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, adaptive=TRUE),
       list(means=snps$means, sds=snps$sds, rowcors=animals$cors, matcor=matcor))
}

#' @export
#' @rdname deprecated
imputation_accuracy1 <- function(truefn, imputefn, nSNPs=NULL, nAnimals=NULL, NAval=9, standardized=TRUE) {
  .Deprecated('imputation_accuracy', package='Siccuracy')
  with(imputation_accuracy(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, adaptive=FALSE),
       list(means=snps$means, sds=snps$sds, rowcors=animals$cors, matcor=matcor))
}