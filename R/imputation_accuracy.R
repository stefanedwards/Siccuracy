
#' @rdname imputation_accuracy
#' @export
imputation_accuracy <- function(true, impute, ...) {
  UseMethod('imputation_accuracy', true)
}


#' Imputation accuracy, aka. correlations
#' 
#' Calculation of column-wise, row-wise, and matrix-wise correlations between
#' matrix in files \code{true} and \code{impute}.
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
#' @param true \emph{True} genotype matrix, or filename.
#' @param impute \emph{Imputed} genotype matrix, or filename.
#' @param ncol Integer, number of SNP columns in files. When \code{NULL}, automagically detected with \code{get_ncols(true)-1}.
#' @param nlines Integer, number of lines in \code{true}. When \code{NULL}, automagically detected with \code{gen_nlines(true)}.
#' @param na Value of missing genotypes.
#' @param standardized Logical, whether to center and scale genotypes by dataset in \code{true}-matrix.
#'        Currently by subtracting column mean and dividing by column standard deviation.
#' @param adaptive Use adaptive method (default) that stores \code{true} in memory and compares rows by ID in first column.
#' @param center Numeric vector of \code{ncol}-length to subtract with for standardization.
#' @param scale Numeric vector of \code{ncol}-length to divide by for standardization.
#' @param p Shortcut for \code{center} and \code{scale} when using allele frequencies. \code{center=2p} and \code{scale=2p(1-p)}.
#' @param excludeIDs Integer vector, exclude these individuals from correlations. \emph{Does not affect calculation of column means and standard deviations.}
#' @param excludeSNPs Integer or logical vector, exclude these columns from correlations. \emph{Does not affect calculation of column means and standard deviations.}
#' @param tol Numeric, tolerance for imputation error when counting correctly imputed genotypes.
#' @return List with following elements:
#' \describe{
#'   \item{\code{matcor}}{Matrix-wise correlation between true and imputed matrix.}
#'   \item{\code{snps}}{Data frame with all snp-wise statistics}
#'   \item{\code{animals}}{Data frame with all animal-wise statistics}
#' }
#' The data frames with statistics consists of columns
#' \describe{
#'   \item{\code{rowID}}{Row ID (\code{$animals} only!).}
#'   \item{\code{means}}{Value subtracted from each column (\code{$snps} only!).}
#'   \item{\code{sds}}{Value used to scale each column (i.e. standard deviations) (\code{$snps} only!).}
#'   \item{\code{cors}}{Pearson correlation between true and imputed genotype.}
#'   \item{\code{correct}}{Number of entries of equal value (within \code{tol})}
#'   \item{\code{true.na}}{Number of entries in that were missing in \code{true} but not \code{impute}.}
#'   \item{\code{imp.na}}{As \code{true.na}, but vice versa.}
#'   \item{\code{both.na}}{Number of entries that were missing in both files.}
#'   \item{\code{correct.pct}}{\code{correct} divided by total number of entries bare missing entries in \code{true}.}
#' }
#' @export
#' @rdname imputation_accuracy
#' @seealso \code{\link{write.snps}} for writing SNPs to a file.
imputation_accuracy.character <- function(true, impute, ncol=NULL, nlines=NULL, na=9, standardized=TRUE, adaptive=TRUE, center=NULL, scale=NULL, p=NULL, excludeIDs=NULL, excludeSNPs=NULL, tol=0.1) {
  stopifnot(file.exists(true))
  stopifnot(file.exists(impute))
  
  standardized <- as.logical(standardized)
  if (is.null(ncol)) ncol <- get_ncols(true)-1
  if (is.null(nlines)) nlines <- get_nlines(true)
  
  m <- as.integer(ncol)
  n <- as.integer(nlines)
  
  stopifnot(m == (get_ncols(impute) - 1))
  
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
  tru_ids <- get_firstcolumn(true)
  imp_ids <- get_firstcolumn(impute)
  ex_ids[!tru_ids %in% imp_ids] <- 1
  if (!is.null(excludeIDs)) {
    ex_ids[tru_ids %in% excludeIDs] <- 1 
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
                  true=as.character(true),
                  imputedfn=as.character(impute),
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
  res$sds[abs(res$sds) < 1e-6 & ((res$colbothna + res$coltruena) == n)] <- NA
  
  if (adaptive & any(is.na(res$rowID))) {
    keep <- !is.na(res$rowID)
    res$rowcors <- res$rowcors[keep]
    res$rowID <- res$rowID[!is.na(res$rowID)]
    res$rowcorrect <- res$rowcorrect[keep]
    res$rowimpna <- res$rowimpna[keep]
    res$rowtruena <- res$rowtruena[keep]
    res$rowbothna <- res$rowbothna[keep]
  }
  
  if (adaptive &  (any((res$colcorrect+res$coltruena+res$colimpna+res$colbothna) > n) | any((res$rowcorrect+res$rowtruena+res$rowimpna+res$rowbothna) > m))) {
    warning('Oh oh. It seems I have counted more elements than SNPs or animals.\n  Check `with(snps, correct+true.na+imp.na+both.na)` and `with(animals, correct+true.na+imp.na+both.na)`.\n  Do you perhaps have repeated IDs?')
  }
  
  if (!is.null(excludeIDs)) res$rowcors[ex_ids==1] <- NA
  
  
  #res[c('means','sds','rowcors','matcor','colcors','rowID')]
  with(res, 
       list(matcor=matcor, snps=data.frame(means, sds, cors=colcors, correct=colcorrect, true.na=coltruena, imp.na=colimpna, both.na=colbothna, correct.pct=colcorrect/(n-coltruena-colbothna-sum(ex_ids))),
            animals=data.frame(rowID, cors=rowcors, correct=rowcorrect, true.na=rowtruena, imp.na=rowimpna, both.na=rowbothna, correct.pct=rowcorrect/(m-rowtruena-rowbothna-sum(ex_snps)))))
}


#' Correlation between two matrices in R
#' @param true Matrix of true genotypes
#' @param impute Matrix of imputed genotypes
#' @param transpose Logical, if SNPs are per row, set to \code{TRUE}.
#' @inheritParams imputation_accuracy.character
#' @rdname imputation_accuracy
#' @export
imputation_accuracy.matrix <- function(true, impute, na=9, standardized=TRUE, center=NULL, scale=NULL, p=NULL, excludeIDs=NULL, excludeSNPs=NULL, tol=0.1, transpose=FALSE) {
  
  standardized <- as.logical(standardized)

  if (transpose) {
    true <- t(true)
    impute <- t(impute)
  }
  
  
  m <- ncol(true)
  n <- nrow(true)
  
  stopifnot(m == ncol(impute))
  
  if ((is.null(rownames(true)) | is.null(rownames(impute))) & nrow(impute) != n) {
    warning('Samples in true or impute are not identified by rownames, and numbers of rows do not match.\nResults may be misleading!')
  }
  
  if (is.null(rownames(true))) rownames(true) <- 1:n
  if (is.null(rownames(impute))) rownames(impute) <- 1:nrow(impute)
  
  # Manage scaling parameters, if needed.
  usermeans <-  (!is.null(p) | !is.null(center) | !is.null(scale))
  
  if (is.null(center)) center=numeric(m)
  if (is.null(scale)) {scale=numeric(m);scale[] <- 1}
  
  if (usermeans) {
    standardized <- TRUE
    
    if (!is.null(p)) {
      stopifnot(length(p)==m)
      center <- 2*p
      scale <- 2*p*(1-p)
    }
    
  } else if (standardized) {
    center <- apply(true, 2, mean, na.rm=TRUE)
    scale <- apply(true, 2, sd, na.rm=TRUE)
  }
  
  # Manage exclusion of samples (ids) or SNPs.
  ex_ids <- rep(FALSE, n)
  if (!is.null(excludeIDs)) {
    ex_ids[rownames(true) %in% excludeIDs] <- TRUE 
  }
  
  ex_snps <- rep(FALSE, m)
  if (!is.null(excludeSNPs)) {
    .is.na <- function(x) {if (is.null(x)) return(logical(0)); is.na(x)}
    if (is.logical(excludeSNPs) & sum(!.is.na(excludeSNPs)) < m) stop('`excludeSNPs` as logical must be same length as SNPs in input files and without NA\'s.')
    
    if (is.numeric(excludeSNPs)) {
      ex_snps[as.integer(excludeSNPs)] <- TRUE
    }
  }
  
  if (any(ex_ids) | any(ex_snps)) {
    true <- true[!ex_ids, !ex_snps]
    impute <- impute[rownames(true), !ex_snps]
    m <- ncol(true)
    n <- nrow(true)
    center <- center[!ex_snps]
    scale <- scale[!ex_snps]
  } else {
    impute <- impute[rownames(true), ]
  }
  
  
  ## Calculate statistics
  #if (standardized) {
    matcor <- cor(as.vector(scale(true, center, scale)), as.vector(scale(impute, center,scale)), use = 'complete.obs')
  #} else {
  #  matcor <- cor(as.vector(true), as.vector(impute), use = 'complete.obs')
  #}
  
  col.stats <- sapply(1:ncol(true), function(ii) vect.stat(true[,ii], impute[,ii], tol=tol, m=center[ii], s=scale[ii]))
  col.stats <- as.data.frame(t(col.stats))
  colnames(col.stats) <- c('cors','correct','true.na','imp.na','both.na')
  col.stats$correct.pct <- with(col.stats, correct/(n - true.na - both.na))
  col.stats <- cbind.data.frame(means=center, sds=scale, col.stats, stringsAsFactors=FALSE)
  
  row.stats <- sapply(1:nrow(true), function(ii) vect.stat(true[ii,], impute[ii,], tol=tol, m=center, s=scale))
  row.stats <- as.data.frame(t(row.stats))
  colnames(row.stats) <- c('cors','correct','true.na','imp.na','both.na')
  row.stats$correct.pct <- with(row.stats, correct/(m - true.na - both.na))
  row.stats <- cbind.data.frame(rowID=rownames(true), row.stats, stringsAsFactors=FALSE)
  
  #if (anyNA(scale)) {
  #  i <- which(is.na(scale))
  #  col.stats$true.na <- 
  #}
  
  list(matcor=matcor, snps=col.stats, animals=row.stats)
}

# private function between two columns or rows
# calculates correlation,
# values equal to each other
vect.stat <- function(x, y, tol, m, s) {
  suppressWarnings(cor <- cor((x-m)/s, (y-m)/s, use="pairwise.complete.obs"))
  correct <- sum(round(abs(x - y), 4) <= tol, na.rm=TRUE)
  x <- is.na(x)  & !is.nan(x)
  y <- is.na(y)  & !is.nan(x)
  x.na <- sum(x & !y)
  y.na <- sum(y & !x)
  both.na <- sum(x & y)
  return(c(cor=cor, correct=correct, x.na=x.na, y.na=y.na, both.na=both.na))
}

#' \code{imputation_accuracy1} and \code{imputation_accuracy3} has been replaced by \code{\link{imputation_accuracy}}.
#' The difference between the two former functions is now covered by the \code{adaptive}-argument of the latter.
#' 
#' @param truefn,imputefn Deprected, was \emph{filename} to files with true and imputed genotype matrix.
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
