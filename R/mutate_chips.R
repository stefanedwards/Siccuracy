#' Column bind files with genotype matrices
#'
#' Does row-wise concatenation of multiple genotype files, ignoring the first column.
#' @param fns Character vector of filenames to concatenate. Max 255 characters per filename due to R-Fortran.
#' @param fnout Filename of resulting file
#' @param nlines Number of lines to read from each file. When \code{NULL} (default), use all lines.
#' @param ncols Number of columns to read from each file. When \code{NULL} (default), use all columns.
#' @param skiplines Integer, number of lines to skip before outputting. \code{nlines} counts towards total number of lines, i.e. both skipped and outputted.
#' @param excludeids Integer vector of first column id's to \emph{exclude} from the output. \code{nlines} also counts towards excluded lines.
#' @param format Character of Fortran formats. Defaults to \code{F5.2}; use \code{I2} for integers, or numeric() or integer() for the former two.
#' @return Exit status, invisible.
#' @export
# subroutine rowconcatenate(files, fns, fnout, nlines, ncols, result)
cbind_SNPs <- function(fns, fnout, nlines=NULL, ncols=NULL, skiplines=0, excludeids=integer(0), format='5.2', int=TRUE) {
  if (!all(file.exists(fns))) {
    stop('Not all files in fns exists.')
  }
  
  if (is.null(nlines)) nlines <- get_nlines(fns[1])
  if (is.null(ncols)) ncols <- sapply(fns, get_ncols) - 1
  if (length(ncols) < length(fns)) rep(ncols, length.out=length(fns))
  
  if (is.integer(format)) {
    format <- 'I2'
  } else if (is.numeric(format)) {
    format <- 'F5.2'
  }
  format <- gsub('^[[:alpha:]]','',format)[1]
  
  
  ftemp <- tempfile()
  write.table(fns, file = ftemp, row.names = FALSE, col.names=FALSE)
  
  res <- .Fortran('cbindsnpsrwrapper',
                  files=length(fns),
                  fnin=as.character(ftemp),
                  fnout=as.character(fnout),
                  nlines=as.integer(nlines),
                  ncols=as.integer(ncols),
                  skiplines=as.integer(skiplines),
                  idlength=length(excludeids),
                  excludeids=as.integer(excludeids),
                  result=integer(1),
                  lenfmt=nchar(format),
                  uesrfmt=as.character(format),
                  int=as.integer(int),
                  PACKAGE='Siccuracy')
  file.remove(ftemp)
  invisible(res$result==0)
}

#' \code{rowconcatenate} has been deprecated.
#' 
#' @export
#' @rdname cbind_SNPs
#' @inheritParams cbind_SNPs
rowconcatenate <- function(fns, fnout, nlines=NULL, ncols=NULL, skiplines=0, excludeids=integer(0)) {
  .Deprecated('cbind_SNPs', package='Siccuracy')
  cbind_SNPs(fns=fns, fnout=fnout, nlines=nlines, ncols=ncols, skiplines=skiplines, excludeids=excludeids)
}