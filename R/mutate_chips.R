#' Column bind files with genotype matrices
#'
#' Does row-wise concatenation of multiple genotype files, ignoring the first column.
#' @param fns Character vector of filenames to concatenate. Max 255 characters per filename due to R-Fortran.
#' @param fnout Filename of resulting file
#' @param nlines Number of lines to read from each file. When \code{NULL} (default), use all lines.
#' @param ncols Number of columns to read from each file. When \code{NULL} (default), use all columns.
#' @param skiplines Integer, number of lines to skip before outputting. \code{nlines} counts towards total number of lines, i.e. both skipped and outputted.
#' @param excludeids Integer vector of first column id's to \emph{exclude} from the output. \code{nlines} also counts towards excluded lines.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param format Character, Fortran edit descriptors for output. See \link{parse.format}.
#' @return Exit status, invisible.
#' @export
# subroutine rowconcatenate(files, fns, fnout, nlines, ncols, result)
cbind_SNPs <- function(fns, fnout, ncols=NULL, nlines=NULL, skiplines=0, excludeids=integer(0), format=NULL, int=TRUE) {
  if (!all(file.exists(fns))) {
    stop('Not all files in fns exists.')
  }
  
  if (is.null(nlines)) nlines <- get_nlines(fns[1])
  if (is.null(ncols)) ncols <- sapply(fns, get_ncols) - 1
  
  format <- parse.format(format, int)
  
  ftemp <- tempfile()
  write.table(fns, file = ftemp, row.names = FALSE, col.names=FALSE)
  
  res <- .Fortran('cbindsnpsrwrapper',
                  files=length(fns),
                  fnin=as.character(ftemp),
                  fnout=as.character(fnout),
                  nlines=as.integer(nlines),
                  ncols=as.integer(ncols),
                  skiplines=as.integer(skiplines),
                  idlength=as.integer(length(excludeids)),
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
#' Use \code{\link{cbind_SNPs}} instead.
#' 
#' @export
#' @rdname deprecated
#' @inheritParams cbind_SNPs
rowconcatenate <- function(fns, fnout, nlines=NULL, ncols=NULL, skiplines=0, excludeids=integer(0)) {
  .Deprecated('cbind_SNPs', package='Siccuracy')
  cbind_SNPs(fns=fns, fnout=fnout, nlines=nlines, ncols=ncols, skiplines=skiplines, excludeids=excludeids, int=FALSE)
}

# #' @export
# #' @rdname cbind_SNPs
# cbind_SNPS <- function(...) {cbind_SNPs(...)}
# # The difference being the last s is upper/lowercase. The latter being grammatic correct.

# Combine chips ####

#' Combines two SNP chip files
#'
#' Combines two SNP chip files into one. The two files may be the same, and 
#' the result will be 
#'
#' We distinguish between \code{HD} and \code{LD} (high and low density) SNP chips,
#' were HD chips have priority over LD.
#'
#' @param hdid IDs of HD genotyped individuals. Corresponds to (subset of) first column of \code{hdfn}.
#' @param ldid IDs of LD genotyped individuals. Corresponds to (subset of) first column of \code{ldfn}.
#' @param hdpos Integer positions of HD genotyped SNPs. Length must correspond to number of columns in \code{hdfn}, excluding ID column (first column).
#' @param ldpos Integer positions of LD genotyped SNPs. Length must correspond to number of columns in \code{hdfn}, excluding ID column (first column).
#' @param hdfn Filename of HD genotypes.
#' @param ldfn Filename of LD genotypes.
#' @param fnout Filename to write merged genotypes to.
#' @param outpos Integer vector of collective SNP positions. Default to sorted, union of \code{hdpos} and \code{ldpos}. Make it anything else and you get?
#' @param na Missing values.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param format Character, Fortran edit descriptors for output. See \link{parse.format}.
#' @seealso \code{\link{get_ncols}}
#' @export
rbind_SNPs <- function(hdid,ldid, hdpos, ldpos, hdfn, ldfn, fnout, outpos=NULL, na=9, format=NULL, int=TRUE) {
  
  if (!file.exists(hdfn)) stop('`hdid` file was not found.')
  if (!file.exists(ldfn)) stop('`ldid` file was not found.')
  
  format <- parse.format(format, int)
  
  hdid <- as.integer(hdid)
  ldid <- as.integer(ldid)
  
  ldid <- setdiff(ldid, hdid)
  
  if (is.null(outpos)) outpos <- sort(unique(c(as.integer(hdpos),as.integer(ldpos))))
  hdpos <- match(hdpos, outpos)
  ldpos <- match(ldpos, outpos)
  
  
  #subroutine rbind_SNPs(fnhd, fnld, fnout, hdcols, ldcols, outcols, &
  #                          nhd, hdid, nld, ldid, hdpos, ldpos,  missing, lenfmt, userfmt, status, asint
  res <- .Fortran('rbindsnps', 
                  fnhd=as.character(hdfn), 
                  fnld=as.character(ldfn), 
                  fnout=as.character(fnout),
                  hdcols=as.integer(length(hdpos)),
                  ldcols=as.integer(length(ldpos)), 
                  outcols=as.integer(length(outpos)),
                  nhd=as.integer(length(hdid)), 
                  hdid=as.integer(hdid),
                  nld=as.integer(length(ldid)), 
                  ldid=as.integer(ldid),
                  hdpos=as.integer(hdpos), 
                  ldpos=as.integer(ldpos),
                  missing=as.integer(na), 
                  lenfmt=nchar(format), 
                  userfmt=as.character(format),
                  status=integer(1),
                  int=as.integer(int))
  res
}

# #' @export
# #' @noRd
# rbind_SNPS <- function(...) rbind_SNPs(...)
# #' @export
# #' @noRd
# rbind_snps <- function(...) rbind_SNPs(...)

#' \code{mergeChips} has been deprecated.
#' Use \code{\link{rbind_SNPs}} instead.
#' 
#' @export
#' @rdname deprecated
#' @inheritParams rbind_SNPs
mergeChips <- function(hdid,ldid, hdpos, ldpos, hdfn, ldfn, fnout, outpos=NULL, missing=9) {
  .Deprecated('rbind_SNPs', package='Siccuracy')
}