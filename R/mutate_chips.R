# Combine chromosome files #####

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
#' @param format Character, Fortran edit descriptors for output. See \link{parseformat}.
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

#' Combines two SNP chip files and provides masking.
#'
#'
#' We distinguish between \code{HD} and \code{LD} (high and low density) SNP chips,
#' were HD chips have priority over LD.
#' 
#' \strong{\code{hdpos} and \code{ldpos}:}
#' The output file is created with \code{outcol} columns. 
#' First, rows from \code{hdfn} are added to the file with positions described in \code{hdpos}. The first column is ignored as this is the ID column.
#' If \code{hdpos=c(2,3,NA,NA,4)}, the first two columns of \code{hdfn} are used for column 2 and 3 in the output, and the last column of \code{hdfn} is used as the 4th.
#' The 3rd and 4th column in \code{hdfn} are ignored. This allows masking of genotypes.
#' After \code{hdfn} is processed for rows in \code{hdid}, the same is repeated for \code{ldfn}. There is no restriction for \code{hdfn=ldfn}.
#' The columns of the output file is prepopulated with \code{na}.
#' 
#'
#' @param hdid IDs of HD genotyped individuals. Corresponds to (subset of) first column of \code{hdfn}.
#' @param ldid IDs of LD genotyped individuals. Corresponds to (subset of) first column of \code{ldfn}.
#' @param hdpos Integer vector of where columns in \code{hdfn} are positioned in output file. Coerced to same length as columns in \code{hdfn}.
#' @param ldpos Integer vector of where columns in \code{ldfn} are positioned in output file. Coerced to same length as columns in \code{ldfn}.
#' @param hdfn Filename of HD genotypes.
#' @param ldfn Filename of LD genotypes.
#' @param fnout Filename to write merged genotypes to.
#' @param outcol Integer, number of columns in output file. When \code{NULL} (default), uses max value of \code{hdpos} and \code{ldpos}.
#' @param na Missing values.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param format Character, Fortran edit descriptors for output. See \link{parseformat}.
#' @seealso \code{\link{get_ncols}}
#' @return Number of lines written, with attribute 'stat' specifying last IO status.
#' @export
rbind_SNPs <- function(hdid,ldid, hdpos, ldpos, hdfn, ldfn, fnout, outcol=NULL, na=9, format=NULL, int=TRUE) {
  
  if (!file.exists(hdfn)) stop('`hdid` file was not found.')
  if (!file.exists(ldfn)) stop('`ldid` file was not found.')
  
  format <- parse.format(format, int)
  
  hdid <- as.integer(hdid)
  ldid <- as.integer(ldid)

  hdcols <- get_ncols(hdfn)-1
  ldcols <- get_ncols(ldfn)-1
  
  hdpos <- as.integer(hdpos[1:hdcols])
  ldpos <- as.integer(ldpos[1:ldcols])

  if (is.null(outcol)) outcol <- max(hdpos, ldpos, na.rm=TRUE)
  
  hdpos[hdpos > outcol] <- 0
  ldpos[ldpos > outcol] <- 0
  
  hdpos[is.na(hdpos)] <- 0
  ldpos[is.na(ldpos)] <- 0
  
  format <- parse.format(format, int)
  
  #subroutine rbindsnps(fnhd, fnld, fnout, hdcols, ldcols, outcols, &
  #                       nhd, hdid, nld, ldid, hdpos, ldpos,  missing, lenfmt, userfmt, asint, stat)
  res <- .Fortran('rbindsnps', 
                  fnhd=as.character(hdfn), fnld=as.character(ldfn), fnout=as.character(fnout),
                  hdcols=as.integer(hdcols), ldcols=as.integer(ldcols),
                  outcols=as.integer(outcol), 
                  nhd=as.integer(length(hdid)), hdid=as.integer(hdid),
                  nld=as.integer(length(ldid)), ldid=as.integer(ldid),
                  hdpos=as.integer(hdpos), ldpos=as.integer(ldpos),
                  missing=as.integer(na),
                  lenfmt=as.integer(nchar(format)), userfmt=as.character(format), asint=as.integer(int),
                  stat=integer(1), n=integer(1))
  structure(res$n, stat=res$stat)
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
#' @param missing Missing value.
mergeChips <- function(hdid,ldid, hdpos, ldpos, hdfn, ldfn, fnout, outpos=NULL, missing=9) {
  .Deprecated('rbind_SNPs', package='Siccuracy')
}


# Mask a single SNP chip ####

#' Masks entries of a single SNP chip file.
#' 
#' This function runs through a single file and masks specified columns of specified rows.
#' It is a simplified version of \code{\link{rbind_SNPs}}.
#' 
#' If \code{maskIDs} or \code{maskSNPs} are zero-length vectors, \code{dropIDs} and \code{dropSNPs} is still performed.
#' 
#' @param fn Input filename.
#' @param outfn Output filename.
#' @param maskIDs IDs of rows to mask.
#' @param maskSNPs Integer indices of columns in \code{fn} to mask, i.e. before considering \code{dropSNPs}.
#' @param dropIDs IDs to exclude from output.
#' @param dropSNPs Integer indices of columns in \code{fn} to exclude from output.
#' @param na Value to use for masking.
#' @param ncol Integer, number of SNP columns in \code{fn}. When \code{NULL} (default), automagically detected with \code{get_ncols(fn)-1}.
#' @param nlines Integer, number of lines to process.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param format Character, Fortran edit descriptors for output. See \link{parseformat}.
#' @export
mask_SNPs <- function(fn, outfn, maskIDs, maskSNPs, dropIDs=NULL, dropSNPs=NULL, na=9, ncol=NULL, nlines=NULL, int=TRUE, format=NULL) {
  if (!file.exists(fn)) stop('Input file was not found.')
  
  format <- parse.format(format, int)
  
  maskIDs <- as.integer(maskIDs)
  maskSNPs <- as.integer(maskSNPs)
  dropIDs <- as.integer(dropIDs)
  dropSNPs <- as.integer(dropSNPs)
  
  IDs <- get_firstcolumn(fn)
  if (!is.null(nlines)) IDs <- IDs[1:nlines]
  maskIDs <- IDs %in% maskIDs
  dropIDs <- IDs %in% dropIDs
  
  if (is.null(nlines)) nlines <- get_nlines(fn)
  if (is.null(ncol)) ncol <- get_ncols(fn)-1
  
  maskSNPs <- c(1:ncol) %in% maskSNPs
  dropSNPs <- c(1:ncol) %in% dropSNPs
  
  stopifnot(sum(dropSNPs) < ncol)
  
  res <- .Fortran('masksnps', fn=as.character(fn), outfn=as.character(outfn), ncol=as.integer(ncol), nlines=as.integer(nlines),
                  imaskIDs=as.integer(maskIDs), imaskSNPs=as.integer(maskSNPs), 
                  idropIDs=as.integer(dropIDs), idropSNPs=as.integer(dropSNPs),
                  na=as.numeric(na), userfmt=as.character(format), lenuserfmt=nchar(format), asint=as.integer(int),
                  stat=integer(1))
  res
}