# Combine chromosome files #####

#' Column bind files with genotype matrices
#'
#' Does row-wise concatenation of multiple genotype files, 
#' ignoring the first column in subsequent files.
#' \strong{NB!} Assumes rows are ordered identically in all files.
#' 
#' @param fns Character vector of filenames to concatenate. Max 255 characters per filename due to R-Fortran.
#' @param fnout Filename of resulting file
#' @param nlines Number of lines to read from each file. When \code{NULL} (default), use all lines.
#' @param ncols Number of columns to read from each file. When \code{NULL} (default), use all columns.
#' @param skiplines Integer, number of lines to skip before outputting. \code{nlines} counts towards total number of lines, i.e. both skipped and outputted.
#' @param excludeids Integer vector of first column id's to \emph{exclude} from the output. \code{nlines} also counts towards excluded lines.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param format Character, Fortran edit descriptors for output. See \link{parseformat}.
#' @return Exit status, invisible.
#' @rdname cbind_snp_files
#' @seealso \code{\link{rbind_snp_files}}
#' @export
cbind_snp_files <- function(fns, fnout, ncols=NULL, nlines=NULL, skiplines=0, excludeids=integer(0), format=NULL, int=TRUE) {
  if (!all(file.exists(fns))) {
    stop('Not all files in fns exists.')
  }
  
  if (is.null(nlines)) nlines <- get_nlines(fns[1])
  if (is.null(ncols)) ncols <- sapply(fns, get_ncols) - 1
  
  format <- parse.format(format, int)
  
  ftemp <- tempfile()
  utils::write.table(fns, file = ftemp, row.names = FALSE, col.names=FALSE, quote=FALSE)
  
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

#' @export
#' @rdname cbind_snp_files
#' @inheritDotParams cbind_snp_files
cbind_SNP_files <- function(...) { cbind_snp_files(...) }


#' Combines rows of two SNP chip files and provides masking.
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
#' @seealso \code{\link{get_ncols}}, \code{\link{cbind_snp_files}}, 
#'           or \code{\link{mask_snp_file}} for more flexible masking (albeit on a single file).
#' @return Number of lines written, with attribute 'stat' specifying last IO status.
#' @rdname rbind_snp_files
#' @export
rbind_snp_files <- function(hdid,ldid, hdpos, ldpos, hdfn, ldfn, fnout, outcol=NULL, na=9, format=NULL, int=TRUE) {
  
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

#' @rdname rbind_snp_files
#' @inheritDotParams rbind_snp_files
#' @export
rbind_SNP_files <- function(...) {rbind_snp_files(...)}

# Mask a single SNP chip ####

#' Masks entries of a single SNP chip file.
#' 
#' This function runs through a \emph{single} file and masks specified columns of specified rows.
#' It is a simplified version of \code{\link{rbind_snp_files}}.
#' 
#' The \code{masking} argument controls which samples / loci are masked as missing with \code{na}.
#' It accepts three different objects: 
#' 1) an integer vector of loci to mask, 
#' 2) a list with first element containing integer vector of IDs for samples to mask and second element an integer vector of loci to mask for those samples,
#' or 
#' a list of 2).
#' For 1), the given loci are masked for \emph{all} samples. For 2) it is limited to those IDs given.
#' By default, the integer vector of loci to mask are mapped to the input file. To map them to output file, set \code{snpsinnew=TRUE}.
#' 
#' 
#' @param fn Input filename.
#' @param outfn Output filename.
#' @param snps Vector of column indicies to use in output file. Defaults to all.
#' @param masking List that specifies how to mask. Positions given here refer to columns in \code{fn}. See details.
#' @param snpsinnew Logical, when \code{TRUE}, positions in \code{masking} are mapped to the output file instead of input file.
#' @param dropIDs IDs to exclude from output.
#' @param na Value to use for masking.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param format Character, Fortran edit descriptors for output. See \link{parseformat}.
#' @export
#' @seealso \code{\link{rbind_snp_files}}
#' @rdname mask_snp_file
#' @return Invisible list of vectors sent to Fortran subtroutine.
#' @examples
#' 
#' SNPs <- Siccuracy:::make.true(9, 12)
#' snpfile <- tempfile()
#' write.snps(SNPs, snpfile)
#' 
#' masking = list(
#' list(1:3, c(1,3,5,7,9,11)),
#' list(4:6, c(2,4,6,8,10,12)),
#' list(7:9, c(1,3,9,10,12))
#' )
#' 
#' fn <- tempfile()
#' res <- mask_SNPs(snpfile, fn, masking=masking, na=9)
#' m <- read.snps(fn)
#' 
mask_snp_file <- function(fn, outfn, masking=NULL, snps=NULL, snpsinnew=FALSE, dropIDs=NULL, na=9, int=TRUE, format=NULL) {
  if (!file.exists(fn)) stop('Input file was not found.')
  
  format <- parse.format(format, int)
  dropIDs <- as.integer(dropIDs)
  
  IDs <- get_firstcolumn(fn)
  n <- length(IDs)
  m <- get_ncols(fn) - 1
  
  if (is.null(snps)) {
    snps <- 1:m
  } else {
    snps <- snps[snps >= 0 & snps <= m]
  }
  
  # Handling masking
  maskIDs <- integer(n)
  maskstart <- integer(0)
  maskend <- integer(0)
  maskSNPs <- integer(0)
  if (is.list(masking)) {
    if (length(masking) == 2) {
      if (!is.list(masking[[1]]) & !is.list(masking[[2]])) {
        masking <- list(masking)
      }
    }
    for (i in 1:length(masking)) {
      maskIDs[match(masking[[i]][[1]], IDs)] <- i
      masking[[i]][[2]] <- masking[[i]][[2]][masking[[i]][[2]] <= m & masking[[i]][[2]] > 0]
    }
    maskSNPs <- unlist(sapply(masking, function(x) x[[2]], simplify = TRUE))
    masklengths <- unlist(sapply(masking, function(x) length(x[[2]]), simplify = TRUE))
    maskend <- cumsum(masklengths)
    maskstart <- c(1, maskend[-length(maskend)]+1)

  } else if (is.vector(masking)) {
    # Ordinary vector with columns to 
    masking <- as.integer(masking)
    masking <- masking[masking <= m & masking > 0]
    maskIDs[] <- 1
    maskstart <- 1
    maskSNPs <- masking
    maskend <- length(maskSNPs)
  } else  if (is.null(masking)) {
    # Do nothing
  } else {
    stop('`mask_SNPs` does not understand how to deal with the given `masking`, as it is neither a list or integer vector.')
  }
  
  # Convert indices to input columns.
  if (snpsinnew) {
    maskSNPs <- snps[maskSNPs]
  }

  dropIDs <- IDs %in% dropIDs
  
  #subroutine masksnps(fn, outfn, ncols, nlines, na, userfmt, lenuserfmt, asint, stat
  res <- .Fortran('masksnps2', fn=as.character(fn), outfn=as.character(outfn), ncol=as.integer(m), nlines=as.integer(n),
                  na=as.numeric(na), userfmt=as.character(format), lenuserfmt=nchar(format), asint=as.integer(int), stat=integer(1),
                  dropIDs=as.integer(dropIDs),
                  snpslength=as.integer(length(snps)), snps=as.integer(snps), 
                  maskIDs=as.integer(maskIDs), maps=as.integer(length(maskstart)), maskstart=as.integer(maskstart), maskend=as.integer(maskend),
                  masklength=as.integer(length(maskSNPs)), maskSNPs=as.integer(maskSNPs))
  invisible(res)
}
#' @rdname mask_snp_file
#' @inheritParams mask_snp_file
#' @export
mask_SNPs <- function(fn, outfn, masking=NULL, snps=NULL, snpsinnew=FALSE, dropIDs=NULL, na=9, int=TRUE, format=NULL) {
  mask_snp_file(fn=fn, outfn=outfn, masking=masking, snps=snps, snpsinnew=snpsinnew, dropIDs=dropIDs, na=na, int=int, format=format)
}
