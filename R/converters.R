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
#' Use \code{format} to change from default 5 characters width and 2 decimals (\code{"5.2"}). 
#' NB! Function does not test for validity of this argument, so change with caution!
#'
#' \strong{Missing values:} Values after summing less than 0 or greater than 2 are assumed as missing and replaced with \code{na}.
#'
#' @param phasefn Filename of input file, every two rows are for same animal.
#' @param genofn Filename of intended output.
#' @param ncol Integer, number of SNP columns in files. When \code{NULL}, automagically detected with \code{get_ncols(phasefn)-1}.
#' @param nlines Integer, maximum number of pairs of lines to convert.
#' @param na Missing value.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param format Character, Fortran edit descriptors for output. See \link{parse.format}.
#' @return Number of rows written.
#' @export
convert_phases <- function(phasefn, genofn, ncol=NULL, nrow=NULL, na=9, int=TRUE, format=NULL) {
  stopifnot(file.exists(phasefn))
  
  if (is.null(ncol)) ncol <- get_ncols(phasefn)-1
  if (is.null(nrow)) nrow <- 0
  
  format <- parse.format(format, int)
  
  # subroutine convert_phase(phasefn, genofn, ncol, nrow, na, int, lenfmt, userfmt)
  res <- .Fortran('convert_phase', phasefn=as.character(phasefn), genofn=as.character(genofn), ncol=as.integer(ncol), nrow=as.integer(nrow), 
                  na=as.integer(na), lenfmt=as.integer(nchar(format)), userfmt=as.character(userfmt))
  res$nrow
}

#' Deprecated functions
#' 
#' \code{phasotogeno} and \code{phasotogeno_int} has been replaced by \code{\link{convert_phases}}. 
#' The difference between the two former functions is covered by the \code{int} argument in the new function.
#' 
#' @param phasefn Filename of input file, every two rows are for same animal.
#' @param genofn Filename of intended output.
#' @param ncol Number of columns to read, if \code{NULL} (default), number is estimated from file.
#' @param nrow Number of rows to maximally read from \code{phasefn}. If \code{NULL}, no limit is used.
#' @return Number of rows written.
#' @export
#' @rdname depcreated
phasotogeno <- function(phasefn, genofn, ncol=NULL, nrow=NULL) {
  .Deprecated('convert_phases', package='Siccuracy')
  convert_phases(phasefn, genofn, ncol, nrow, int=FALSE)
}
#' @export
#' @rdname depcreated
#' @inheritParams phasotogeno
phasotogeno_int <- function(phasefn, genofn, ncol=NULL, nrow=NULL) {
  .Deprecated('convert_phases', package='Siccuracy')
  convert_phases(phasefn, genofn, ncol, nrow, int=TRUE)
}



# Convert plink A format #########################

#' Convert PLINK recoded A to SNP file.
#' 
#' Facilitates converting a PLINK binary file to simplified SNP file format.
#' Requires using PLINK to recode it to the \code{A} format by using command line \code{plink -bfile <file stem> --recode A}.
#' This function then swiftly strips of first 6 columns (family ID, sample ID, paternal ID, maternal ID, sex, phenotypic record) 
#' and inserts an integer-based ID column. \code{NA}'s outputted from PLINK are replaced with \code{naval} argument.
#' 
#' The new integer IDs can be supplied. If not, they will be made for you.
#' \code{newID} may be an integer vector and will be used as is.
#' If data.frame with columns \code{famID}, \code{sampID}, and \code{newID}, they will be reordered to match input file.
#' 
#' @param rawfn Plink output filename. Most likely \code{plink.raw} if PLINK command line argument \code{--output} is not used.
#' @param outfn Filename of new file.
#' @param newID Integer scalar (default \code{0}) for automatically assigning new IDs. See description for more. 
#' @param ncol Integer,number of SNP columns in \code{rawfn}. When \code{NULL}, automagically detected with \code{get_ncols(rawfn)-6}.
#' @param nlines Number of lines to process.
#' @param na Missing value, 
#' @return Data.frame with columns \code{famID}, \code{sampID}, and \code{newID}.
#' @references 
#' \itemize{
#'  \item PLINK. Purcell and Chang. \url{https://www.cog-genomics.org/plink2}
#'  \item \href{http://www.gigasciencejournal.com/content/4/1/7}{Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4.}
#' }
#' @export
convert_plinkA <- function(rawfn, outfn, newID=0, ncol=NULL, nlines=NULL, na=9) {
  stopifnot(file.exists(rawfn))
  
  if (is.data.frame(newID)) stopifnot(all(c('famID','sampID','newID') %in% names(newID)))
  #if (is.null(ncol)) ncol <- get_ncols(rawfn) - 6
  
  firstline <- scan(rawfn, what=character(), nlines=1, quiet=TRUE)
  if (is.null(ncol)) ncol <- length(firstline) - 6
  header <- all(firstline[1:3] == c('FID','IID','PAT'))
  
  firstcols <- get_firstcolumn(rawfn, class=list('character','character'), col.names=c('famID','sampID'), header=header)
  if (is.null(nlines)) {
    if (length(newID) == 1) {
      nlines <- nrow(firstcols)
    } else if (is.data.frame(newID)) {
      nlines = nrow(newID)
    } else {
      nlines = length(newID)
    }
  }
  
  if (length(newID) == 1) {
    .newID <- firstcols
    .newID$newID <- 1:nrow(.newID) + newID
  } else if (is.atomic(newID)) {
    .newID <- cbind(firstcols[1:nlines,], newID[1:nlines])
  } else {
    newID <- do.call(data.frame, append(lapply(newID, as.character), list(stringsAsFactors=FALSE)))
    .newID <- merge(firstcols, newID, by=c('famID','sampID'), sort=FALSE, all.x=TRUE, all.y=FALSE, stringsAsFactors=FALSE)
  }
  
  nlines <- min(nlines, nrow(.newID))
  newID <- .newID[1:nlines,]
  newID$newID <- as.integer(newID$newID)
  
  #subroutine convert_plinkA(rawfn, outputfn, newID, ncol, nrow, naval, stat) 
  res <- .Fortran('convertplinka', rawfn=as.character(rawfn), outputfn=as.character(outfn), newID=as.integer(newID$newID),
                  ncol=as.integer(ncol), nrow=as.integer(nlines), naval=as.integer(na), header=as.integer(header),
                  stat=integer(1), PACKAGE='Siccuracy', NAOK=TRUE)
  
  newID
}


# Convert PLINK binary files ########
# Our understanding of the different filetypes.
# The binary bim/bed/ped format lists for each loci (in bm-file) the major and minor allele.
# The corresponding bed file counts the *minor* allele.
#
# Using `PLINK --recode A` option probably counts the *major* allele; it is llisted in the header. 

#' Converts PLINK binary format to flat format
#' 
#' @details 
#' The new integer IDs can be supplied. If not, they will be made for you.
#' \code{newID} may be an integer vector and will be used as is.
#' If data.frame with columns \code{famID}, \code{sampID}, and \code{newID}, they will be reordered to match input file.
#' 
#' \strong{\code{method}} \emph{simple} stores entire genotype matrix \emph{in memory}, as PLINK binary files are stored in locus-major mode, 
#' i.e. first \eqn{m} bits store first locus for all \eqn{n} animals.
#' Since we are interested in writing out all \eqn{m} loci for each animal, for efficiency we need to read the entire file.
#' \emph{As an alternate} ...
#' 
#' @param bfile Filename of PLINK binary files, i.e. without extension.
#' @param outfn Filename of new file.
#' @param na Missing value.
#' @param newID Integer scalar (default \code{0}) for automatically assigning new IDs. See description for more. 
#' @param nlines Number of lines to process.
#' @param fam If binary files have different stems, specify each of them with \code{fam}, \code{bim}, \code{bed}, and set \code{bfile=NULL}.
#' @param bim See \code{fam}.
#' @param bed See \code{fam}.
#' @param countminor Logical: Should the output count minor allele (default), or major allele as \code{plink --recode A}.
#' @param maf Numeric, restrict SNPs to SNPs with this frequency. 
#' @param extract Extract only these SNPs.
#' @param exclude Do not extract these SNPs.
#' @param keep Keep only these samples.
#' @param Remove Removes these samples from output.
#' @param method Character, which of following methods to use: \code{simple} (store entire matrix in memeory) or ??
#' @export
#' @references 
#' \itemize{
#'  \item PLINK v. 1.07 BED file format: \url{http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml}
#'  \item Shaun Purvell and Christopher Chang. \emph{PLINK v. 1.90} \url{https://www.cog-genomics.org/plink2}
#'  \item Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) \href{http://www.gigasciencejournal.com/content/4/1/7}{Second-generation PLINK: rising to the challenge of larger and richer datasets.} \emph{GigaScience}, 4.
#' }
#
convert_plink <- function(bfile, outfn, na=9, newID=0, nlines=NULL, fam=NULL, bim=NULL, bed=NULL, countminor=TRUE, maf=0.0, extract=NULL, exclude=NULL, keep=NULL, remove=NULL, method='simple') {
  
  # Get filenames
  if (!(is.null(bfile) | is.na(bfile))) {
    bfile <- sub_ext(bfile, c('fam','bim','bed'))
    fam <- bfile[1]
    bim <- bfile[2]
    bed <- bfile[3]
  }
  stopifnot(all(file.exists(fam, bim, bed)))
  
  # Make up new IDs, as in convert_plinkA:
  if (is.data.frame(newID)) stopifnot(all(c('famID','sampID','newID') %in% names(newID)))
  #if (is.null(ncol)) ncol <- get_ncols(rawfn) - 6
  
  firstline <- scan(fam, what=character(), nlines=1, quiet=TRUE)
  if (is.null(ncol)) ncol <- length(firstline) - 6
  header <- all(firstline[1:3] == c('FID','IID','PAT'))
  
  firstcols <- get_firstcolumn(fam, class=list('character','character'), col.names=c('famID','sampID'), header=header)
  if (is.null(nlines)) {
    if (length(newID) == 1) {
      nlines <- nrow(firstcols)
    } else if (is.data.frame(newID)) {
      nlines = nrow(newID)
    } else {
      nlines = length(newID)
    }
  }
  
  if (length(newID) == 1) {
    .newID <- firstcols
    .newID$newID <- 1:nrow(.newID) + newID
  } else if (is.atomic(newID)) {
    .newID <- cbind(firstcols[1:nlines,], newID[1:nlines])
  } else {
    newID <- do.call(data.frame, append(lapply(newID, as.character), list(stringsAsFactors=FALSE)))
    .newID <- merge(firstcols, newID, by=c('famID','sampID'), sort=FALSE, all.x=TRUE, all.y=FALSE, stringsAsFactors=FALSE)
  }

  nlines <- min(nlines, nrow(.newID))
  newID <- .newID[1:nlines,]
  newID$newID <- as.integer(newID$newID)
  
  # Number of columns:
  ncol <- get_nlines(bim)
  
  
  # Detect subroutine:
  use.method <- pmatch(method, c('simple','lowmem'))
  
  if (use.method == 1) {
    #subroutine readplinksimple(bed, fnout, ncol, nlines, na, newID, status)
    res <- .Fortran('readplinksimple', bed=as.character(bed), fnout=as.character(outfn[1]), 
                    ncol=as.integer(ncol), nlines=as.integer(nlines), na=as.integer(na), newID=as.integer(newID$newID), minor=as.integer(countminor), 
                    maf=as.numeric(maf), extract=as.integer(rep(1, ncol)), keep=integer(nlines), 
                    status=as.integer(0))
    ret <- newID
  } else if (use.method == 2) {
    # The not-so-simple complex way to do stuff. Boy, do we get to have fun now!
    
    
    
  }  
  #res$newID=newID
  #res
}