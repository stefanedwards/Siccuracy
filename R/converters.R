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
  if (is.null(ncol)) ncol <- get_ncols(rawfn) - 6
  
  firstcols <- get_firstcolumn(rawfn, class=list('character','character'), col.names=c('famID','sampID'))
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
  
  newID <- .newID[1:nlines,]
  
  #subroutine convert_plinkA(rawfn, outputfn, newID, ncol, nrow, naval, stat) 
  res <- .Fortran('convertplinka', rawfn=as.character(rawfn), outputfn=as.character(outfn), newID=as.integer(newID$newID),
                  ncol=as.integer(ncol), nrow=as.integer(nlines), naval=as.integer(naval),
                  stat=integer(1), PACKAGE='Siccuracy', NAOK=TRUE)
  
  newID
}