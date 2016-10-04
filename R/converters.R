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
#' Use \code{numeric.format} to change from default 5 characters width and 2 decimals (\code{"5.2"}). 
#' NB! Function does not test for validity of this argument, so change with caution!
#'
#' \strong{Missing values:} Values after summing less than 0 or greater than 2 are assumed as missing and replaced with \code{naval}.
#'
#' @param phasefn Filename of input file, every two rows are for same animal.
#' @param genofn Filename of intended output.
#' @param ncol Number of SNPs in file, if \code{NULL} (default), number is estimated from file.
#' @param nrow Number of rows to maximally read from \code{phasefn}. If \code{NULL}, no limit is used.
#' @param naval Missing value, default \code{9}.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param numeric.format Character describing \code{<width>.<decimals>}. Default \code{'5.2'}. Has no effect when \code{int=TRUE}.
#' @return Number of rows written.
#' @export
convert_phases <- function(phasefn, genofn, ncol=NULL, nrow=NULL, int=TRUE, naval=9, numeric.format='5.2') {
  stopifnot(file.exists(phasefn))
  
  if (is.null(ncol)) ncol <- get_ncols(phasefn)-1
  if (is.null(nrow)) nrow <- 0
  
  subroutine <- ifelse(int, 'convert_phase_int', 'convert_phase_num')
  res <- .Fortran(subroutine, phasefn=as.character(phasefn), genofn=as.character(genofn), ncol=as.integer(ncol), nrow=as.integer(nrow), 
                  naval=as.integer(naval), numfmt=as.character(numeric.format), lennumfmt=as.integer(nchar(numeric.format)))
  res$nrow
}

#' @param phasefn Filename of input file, every two rows are for same animal.
#' @param genofn Filename of intended output.
#' @param ncol Number of columns to read, if \code{NULL} (default), number is estimated from file.
#' @param nrow Number of rows to maximally read from \code{phasefn}. If \code{NULL}, no limit is used.
#' @return Number of rows written.
#' @export
#' @rdname depcreated
phasotogeno <- function(phasefn, genofn, ncol=NULL, nrow=NULL) {
  .Deprecated('imputation_accuracy', package='Siccuracy')
  convert_phases(phasefn, genofn, ncol, nrow, int=FALSE)
}
#' @export
#' @rdname depcreated
#' @inheritParams phasotogeno
phasotogeno_int <- function(phasefn, genofn, ncol=NULL, nrow=NULL) {
  .Deprecated('imputation_accuracy', package='Siccuracy')
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
#' @param rawfn Output file from plink. Most likely \code{plink.raw} if PLINK command line argument \code{--output} is not used.
#' @param outfn Filename of new file.
#' @param NewID Integer scalar (default \code{0}) for automatically assigning new IDs. See description for more. 
#' @param ncol Number of \emph{SNP} columns in \code{rawfn}. Is total number of columns minus 6.
#' @param nlines Number of lines to process.
#' @param naval Integer scalar to use for \code{NA} values.
#' @return Data.frame with columns \code{famID}, \code{sampID}, and \code{newID}.
#' @references 
#' \itemize{
#'  \item PLINK. Purcell and Chang. \url{https://www.cog-genomics.org/plink2}
#'  \item \href{http://www.gigasciencejournal.com/content/4/1/7}{Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4.}
#' }
#' @export
convert_plinkA <- function(rawfn, outfn, newID=0, ncol=NULL, nlines=NULL, naval=9) {
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