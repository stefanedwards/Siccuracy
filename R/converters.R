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
#' @param fn Filename of input file, every two rows are for same animal.
#' @param outfn Filename of intended output.
#' @param ncol Integer, number of SNP columns in files. When \code{NULL}, automagically detected with \code{get_ncols(phasefn)-1}.
#' @param nlines Integer, maximum number of pairs of lines to convert.
#' @param na Missing value.
#' @param int Logical (default \code{TRUE}), read and write integers.
#' @param format Character, Fortran edit descriptors for output. See \link{parseformat}.
#' @return Number of rows written.
#' @export
convert_phases <- function(fn, outfn, ncol=NULL, nlines=NULL, na=9, int=TRUE, format=NULL) {
  stopifnot(file.exists(fn))
  
  if (is.null(ncol)) ncol <- get_ncols(fn)-1
  if (is.null(nlines)) nlines <- 0
  
  format <- parse.format(format, int)
  
  # subroutine convert_phase(phasefn, genofn, ncol, nrow, na, int, lenfmt, userfmt)
  res <- .Fortran('convert_phase', phasefn=as.character(fn), genofn=as.character(outfn), ncol=as.integer(ncol), nrow=as.integer(nlines), 
                  na=as.integer(na), int=as.integer(int), lenfmt=as.integer(nchar(format)), userfmt=as.character(format))
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
#' @rdname deprecated
phasotogeno <- function(phasefn, genofn, ncol=NULL, nrow=NULL) {
  .Deprecated('convert_phases', package='Siccuracy')
  convert_phases(phasefn, genofn, ncol, nrow, int=FALSE)
}
#' @export
#' @rdname deprecated
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
#' @seealso 
#' \code{\link{convert_plink}} is a direct conversion that does not rely on PLINK.
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
#' \emph{lowmem} breaks the loci into smaller chunks (e.g. by chromosome), writes each chunk to a file, and merges them back as with \code{\link{cbind_SNPs}}.
#' \emph{dryrun} does not call the Fortran subroutine, but returns the treated arguments that would have been sent to the subroutine.
#'
#' For \code{method='lowmem'} use argument \code{fragment} to indicate how the loci are subdivided. 
#' When \code{fragment='chr'} (case unsensitive), loci are split according to 1st column of .bim file.
#' If \code{fragment} is a scalar integer, loci are split into this number of blocks.
#' If an integer vector of same length as \code{ncol}, it directly specifies which block a locus is sent to. \code{max(fragment)} specifies the number of blocks.
#' 
#'   
#' @section Filtering loci or samples:
#' Filters on loci or samples can be employed in a number of ways; filtering on loci and samples are handled independently.
#' Inclusion criteria (\code{extract} and \code{keep}) reduces the output to only those loci or samples that pass the criteria.
#' Exclusion criteria (\code{exclude} and \code{remove}) are applied \emph{after} inclusion criteria, and \emph{reduces} the output further.
#'
#' \code{extract} and \code{exclude} can be any combination of:
#' \describe{
#'  \item{Logical}{Vector of same length as loci in input file.}
#'  \item{Integer or numeric}{Indicates positional which loci to include or exclude. Numeric vectors are coerced to integer vectors.}
#'  \item{Character}{Matched against probe IDs, i.e. 2nd column of .bim file.}
#' }
#' For restricting the output to certain chromosomes, use \code{extract_chr}. The output is the intersect of \code{extract} and \code{exctract_chr}.
#' 
#' \code{keep} and \code{remove} are as \code{exctract} and \code{exclude} above, can be a combination of, and can additionally be:
#' \describe{
#'  \item{Character}{Matched against both famID or sampID, i.e. 1st and 2nd column of .fam file.}
#'  \item{List with named elements \code{famID} and/or \code{sampID}}{The named elements are matched against, respectively, the 1st and 2nd column of the .fam file.}
#' }
#' 
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
#' @param chr Vector of chromosomes to limit output to.
#' @param extract Extract only these SNPs, see Details.
#' @param exclude Do not extract these SNPs, see Details.
#' @param extract_chr Extract only these chromosomes, see Details.
#' @param keep Keep only these samples, see Details.
#' @param remove Removes these samples from output, see Details.
#' @param fragments \code{"chr"} or integer vector. Only used when \code{method='lowmem'}.
#' @param remerge Logical, whether to re-merge fragmented blocks. Only used when \code{method='lowmem'}.
#' @param fragmentfns Character vector or function for producing filenames.
#' @param method Character, which of following methods to use: \code{simple}, \code{lowmem}, or \code{drymem}. See Details.
#' @export
#' @references 
#' \itemize{
#'  \item PLINK v. 1.07 BED file format: \url{http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml}
#'  \item Shaun Purvell and Christopher Chang. \emph{PLINK v. 1.90} \url{https://www.cog-genomics.org/plink2}
#'  \item Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) \href{http://www.gigasciencejournal.com/content/4/1/7}{Second-generation PLINK: rising to the challenge of larger and richer datasets.} \emph{GigaScience}, 4.
#' }
#' @seealso 
#' \code{convert_plink} is a direct conversion that does not rely on PLINK.
#' See the alternate \code{\link{convert_plinkA}} which re-formats the output from \code{plink --recode A}.
convert_plink <- function(bfile, outfn, na=9, newID=0, nlines=NULL, fam=NULL, bim=NULL, bed=NULL, countminor=TRUE, maf=0.0, chr=NULL, extract=NULL, exclude=NULL, extract_chr=NULL, keep=NULL, remove=NULL, method='simple', fragments="chr", remerge=TRUE, fragmentfns=NULL) {
  
  # Get filenames
  if (!(is.null(bfile) | is.na(bfile))) {
    bfile <- sub_ext(bfile, c('fam','bim','bed'))
    fam <- bfile[1]
    bim <- bfile[2]
    bed <- bfile[3]
  }
  stopifnot(all(file.exists(fam, bim, bed)))

  # Detect subroutine:
  use.method <- pmatch(method, c('simple','lowmem','dryrun'))

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
    .newID <- firstcols
    .newID$newID <- 0
    .newID$newID[1:length(newID)] <- newID
  } else {
    newID <- do.call(data.frame, append(lapply(newID[,c('famID','sampID','newID')], as.character), list(stringsAsFactors=FALSE)))
    .newID <- merge(firstcols, newID, by=c('famID','sampID'), sort=FALSE, all.x=TRUE, all.y=FALSE, stringsAsFactors=FALSE)
  }

  nlines <- nrow(.newID)
  newID <- .newID[1:nlines,]
  newID$newID <- as.integer(newID$newID)
  
  
  # Handle samples
  .keep <- newID$newID != 0
    
  if (!is.null(keep) | !is.null(remove)) {
    if (is.logical(keep)) {
      if (sum(!is.na(keep)) < nlines) stop('`keep` as logical must be at least same length as samples in input file, without NA\'s, and as `newID` (if not scalar).')
    }
    if (is.logical(remove)) {
      if (sum(!is.na(remove)) < nlines) stop('`keep` as logical must be at least same length as samples in input file, without NA\'s, and as `newID` (if not scalar).')
    }

    if (is.numeric(keep)) keep <- as.integer(keep)
    if (is.numeric(remove)) remove <- as.integer(remove)
   
    if (!is.logical(keep) & !is.integer(keep) & !is.list(keep)) keep <- as.character(keep)
    if (!is.logical(remove) & !is.integer(remove) & !is.list(remove)) remove <- as.character(remove)
    
    if (is.logical(keep)) {
      .keep <- keep
    } else if (length(keep) > 0) {
      .keep[] <- FALSE
      if (is.integer(keep)) {
        .keep[keep] <- TRUE
      } else if (is.character(keep)) {
        .keep[newID$famID %in% keep] <- TRUE
        .keep[newID$sampID %in% keep] <- TRUE
      } else if (is.list(keep)) {
        .keep[newID$famID %in% keep$famID] <- TRUE
        .keep[newID$sampID %in% keep$sampID] <- TRUE
      }
    }
    
    if (is.logical(remove)) {
      .keep <- .keep & !remove
    } else if (length(remove) > 0) {
      if (is.integer(remove)) {
        .keep[remove] <- FALSE
      } else if (is.character(remove)) {
        .keep[newID$famID %in% remove] <- FALSE
        .keep[newID$sampID %in% remove] <- FALSE
      } else if (is.list(remove)) {
        .keep[newID$famID %in% remove$famID] <- FALSE
        .keep[newID$sampID %in% remove$sampID] <- FALSE
      }      
    }
  }
  newID$keep <- .keep
  
  # Handle SNPS
  # Number of columns:
  ncol <- get_nlines(bim)
  .extract <- rep(TRUE, ncol)
  snps <- NULL
  
  if (!is.null(chr)) {
    snps <- get_firstcolumn(bim, class=c('character','character'), col.names=c('chr','rs'), stringsAsFactors=FALSE)
    .extract <- snps$chr %in% as.character(chr)
  }
  
  # Decide upon exclude/include
  if (!is.null(extract) | !is.null(exclude) | !is.null(extract_chr)) {
    .is.na <- function(x) {if (is.null(x)) return(logical(0)); is.na(x)}
    if (is.logical(extract) & sum(!.is.na(extract)) < ncol) stop('`extract` as logical must be same length as SNPs in input file and without NA\'s.')
    if (is.logical(exclude) & sum(!.is.na(exclude)) < ncol) stop('`exclude` as logical must be same length as SNPs in input file and without NA\'s.')
    
    if (is.numeric(extract)) extract <- as.integer(extract)
    if (is.numeric(exclude)) exclude <- as.integer(exclude)
    
    if (!is.logical(extract) & !is.integer(extract)) extract <- as.character(extract)
    if (!is.logical(exclude) & !is.integer(exclude)) exclude <- as.character(exclude)

    if ((is.character(extract) | is.character(exclude) | !is.null(extract_chr)) & is.null(snps)) {
      snps <- get_firstcolumn(bim, class=c('character','character'), col.names=c('chr','rs'), stringsAsFactors=FALSE)
    }    
  
    if (is.logical(extract)) {
      .extract <- extract
    } else if (length(extract) > 0) {
      .extract[] <- FALSE
      if (is.integer(extract)) {
        .extract[extract] <- TRUE
      } else {
        .extract[snps$rs %in% extract] <- TRUE
      }
    }
    
    if (!is.null(extract_chr)) {
      extract_chr <- as.character(extract_chr)
      .extract <- .extract & snps$chr %in% extract_chr
    }
    
    if (is.logical(exclude)) {
      .extract <- .extract & !exclude
    } else if (length(exclude) > 0) {
      if (is.integer(exclude)) {
        .extract[exclude] <- FALSE
      } else {
        .extract[snps$rs %in% exclude] <- FALSE
      }
    }
  }
  
  
  ##  Fragments
  if (use.method != 1) {
    if (is.character(fragments) & length(fragments) == 1 & tolower(fragments) == 'chr') {
      if (is.null(snps)) snps <- get_firstcolumn(bim, class=c('character','character'), col.names=c('chr','rs'), stringsAsFactors=FALSE)
      fragments <- cumsum(!duplicated(snps$chr))
    }
    if (!is.integer(fragments)) fragments <- as.integer(fragments)
    if (length(fragments) == 1) {
      fragments <- rep(1:fragments, each=ceiling(ncol/fragments))[1:ncol]
    }
    if (length(fragments) != ncol) stop('Argument `fragments` should be scalar or same length as number of loci.')
    
    
    if (is.function(fragmentfns)) {
      n <- length(formals(fragmentfns))
      if (n == 0) {
        tmpfiles <- replicate(max(fragments), fragmentfns(), simplify=TRUE)
      } else if (n == 1) {
        tmpfiles <- sapply(1:max(fragments), fragmentfns)
      } else {
        tmpfiles <- sapply( 1:max(fragments), fragmentfns, max(fragments))
      }
    } else {
      tmpfiles <- as.character(fragmentfns)
    }
    tmpfiles <- c(tmpfiles, as.character(replicate(2*max(fragments)-length(tmpfiles), tempfile())))
    
  }  
  
  if (use.method == 1) {
    #subroutine readplinksimple(bed, fnout, ncol, nlines, na, newID, minor, maf, extract, keep, status)
    res <- .Fortran('readplinksimple', bed=as.character(bed), fnout=as.character(outfn), 
                    ncol=as.integer(ncol), nlines=as.integer(nlines), na=as.integer(na), newID=as.integer(newID$newID), minor=as.integer(countminor), 
                    maf=as.numeric(maf), extract=as.integer(.extract), keep=as.integer(.keep), 
                    status=as.integer(0))
  } else if (use.method == 2) {
    # The not-so-simple complex way to do stuff. Boy, do we get to have fun now!
    tmpfile <- tempfile()
    writeLines(tmpfiles, tmpfile)
    
    #subroutine convertplinkrwrapper(listfn, n, remerge, fragments, &
    #                       bed, fnout, ncol, nlines, na, newID, minor, maf, extract, keep, status)
    res <- .Fortran('convertplinkrwrapper', listfn=as.character(tmpfile), n=as.integer(length(tmpfiles)), remerge=as.integer(remerge), fragments=as.integer(fragments),
                    bed=as.character(bed), fnout=as.character(outfn), 
                    ncol=as.integer(ncol), nlines=as.integer(nlines), na=as.integer(na), newID=as.integer(newID$newID),
                    minor=as.integer(countminor), maf=as.numeric(maf), extract=as.integer(.extract), keep=as.integer(.keep),
                    status=as.integer(0))
    
    res$listfn <- NULL
    res$n <- NULL
    res$fragmentfns <- tmpfiles
  } else if (use.method == 3) {
    res <- list(bed=as.character(bed), fnout=as.character(outfn), 
                ncol=as.integer(ncol), nlines=as.integer(nlines), na=as.integer(na), newID=as.integer(newID$newID), minor=as.integer(countminor), 
                maf=as.numeric(maf), extract=as.integer(.extract), keep=as.integer(.keep),
                fragments=as.integer(fragments), fragmentfns=as.character(tmpfiles))
  } 
  res$newID=newID
  res
}