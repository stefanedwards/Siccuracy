# SHAPEIT haps/sample file formats:
# http://www.shapeit.fr/pages/m02_formats/hapssample.html
# SAMPLE file:
# three columns, famID, sampleID and #missing
# space delimited
# first row is ID_1 ID_2 missing
# second row is 0 0 0
# (?)
# following rows are samples
#
#
# HAPS file is:
# 1. Chromosome number [integer]
# 2. SNP ID [string]
# 3. SNP Position [integer]
# 4. First allele [string] (for first individual)
# 5. Second allele [string] (for first individual)
# 5+i:    Code of first allele for ith individual (ARRAYS START AT ZERO!).
# 5+i+1:  Code of second allele for ith individual (-//-)
# (5+i, 5+i+1) are the first and second allele for the ith individual (ARRAYS START AT ZERO!). Heh.


#' SHAPEIT haps/sample format
#' 
#' \code{read.haps} loads SHAPEIT haps/sample into a \code{haps} object.
#' 
#' \code{read.haps} accepts compressed files (.gz and .bz2, see \link[base]{connections}), 
#' but then both filename arguments must be given.
#' 
#' @param haps Filename of .haps file.
#' @param sample Filename of .sample file. 
#'               If \code{NULL} (default), uses filestem of \code{haps} appended with \code{.sample}.
#' @return \code{read.haps} returns a list with three entries: 
#' \describe{
#'   \item{\code{samples}:}{Data fram with contents of .sample file (three columns).}
#'   \item{\code{haps}:}{Matrix with unphased genotypes of each allele, with each loci per row and indivduals in tuples of columns.}
#'   \item{\code{map}:}{Data frame with map of loci, i.e. chromosome, snp ID, chromosomal position, and the two alleles.}
#' }
#' @rdname haps-format
#' @name SHAPEIT format
#' @seealso \code{\link{extract.gt}} for vcfR-style extraction of genotypes
#' @export
read.haps <- function(haps, sample=NULL) {
  if (is.null(sample)) {
    sample <- paste0(tools::file_path_sans_ext(haps), '.sample')
  }
  if (tools::file_ext(haps) == '') haps <- paste0(haps, '.haps')
  
  snps <- read.table(haps, header=FALSE, as.is=TRUE, row.names=NULL, check.names = FALSE)
  samples <- read.table(sample, header=TRUE, row.names=NULL, as.is=TRUE)[-1,]
  map <- snps[,1:5]
  colnames(map)[1:5] <- c('chromosome','snpID','pos','A1','A2')
  
  stopifnot(ncol(snps) == 5+2*nrow(samples))
  
  same <- with(samples, ID_1 == ID_2)
  IDs <- samples$ID_1
  IDs[!same] <- with(samples[!same,], paste(ID_1, ID_2, sep='_'))
  rownames(samples) <- IDs
  
  cols <- rep(NA, nrow(samples)*2)
  cols[seq.int(1, length(cols), by=2)] <- paste(IDs, 1, sep='_')
  cols[seq.int(2, length(cols), by=2)] <- paste(IDs, 2, sep='_')
  
  snps <- as.matrix(snps[,-c(1:5)])
  storage.mode(snps) <- 'integer'
  dimnames(snps) <- list(map[,2], cols)
  
  res <- list(samples=samples, haps=snps, map=map)
  res <- as.haps(res)
  stopifnot(is.haps(res))
  
  return(res)
}

#' Converts and tests if is object is of type \code{haps}
#' 
#' \code{haps} class requires a list object with three elements, 
#' \code{haps}, \code{sample}, and \code{map}, with corresponding dimensions.
#' 
#' \code{is.haps} doesn't actually test if \code{x} inherits from class \code{haps},
#' just that the contents matches.
#' 
#' @param x Object to test or convert
#' @rdname haps-format
is.haps <- function(x) {
  a <- all(c('samples','haps','map') %in% names(x))
  if (a) {
    b <- 2*nrow(x$samples) == ncol(x$haps) & 
         nrow(x$haps) == nrow(x$map) &
         ncol(x$samples) >= 3
    c <- all(c('ID_1','ID_2','missing') %in% names(x$samples))
  } else {
    b <- FALSE
  }
  return(a & b)
}

#' @export
#' @inheritParams is.haps
#' @rdname haps-format
as.haps <- function(x) {
  if (is.haps(x)) {
    class(x) <- c('haps', class(x))
    
    if (is.null(rownames(x$samples))) {
      same <- with(x$samples, ID_1 == ID_2)
      IDs <- x$samples$ID_1
      IDs[!same] <- with(x$samples[!same,], paste(ID_1, ID_2, sep='_'))
      rownames(x$samples) <- IDs
    }
  }
  x
}



#' \code{write.haps} writes haps object to haps/sample file formats
#' 
#' @rdname haps-format
#' @export
write.haps <- function(x, haps, sample) {
  stopifnot(is.haps(x))
  
  write.table(rbind(c('0','0','0'), x$sample), sample, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=' ')
  
  write.table(cbind(x$map, x$haps), haps, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=' ')
}




# Write files ------------------


#' @rdname write.snps
#' @inheritParams write.snps.matrix
# @param phased Logical, when \code{FALSE} (default), collapses all loci into single-value
#               where haplotype information is lost.
#' @param newID ... New ID, see \code{\link{convert_plinkA}}, but with some caveats here. I don't know.
#' @export
write.snps.haps <- function(x, file, row.names=TRUE, na='9', newID=0, ...) {
  # if (phased) {
  #   gt <- extract.phased(x, as.numeric=TRUE)
  # } else {
    gt <- extract.gt(x, as.numeric=TRUE)
#  }
  
  # Handle IDs - code based of that for convert_plinkA
  if (is.data.frame(newID)) {
    if (all(!c('ID','newID') %in% names(newID))) {
      stopifnot(all(c('ID_1','ID_2','newID') %in% names(newID)))
      
      same <- with(newID, ID_1 == ID_2)
      newID$ID <- newID$ID_1
      newID$ID[!same] <- with(newID, paste(ID_1, ID_2, sep='_'))
      rownames(newID) <- newID$ID
    }
  }
  
  
  firstcols <- x$samples
  if (length(newID) == 1) {
    .newID <- firstcols
    .newID$newID <- 1:nrow(.newID) + newID
  } else if (is.atomic(newID)) {
    nlines <- min(nrow(firstcols), length(newID))
    .newID <- cbind(firstcols[1:nlines,], newID[1:nlines])
  } else {
    newID <- do.call(data.frame, append(lapply(newID, as.character), list(stringsAsFactors=FALSE)))
    .newID <- merge(firstcols, newID, by=0, sort=FALSE, all.x=TRUE, all.y=FALSE, stringsAsFactors=FALSE)
  }
  
  .newID <- .newID[match(colnames(gt), rownames(.newID)) ,]
  if (anyNA(.newID$newID)) {
    r <- is.na(.newID$newID)
    .newID$newID[r] <- 1:sum(r) + max(.newID$newID[!r])
  }
  colnames(gt) <- as.character(.newID$newID)

  
  write.snps.matrix(t(gt), file=file, row.names=TRUE, na=na, ...)
  .newID
}


# Imputation accuracy -------------------

#' @inheritParams imputation_accuracy.matrix
#' @export
#' @rdname imputation_accuracy
imputation_accuracy.haps <- function(true, impute, standardized=TRUE, center=NULL, scale=NULL, p=NULL, tol=0.1, ..., excludeIDs=NULL, excludeSNPs=NULL) {
  true <- extract.snps.haps(true, ...)
  impute <- extract.snps.haps(impute, ...)
  imputation_accuracy(true, impute, 
                      standardized = standardized,
                      center=center,
                      scale=scale,
                      p=p,
                      excludeIDs=excludeIDs,
                      excludeSNPs=excludeSNPs,
                      tol=tol,
                      transpose=FALSE)
}
