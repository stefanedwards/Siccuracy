# Extract.gt --------------

#' Extract genotypes from SNP objects.
#' 
#' gt genotypes are stored with loci per row and samples per column.
#' For the transposed, see \code{\link{extract.snps}}.
#'
#' @export
#' @param x Object to extract genotypes from.
#' @param as.numeric Logical, when \code{TRUE} (default), return numeric matrix
#'        of SNP counts. When \code{FALSE}, character vector with both alleles 
#'        separated by '/'.
#' @param ... Arguments passed on to \code{vcfR}'s \code{\link[vcfR]{extract.gt}}.
#' @rdname extract.gt
#' @name extract.gt
#' @seealso \code{\link{extract.snps}}
extract.gt <- function(x, as.numeric, ...) {
  UseMethod('extract.gt', x)
}

#' @rdname extract.gt
#' @export
extract.gt.vcfR <- function(x, as.numeric, ...) {
  vcfR::extract.gt(x, as.numeric=as.numeric, ...)
}

#' @return \code{extract.gt} returns a matrix with loci per row and samples per column.
#' @rdname extract.gt
#' @export
extract.gt.haps <- function(x, as.numeric=TRUE, ...) {
  m <- ncol(x$haps)
  if (as.numeric) {
    snps <- x$haps[,seq.int(1,m,by=2)] + x$haps[,seq.int(2,m,by=2)]
    dimnames(snps) <- list(rownames(x$haps), rownames(x$samples))
  } else {
    alleles <- as.matrix(x$map[,4:5])
    snps <- paste(apply(x$haps[,seq.int(1,m,by=2)], 2, map2allele, alleles=alleles),
                  apply(x$haps[,seq.int(2,m,by=2)], 2, map2allele, alleles=alleles),
                  sep='/')
    snps <- matrix(snps, ncol=m/2, dimnames=list(rownames(x$haps), rownames(x$samples)))
  }
  snps
}

map2allele <- function(column, alleles) { alleles[cbind(1:nrow(alleles), column+1)] }

#' @inheritParams extract.gt
#' @param as.integer Rounds genotype dosages to whole integers.
#' @rdname extract.gt
#' @export
extract.gt.oxford <- function(x, as.numeric=TRUE, ..., as.integer=FALSE) {
  if (!as.numeric) 
    warning('Oxford format does not support non-numeric coding. `as.numeric` ignored.')
  
  snps <- apply(x$probs, 1:2, reduce.oxford)
  if (as.integer) {
    snps <- round(snps)
    storage.mode(snps) <- 'integer'
  }
  
  snps
}
reduce.oxford <- function(x) {(x[2]*1 + x[3]*2) / sum(x)}



#' #' Extracts 
#' #' 
#' #' @inheritParams extract.snps
#' #' @return \code{extract.phased} returns matrix with individuals by row, with two rows
#' #'   per individual for first and second allele.
#' extract.phased <- function(x) {
#'   snps <- t(x$haps)
#'   rownames(snps) <- rep(rownames(x$samples), each=2)
#'   snps
#' }

# Extract SNPs -----------------------
# SNPs are set in columns and individuals by row.

#' Extract genotypes from SNP objects, in AlphaImpute format
#'
#' The SNP format lists samples as rows and loci as columns. 
#' it is simply \code{gt} transposed.
#' The SNP format is \emph{always} numeric.
#' 
#' @param x Object to extract genotypes from.
#' @rdname extract.snps
#' @name extract.snps
#' @return \code{extract.snps} returns matrix with indivduals by row, with numeric genotypes.
#' @export
extract.snps <- function(x, as.integer) {
  UseMethod('extract.snps', x)
}

#' @rdname extract.snps
#' @export
extract.snps.haps <- function(x, as.integer=FALSE) {
  snps <- t(extract.gt.haps(x, is.numeric=TRUE))
  if (as.integer) {
    snps <- round(snps)
    storage.mode(snps) <- 'integer'
  }
  snps
}

#' @inheritParams extract.snps
#' @inheritParams extract.gt.vcfR
#' @rdname extract.snps
#' @export
extract.snps.vcfR <- function(x, as.integer=FALSE, ...) {
  snps <- t(extract.gt.vcfR(x, is.numeric=TRUE, ...))
  if (as.integer) {
    snps <- round(snps)
    storage.mode(snps) <- 'integer'
  }
  snps
}


#' @inheritParams extract.snps.oxford
#' @rdname extract.snps
#' @export
extract.snps.oxford <- function(x, as.integer=FALSE) {
  t(extract.gt.oxford(x, as.integer=as.integer))
}