# Extract.gt --------------

#' Extract genotypes from SNP objects.
#' 
#' gt genotypes are stored with loci per row and samples per column.
#' For the transposed, see \code{\link{extract.snps}}.
#'
#' @export
#' @rdname extract.gt
#' @name extract.gt
#' @seealso \code{\link{extract.snps}}
extract.gt <- function(x, ...) {
  UseMethod('extract.gt', x)
}

#' @rdname extract.gt
#' @export
extract.gt.vcfR <- function(x, ...) {
  vcfR::extract.gt(x, ...)
}

#' @param x Object to extract genotypes or phases from.
#' @param is.numeric Logical, when \code{TRUE} (default), return numeric matrix
#'        of SNP counts. When \code{FALSE}, character vector with both alleles 
#'        separated by '/'.
#' @return \code{extract.gt} returns a matrix with loci per row and samples per column.
#' @rdname extract.gt
#' @export
extract.gt.haps <- function(x, is.numeric=TRUE, ...) {
  m <- ncol(x$haps)
  if (is.numeric) {
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

#' Extracts 
#' 
#'
#' @return \code{extract.phased} returns matrix with individuals by row, with two rows
#'   per individual for first and second allele.
extract.phased <- function(x) {
  snps <- t(x$haps)
  rownames(snps) <- rep(rownames(x$samples), each=2)
  snps
}


#' Extract genotypes from SNP objects, in AlphaImpute format
#' @rdname extract.snps
#' @name extract.snps
#' @export
extract.snps <- function(x, ...) {
  UseMethod('extract.snps', x)
}


#' @rdname extract.snps
#' @export
#' @return \code{extract.snps} returns matrix with indivduals by row, with numeric genotypes.
extract.snps.haps <- function(x) {
  t(extract.gt.haps(x, is.numeric=TRUE))
}

#' @rdname extract.snps
#' @export
extract.snps.vcfR <- function(x) {
  t(extract.gt.vcfR(x, is.numeric=TRUE))
}