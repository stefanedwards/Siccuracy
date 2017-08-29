# Oxford format:
# A standard? Who knows.


#' Oxford format
#' 
#' Oxford format is a flat, single file, with a single loci per row.
#' The columns are chromosome, SNP ID, (chromosomal) position, 
#' reference allele and alternative allele. 
#' Following triplets of columns denote each sample's genotype dosage / probability.
#' 
#' @param file Filename to read from or write to.
#' @param ... Parameters forwarded to \code{\link[utils]{read.table}} or \code{\link[utils]{write.table}},
#' @export
#' @rdname oxford
#' @name Oxford format
#' @return \code{read.oxford} returns an \code{oxford} object, i.e. a list with two entries: 
#' \describe{
#'   \item{\code{map}:}{Data frame with map of loci, i.e. chromosome, snp ID, chromosomal position, and the two alleles.}
#'   \item{\code{probs}:}{3 dimensional array; rows (d1) are loci, columns (d2) are samples, in three layers comprising
#'                        gene dosage of homozygote for reference allele, heterozygote, and homozygote for alternative allele.}
#' }
#' @import abind abind
read.oxford <- function(file, ...) {
  args <- merge.list(list(...), list(as.is=TRUE, header=FALSE, row.names=NULL))
  args$file <- file
  oxford <- do.call(read.table, args)# read.table(file, as.is=TRUE, header=FALSE, row.names=NULL, ...)
  
  map <- oxford[,1:5]
  colnames(map) <- c('chr','snpID','pos','A1','A2')
  
  probs <- as.matrix(oxford[,-c(1:5)])
  storage.mode(probs) <- 'numeric'
  
  stopifnot(ncol(probs) %% 3 == 0)
  
  cols <- (1:ncol(probs)) %% 3
  probs <- abind(probs[,cols == 1], probs[,cols==2], probs[,cols==0], along=3)
  dimnames(probs) <- list(map$snpID, 1:ncol(probs), c('hom1','het','hom2'))
  
  res <- list(map=map, probs=probs)
  class(res) <- c('oxford','list')
  stopifnot(is.oxford(res))

  res
}

#' @param x Oxford object.
#' @rdname oxford
#' @export
is.oxford <- function(x) {
  !is.null(x$probs) && 
    !is.null(x$map) && 
    is.data.frame(x$map) &&
    is.array(x$probs) && 
    length(dim(x$probs)) == 3
}

#' @rdname oxford
#' @export
write.oxford <- function(x, file, ...) {
  args <- merge.list(list(...), list(as.is=TRUE, header=FALSE, row.names=NULL, file=file))
  args$x <- cbind(x$map, apply(x$probs, 1:2, paste, collapse=' '))
  do.call(write.table, args)
}

#' @inheritParams extract.gt
#' @param as.integer Rounds genotype dosages to whole integers.
#' @rdname extract.gt
#' @export
extract.gt.oxford <- function(x, as.integer=FALSE) {
  snps <- apply(x$probs, 1:2, reduce.oxford)
  if (as.integer) {
    snps <- round(snps)
    storage.mode(snps) <- 'integer'
  }
  
  snps
}
reduce.oxford <- function(x) {(x[2]*1 + x[3]*2) / sum(x)}

#' @inheritParams extract.gt.oxford
#' @rdname extract.snps
#' @export
extract.snps.oxford <- function(x, as.integer=FALSE) {
  t(extract.gt.oxford(x, as.integer=as.integer))
}

#' @rdname write.snps
#' @inheritParams write.snps.matrix
#' @param row.names If genotype matrix is "raw" and has first column with animals IDs, set this to \code{FALSE}.
#' @export
write.snps.oxford <- function(x, file, row.names=TRUE, na='9', ...) {
  write.snps.matrix(extract.snps.oxford(x, ...), file=file, row.names=row.names, na=na, ...)
}
  