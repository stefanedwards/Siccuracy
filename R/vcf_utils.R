# VCFR format
# Developed using the vcfR package v. 1.4.0.
# August 2017

#' VCF objects in Siccuracy
#' 
#' To allow fluid interaction with the \emph{Variant Call Format}, 
#' some functions accepts the objects from the \link[=vcfR]{vcfR} package.
#' 
#' The following functions have been implemented 
#' with method dispatching for \link[=vcfR]{vcfR} objects:
#' 
#' \describe{
#'   \item{Imputation_accuracy}{\code{\link[=imputation_accuracy.vcfR]{imputation_accuracy}}}
#'   \item{Writing SNP data}{\code{\link[=write.snps.vcfR]{write.snps.vcfR}}}
#' }
#' 
#' @note 
#' No functions have been implemented for VCF \emph{files}.
#' The functions above all assume the files have been read in with
#' \code{\link[vcfR]{read.vcfR}}.
#' 
#' Direct conversion from VCF formatted files to AlphaImpute may be done by
#' first using PLINK to recode to 'A' format, and then use \code{\link{convert_plinkA}}.
#' For the first step, PLINK may be called as:
#' 
#' \code{plink --vcf <vcf-filename> --recode A --allow-extra-chr --out temp_file}
#' 
#' and then to convert to AlphaImpute format:
#' 
#' \code{convert_plinkA('temp_file.raw', 'alphaimpute_genotypes.txt')}
#' 
#' @name VCF_Siccuarcy
#' @aliases Siccuracy_VCF
#' @author Stefan McKinnon Edwards <sme@@iysik.com>, August 2017
NULL

#' Writes genotypes in VCF class objects to files in \link{AlphaImpute-format}.
#' 
#' @details 
#' The new integer IDs can be supplied. If not, they will be made for you.
#' \code{newID} may be an integer vector and will be used as is.
#' If data.frame with columns \code{sampID} and \code{newID}, they will be reordered to match input file.
#' 
#' @inheritParams write.snps
#' @param element element to extract from vcf genotype, see \code{\link[vcfR]{extract.gt}}.
#' @param as.numeric logical, should the matrix be converted to numerics.
#' @param convertNA logical indicating whether to convert "." to \code{NA} 
#'                  (see \code{na} for which value to use when writing to text file).
#' @param newID Integer scalar (default \code{0}) for automatically assigning 
#'              new IDs. Use \code{NULL} to disable. See description for more.
#' @param ... Arguments forwarded to \code{\link[vcfR]{extract.gt}} and \code{\link{write.snps}}. 
#' @import vcfR
#' @seealso \code{\link{write.snps}}
#' @return Returns data frame with ID mapping, truncated to those elements outputted.
#' @export
write.snps.vcfR <- function(x, 
                            file, 
                            na='9', 
                            newID=0,
                            element='GT', 
                            as.numeric=TRUE, 
                            convertNA=TRUE,
                            ... ) {
  gt <- vcfR::extract.gt(x, element=element, as.numeric=as.numeric, ...)
  # gt has by default rows as markers, and individuals in columns.
  
  # Assign new IDs; 
  # based on write.snps.matrix, which explains the terminolgy of 'firstcols'
  # and 'nlines'
  if (is.data.frame(newID)) stopifnot(all(c('sampID','newID') %in% names(newID)))
  
  firstcols <- data.frame('sampID'=colnames(gt), stringsAsFactors = FALSE)
  
  nlines <- ncol(gt)
  
  if (length(newID) == 1) {
    .newID <- firstcols
    .newID$newID <- 1:nrow(.newID) + as.integer(newID)
  } else if (is.atomic(newID)) {
    .newID <- firstcols
    .newID$newID <- 0
    .newID$newID[1:length(newID)] <- newID
  } else {
    newID$sampID <- as.character(newID$sampID)
    newID$newID <- as.integer(newID$newID)
    .newID <- merge(firstcols, newID, by=c('sampID'), sort=FALSE, all.x=TRUE, all.y=FALSE, stringsAsFactors=FALSE)
  }
  
  .newID <- .newID[match(colnames(gt), .newID$sampID) ,]
  if (anyNA(.newID$newID)) {
    r <- is.na(.newID$newID)
    .newID$newID[r] <- 1:sum(r) + max(.newID$newID[!r])
  }
  colnames(gt) <- as.character(.newID$newID)
  
  write.snps.matrix(t(gt), file=file, na=na, row.names=TRUE, ...)
  return(.newID)
}



#' @inheritParams imputation_accuracy.matrix
#' @param ... Arguments passed on to \code{\link[vcfR]{extract.gt}}.
#' @export
#' @rdname imputation_accuracy
imputation_accuracy.vcfR <- function(true, impute, standardized=TRUE, center=NULL, scale=NULL, p=NULL, excludeIDs=NULL, excludeSNPs=NULL, tol=0.1, ...) {
  true <- vcfR::extract.gt(true, element='GT', as.numeric=TRUE, ...)
  impute <- vcfR::extract.gt(impute, element='GT', as.numeric=TRUE, ...)
  imputation_accuracy(true, impute, 
                      standardized = standardized,
                      center=center,
                      scale=scale,
                      p=p,
                      excludeIDs=excludeIDs,
                      excludeSNPs=excludeSNPs,
                      tol=tol,
                      transpose=TRUE)
}