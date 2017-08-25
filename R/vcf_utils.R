# VCFR format
# Developed using the vcfR package v. 1.4.0.
# August 2017

#' For VCF class objects, write genotype fields to \link{AlphaImpute-format}.
#' 
#' @rdname write.snps
#' @param element element to extract from vcf genotype, see \code{\link[vcfR]{extract.gt}}.
#' @param as.numeric logical, should the matrix be converted to numerics.
#' @param convertNA logical indicating whether to convert "." to \code{NA} 
#'                  (see \code{na} for which value to use when writing to text file).
#' @param newID Integer scalar (default \code{0}) for automatically assigning 
#'              new IDs. Use \code{NULL} to disable. See description for more. 
#' @import vcfR
#' @return For VCF class objects, returns data frame with ID mapping.
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
  
  # Assign new IDs
  if (!is.null(newIDS)) {
    if (is.integer(newID) && length(newID) == 1) {
      ids <- data.frame(vcr=colnames(gt), newID=seq_len(ncol(gt)) + newID)
      colnames(gt) <- ids$newID
    } else if (FALSE) {    
      # TO DO: Other ID mappings!
    }
  }
  
  write.snps.matrix(t(gt), file=file, na=na, ...)
  return(ids)
}