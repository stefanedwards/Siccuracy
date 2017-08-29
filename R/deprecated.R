#' Deprecated functions
#' 
#' @name Siccuracy-deprecated
#' @rdname deprecated
NULL

# cbind ------
#' Deprecated name, use \code{cbind_snp_files}
#' 
#' \code{cbind_SNPs} has been renamed to \code{\link{cbind_snp_files}}.
#' 
#' @rdname deprecated
#' @export
#' @inheritParams  cbind_snp_files
cbind_SNPs <- function(fns, fnout, nlines, ncols, skiplines, excludeids, int) {
  .Deprecated('cbind_snp_files', package='Siccuracy')
  cbind_snp_files(fns=fns, fnout=fnout, nlines=nlines, ncols=ncols, skiplines=skiplines, excludeids=excludeids, int=int)
}

#' Deprecated name, use \code{cbind_snp_files}.
#' 
#' \code{rowconcatenate} has been renamed to \code{\link{cbind_snp_files}}.
#' 
#' @export
#' @rdname deprecated
#' @inheritParams cbind_snp_files
rowconcatenate <- function(fns, fnout, nlines=NULL, ncols=NULL, skiplines=0, excludeids=integer(0)) {
  .Deprecated('cbind_snp_files', package='Siccuracy')
  cbind_snp_files(fns=fns, fnout=fnout, nlines=nlines, ncols=ncols, skiplines=skiplines, excludeids=excludeids, int=int)
}

# Imputation accuracy ------

#' Deprecated names, use \code{imputation_accuracy}.
#' 
#' \code{imputation_accuracy1} and \code{imputation_accuracy3} has been replaced by \code{\link{imputation_accuracy}}.
#' The difference between the two former functions is now covered by the \code{adaptive}-argument of the latter.
#' 
#' @param truefn,imputefn Deprected, was \emph{filename} to files with true and imputed genotype matrix.
#' @param nSNPs Deprecated, use \code{ncol}.
#' @param nAnimals Deprecated, use \code{nlines}.
#' @param NAval Deprecated, use \code{na}.
#' 
#' @export
#' @rdname deprecated
#' @inheritParams imputation_accuracy
imputation_accuracy3 <- function(truefn, imputefn, nSNPs=NULL, nAnimals=NULL, NAval=9, standardized=TRUE) {
  .Deprecated('imputation_accuracy', package='Siccuracy')
  with(imputation_accuracy(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, adaptive=TRUE),
       list(means=snps$means, sds=snps$sds, rowcors=animals$cors, matcor=matcor))
}

#' @export
#' @rdname deprecated
imputation_accuracy1 <- function(truefn, imputefn, nSNPs=NULL, nAnimals=NULL, NAval=9, standardized=TRUE) {
  .Deprecated('imputation_accuracy', package='Siccuracy')
  with(imputation_accuracy(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, adaptive=FALSE),
       list(means=snps$means, sds=snps$sds, rowcors=animals$cors, matcor=matcor))
}


# MergeChips -----

#' Deprecated name, use \code{rbind_snp_files}.
#' 
#' \code{mergeChips} has been renamed to \code{\link{rbind_snp_files}}.
#' 
#' @export
#' @param outpos Integer vector of collective SNP positions. Default to sorted, union of \code{hdpos} and \code{ldpos}. Make it anything else and you get?
#' @rdname deprecated
#' @inheritParams rbind_SNPs
#' @param missing Missing value.
mergeChips <- function(hdid,ldid, hdpos, ldpos, hdfn, ldfn, fnout, outpos=NULL, missing=9) {
  .Deprecated('rbind_snp_files', package='Siccuracy')
}

# Phases -------------
#' Deprecated name, use \code{convert_phases}.
#' 
#' \code{phasotogeno} and \code{phasotogeno_int} have been replaced by 
#' \code{\link{convert_phases}}. 
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



# rbind -----
#' Deprecated name, use \code{rbind_snp_files}
#'
#' \code{rbind_SNPs} has been renamed to \code{\link{rbind_snp_files}}.
#' 
#' @rdname deprecated
#' @export
#' @inheritParams rbind_snp_files
rbind_SNPs <- function(hdid,ldid, hdpos, ldpos, hdfn, ldfn, fnout, outcol=NULL, na=9, format=NULL, int=TRUE) {
  .Deprecated('rbind_snp_files', package='Siccuracy')
  rbind_snp_files(hdid=hdid, ldid=ldid, hdpos=hdpos, ldpos=ldpos, hdfn=hdfn, ldfn=ldfn, fnout=fnout, outcol=outcol, na=na, format=format, int=int)
}