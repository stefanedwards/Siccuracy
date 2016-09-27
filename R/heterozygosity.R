

#' Calculates hetereozygosity stats for a single population
#'
#' Calculate gene frequencies and observed and expected heterozygosites for locus.
#' These routines only work on a genotype file (e.g. as written by \code{\link{write.snps}}) with first column as ID column,
#' and subsequent columns are genotypes coded as \code{0}, \code{1}, or \code{2}.
#'
#' The gene frequencies, \eqn{p} and \eqn{q}, refers to the gene frequencies of genotypes \code{0} and \code{2}, respectively.
#' They are calculated as such, for each column:
#' \deqn{
#'   p = \frac{\# 0 \cdot 2 + \# 1 }{2 \cdot \# rows }
#' }{p = (2 * # homozygotes + # heterozygotes) / (2 * # rows)}
#' \deqn{
#'   q = 1 - p
#' }
#'
#' Observed (\eqn{H_{obs}}{Hobs}) and expected (\eqn{H_{exp}}{Hexp}) heterozygosity are calculated as
#' \deqn{
#'   H_{obs} = \frac{\# 1}{\# rows}
#'   }{Hobs = #1/n}
#' \deqn{
#'   H_{exp} = 2 \cdot p \cdot q
#'   }{Hexp = 2*p*q}
#'
#' In the above, \code{NA} elements are ignored and thus do not count toward either genotype nor number of rows.
#'
#' Finally, an inbreeding coeffiecent would be calculated as:
#' \deqn{
#'   F = \frac{H_{exp} - H_{obs}}{H_{exp}}
#' }{F = (Hexp - Hobs)/Hexp}
#'
#' @export
#' @param fn Filename of genotype matrix (0,1,2) with first column denoting ID.
#' @param ncols Number of SNP columns in \code{fn}; if \code{NULL} (default), number is retrieved with \code{\link{get_ncols}}.
#' @param NAval Integer value for unknown values, which are ignored.
#' @return Data frame with columns
#' \describe{
#'   \item{\code{p}}{Vector of allele frequencies of alleles coded as '0'.}
#'   \item{\code{q}}{Vector of allele frequencies of alleles coded as '1'.}
#'   \item{\code{Hobs}}{Vector of observed heterozygosity (\eqn{H_{obs}}{Hobs}).}
#'   \item{\code{Hexp}}{Vector of expected heterozygosity (\eqn{H_{exp}}{exp}).}
#'   \item{\code{n}}{Vector of observed genotypes for each locus (ignoring NA-values). Alleles are twice this number.}
#' }
#' @references Equations based on \url{http://www.uwyo.edu/dbmcd/molmark/practica/fst.html}.
heterozygosity <- function(fn, ncols=NULL, NAval=9) {
  stopifnot(file.exists(fn))
  
  if (is.null(ncols)) {
    ncols <- get_ncols(fn) - 1
  }
  ncols <- as.integer(ncols)
  # subroutine heterozygosity(fn, ncols, NAval, p, q, Hobs, Hexp, n)
  res <- .Fortran('heterozygosity', fn=as.character(fn), ncols=as.integer(ncols), NAval=as.integer(NAval),
                  p=numeric(ncols), q=numeric(ncols), Hobs=numeric(ncols), Hexp=numeric(ncols), n=integer(ncols))
  
  as.data.frame(res[c('p','q','Hobs','Hexp','n')])
}