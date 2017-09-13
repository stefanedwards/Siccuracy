# Internals

#' Replaces extension of file names. Stolen directly from knitr.
#'
#' Might violate a license by including it here.
#' 
#' @param x Filenames, character vector.
#' @param ext File extensions to apply to \code{x}, without leading dot.
#' @references knitr
#' @author Yihui Xie
#' @keywords internal
#' @noRd
#' @rdname internal_sub_ext
sub_ext <- function (x, ext) {
  i = grep("\\.([[:alnum:]]+)$", x)
  x[i] = sans_ext(x[i])
  paste(x, ext, sep = ".")
}
#' @inheritParams sub_ext
#' @rdname internal_sub_ext
#' @noRd
sans_ext <- function (x, compression = FALSE) {
  if (compression) 
    x <- sub("[.](gz|bz2|xz)$", "", x)
  sub("([^.]+)\\.[[:alnum:]]+$", "\\1", x)
}


# @noRd
# @rdname parseformat
#'  Format descriptors
#'
#' Some functions allow the user to alter the output format by using arguments \code{int} and \code{format}.
#' The latter argument accepts a number of ways to define the output.
#' 
#' The output format is controlled by Fortran's format descriptiors; these usually take the form:
#' \describe{
#'  \item{Integers}{\code{Iw}}
#'  \item{Numeric}{\code{Fw.d}}
#' }
#' where \code{w} describes the width of the column (in characters) and \code{d} the number of decimals.
#' \strong{Note:} When using the \code{I} format code, you must specify \code{int=TRUE}. When using \code{int=TRUE}, real values are rounded to nearest integer.
#' 
#' The type of accepted values for \code{format} are:
#' \describe{
#'  \item{integer}{If a scalar integer, it specifies the width. Other integer vectors default to \code{I2}.}
#'  \item{numeric}{If a scalar numeric, the integer-part is used for width and the first decimal specifies the number of decimales. Other numeric vectors default to \code{F5.2}.}
#'  \item{character}{Used as-is, but prepended if necessary with either \code{I} or \code{F}.}
#' }
#' \strong{Note:} No checks are performed on the width and number of decimals.
#' 
#' @name parseformat
#' 
NULL

#' @noRd
parse.format <- function(format=NULL, int=TRUE) {
  if (is.null(format)) format <- ifelse(int, integer(1), numeric(1))
  if (is.integer(format)) {
    if (length(format) != 1) {
      format <- 'I2'
    } else if (format == 0) {
      format <- 'I2'
    } else {
      format <- as.character(format)
    }
  } else if (is.numeric(format)) {
    if (length(format) != 1) {
      format <- 'F5.2'
    } else if (format == 0.0) {
      format <- 'F5.2'
    } else {
      format <- formatC(format, digits=1, format='f')
    }
  }
  if (!grepl('^[[:alpha:]]+', format)) format <- paste0(ifelse(int, 'I', 'F'), format)

  format  
}

#' Pads a string to 255 characters, before sending off to Fortran
#' 
#' @noRd
as_fortran_character <- function(s) { sprintf('%-255s', s) }


# For testing.
make.true <- function(n,m,types=0:2) {
  true <- matrix(sample(types, size = n*m, replace=TRUE), ncol=m)  # fill true with random 0, 1, or 2.
  # add non-segregating site
  i <- sample.int(m, 1)
  true[2:n,i] <- true[1,i]
  rownames(true) <- as.character(1:nrow(true))
  true
}
make.imputed <- function(true) {
  m <- ncol(true)
  n <- nrow(true)
  imputed <- true
  imputed[sample.int(n*m, floor(n*m*0.5))] <- sample(0:2, size=floor(n*m*0.5), replace=TRUE) # change half the elements
  imputed[sample.int(n*m, floor(n*m*0.1))] <- NA  # some elements are missing.
  imputed
}
make.test <- function(n, m) {
  true <- make.true(n, m)
  imputed <- make.imputed(true)
  truefn <- tempfile('true', fileext = '.txt')
  imputedfn <- tempfile('imputed', fileext='.txt')
  write.snps(true, truefn)
  write.snps(imputed, imputedfn)
  list(true=true, imputed=imputed, truefn=truefn, imputedfn=imputedfn)
}
make.phase <- function(n,m) {
  x <- make.true(2*n, m, types=0:1) 
  rownames(x)[seq.int(2,nrow(x), by=2)] <- rownames(x)[seq.int(1,nrow(x), by=2)]
  x
}