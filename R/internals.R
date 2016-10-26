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


#' Parses formatting string.
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

# For testing.
make.true <- function(n,m) {
  true <- matrix(sample(0:2, size = n*m, replace=TRUE), ncol=m)  # fill true with random 0, 1, or 2.
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
