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