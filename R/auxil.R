#' Auxiliary functions
#'
#' Short-hand functions for retriving file information such as number of columns or rows,
#' or picking the first column (relevant if this contains IDs).
#' @name auxfunc
#' @rdname auxfunc
NULL

#' Gets number of columns from first line of file.
#'
#' \code{get_ncols} returns number of columns.
#'
#' @param file Filename or connection to read from.
#' @return \code{get_ncols}: Number of read fields, including ID columns.
#' @export
#' @rdname auxfunc
get_ncols <- function(file) {
  #on.exit(try(close(f), silent=TRUE))
  #f <- gzfile(fn, 'r')
  s <- scan(file, what=character(), quiet=TRUE, nlines=1)
  #close(f)
  length(s)
}


#' Gets number of rows in as file.
#'
#' \code{get_nlines} returns number of lines.
#'
#' @param fn Filename.
#' @export
#' @return \code{get_nlines}: Number of lines, or \code{NA} on error.
#' @rdname auxfunc
#' @section Error messages:
#' \code{get_nlines} is implemented in Fortran and may return \code{NA} with an error message, if an error is encountered.
#' 
#' IOSTAT errors:
#' \describe{
#'  \item{5010}{The first column most likely contains non-integer or non-numeric elements.}
#' }
get_nlines <- function(fn) {
  stopifnot(file.exists(fn))
  res <- .Fortran('get_nlines', fn=as.character(fn), nlines=integer(1), stat=integer(1))
  if (res$nlines == 0 & res$stat != 0) {
    warning(paste0('get_nlines did not read lines; IOSTAT error ', res$stat, '.'))
    return(structure(NA, code=res$stat))
  }
  res$nlines
}

#' Reads first column(s)
#'
#' \code{get_firstcolumn} is a wrapper for \code{read.table} that allows the user
#' to easily get the first column(s) by padding \code{colClasses} with \code{"NULL"} character elements to
#' skip columns.
#'
#' @export
#' @param class Character vector of classes to use for column(s), see \code{colClasses} in \code{\link[utils]{read.table}}.
#' @param ... Parameters sent to \code{read.table}.
#' @return \code{get_firstcolumn}: Vector with elements of first column, or \code{data.frame} if \code{class} contains multiple not-\code{"NULL"} elements.
#' @rdname auxfunc
get_firstcolumn <- function(fn, class='integer', ...) {
  args <- list(...)
  args$showWarning <- NULL
  args$doError <- NULL
  cols <- get_ncols(fn)
  
  classes <- rep('NULL', times=get_ncols(fn))
  classes[1:length(class)] <- class

  args$file <- fn
  args$colClasses=classes
  
  if (!is.null(args$col.names)) {
    if (length(args$col.names) != cols)
      args$col.names <- c(args$col.names, paste0('X', 1:(cols-length(args$col.names))))
  }
  
  res <- do.call(utils::read.table, args)
  if (ncol(res) == 1) res <- res[,1]
  return(res)
}

NULL

#' Write genotype matrices to file.
#'
#' \code{write.snps} is short hand for \code{write.table} with some default options.
#' For file format see \link{Siccuracy}.
#' 
#' @param x The matrix to write.
#' @param file Filename or connection to write to.
#' @param row.names If genotype matrix is "raw" and has first column with animals IDs, set this to \code{FALSE}.
#' @param na The string to use for \code{NA} in the data.
#' @param ... Passed on to write.table.
#' @seealso \code{\link[utils]{write.table}}, \link[base]{connections}
#' @export
write.snps <- function(x, file, row.names=TRUE, na='9', ...) {
  utils::write.table(x, file, col.names=FALSE, row.names=row.names, quote=FALSE, na=na, ...)
}

#' Read genotype matrix from file.
#'
#' \code{read.snps} is wrapper around \code{\link[base]{scan}} that also converts into native matrix format.
#' 
#' Assumes a file format with no header and first column are IDs. If no ID column, use \code{extractIDs = FALSE}.
#' Usually white-space delimted, but separator can be set with \code{sep} argument as per \code{\link[base]{scan}}.
#' For file format example see \link{Siccuracy}.
#' 
#' If ID columns contains alphabetical elements, use \code{what=character()}. This will however return a matrix. 
#' Use e.g. \code{storage.mode(m) <- 'integer'} to convert to integer-based, keeping all other attributes.
#' 
#'
#' @param file Filename or connection.
#' @param ncols Integer, number of columns in \code{file}, including ID column. When \code{NULL} (default), automagically detected with \code{\link{get_ncols}}.
#' @param na Missing values; entries with this value are replaced with \code{NA}..
#' @param what The \link[base]{typeof} of data to be read, e.g. \code{integer()},  \code{numeric()}, or \code{character()}.
#' @param extractIDs Logical, default \code{TRUE}, trim of first column and use as rownames.
#' @param quiet logical: if \code{FALSE}, \code{scan()} will print a line, saying how many items have been read.
#' @param ... Passed on to \code{\link[base]{scan}}.
#' @return Native \link[base]{matrix}.
#' @export
#' @seealso  \code{\link{write.snps}}, \link[base]{typeof}, \code{\link{get_ncols}}, \code{\link[base]{scan}}.
#' @examples
#'
#' # Make test data
#' tmpfile <- tempfile()
#' write.snps(Siccuracy:::make.true(5, 10), tmpfile)
#' M <- read.snps(tmpfile)
#'
#' # Process the genotypes in chunks:
#' # Setting mode to 'r' is very important to avoid resetting the pointer at file head!
#' ncols <- get_ncols(tmpfile)
#' f <- file(tmpfile, 'r') 
#' while (TRUE) {
#'   M <- read.snps(f, nlines=500, ncols=ncols)
#'   if (nrow(M) == 0) break
#' }
#' close(f)
read.snps <- function(file, ncols=NULL, na=NA, what=integer(), extractIDs=TRUE, quiet=TRUE, ...) {
  if (is.null(ncols) & is.character(file)) {
    ncols <- get_ncols(file)
  }
  if (is.null(ncols)) stop('Cannot automagically detect number of columns to read as input file is a connection, not a character.')

  M <- matrix(scan(file, what=what, quiet=quiet, ...), ncol=ncols, byrow=TRUE)
  if (nrow(M) == 0) return(M)
  if (extractIDs) {
    rownames(M) <- M[,1]
    M <- M[,-1,drop=FALSE]
  }
  if (!is.na(na)) M[M==na] <- NA
  M
}
