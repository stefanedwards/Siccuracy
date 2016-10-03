#' Auxiliary functions
#'
#' Short-hand functions for retriving file information such as number of columns or rows,
#' or picking the first column (relevant if this contains animal IDs).
#' @name auxfunc
#' @rdname auxfunc
NULL

#' Gets number of columns from first line of file.
#'
#' \code{get_ncols} works only with numeric values.
#'
#' @param fn Filename to read from. May be text, gz, or bz2.
#' @return \code{get_ncols}: Number of read fields, including ID columns.
#' @export
#' @rdname auxfunc
get_ncols <- function(fn) {
  on.exit(try(close(f), silent=TRUE))
  f <- gzfile(fn, 'r')
  s <- scan(f, what=character(), quiet=TRUE, nlines=1)
  close(f)
  length(s)
}


#' Gets number of rows in as file.
#'
#' \code{get_nlines} works only with numeric values.
#'
#' @param showWarning Logical, whether to display a warning if reading file is erroneous.
#' @param doError Logical, if error occured while reading file, will it stop
#' @export
#' @return \code{get_nlines}: Number of lines.
#' @rdname auxfunc
get_nlines <- function(fn, showWarning=TRUE, doError=FALSE) {
  stopifnot(file.exists(fn))
  res <- .Fortran('get_nlines', fn=as.character(fn), nlines=integer(1), stat=integer(1))
  if (res$nlines == 0 & res$stat != 0) {
    msg <- paste0('get_nlines did not read lines; stat error ', res$stat, '.')
    if (doError) { stop(msg)
    } else if (showWarning) warning(msg)
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
  
  res <- do.call(read.table, args)
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
#' @param fn Filename or connection of file to write to.
#' @param row.names If genotype matrix is "raw" and has first column with animals IDs, set this to False.
#' @param na Character to write in-place of missing values, defaults '9'.
#' @param ... Passed on to write.table.
#' @seealso \code{\link[utils]{write.table}}, \link[base]{connections}
#' @export
write.snps <- function(x, fn, row.names=TRUE, na='9', ...) {
  write.table(x, fn, col.names=FALSE, row.names=row.names, quote=FALSE, na=na, ...)
}

#' Read genotype matrix from file.
#'
#' \code{read.snps} is wrapper around \code{\link[base]{scan}} that also converts into native matrix format.
#' 
#' Assumes a file format with no header and first column are IDs. If no ID column, use \code{extractIDs = FALSE}.
#' Usually white-space delimted, but separator can be set with \code{sep} argument as per \code{\link[base]{scan}}.
#' For file format example see \link{Siccuracy}.
#'
#' @param file Name of file to read from or connection.
#' @param nlines Integer. If positive, the maximum number of line to read.
#' @param ncols Integer, default \code{NULL}. Number of columns in input file, including ID column. If \code{NULL} and \code{file} is string, automagically detected.
#' @param na If not \code{NA} (default), entries with the value are replaced with \code{NA}.
#' @param what Type of internal storage for matrix. Use \code{integer()} or \code{numeric()}, but expect issues if ID column contains alphabetical components.
#' @param extractIDs Logical, default \code{TRUE}, trim of first column and use as rownames?
#' @param quiet Logical, default \code{TRUE}. If \code{FALSE}, \code{scan} will print a line saying how may items have been read.
#' @param ... Passed on to \code{\link[base]{scan}}.
#' @return Native \link[base]{matrix}.
#' @export
#' @seealso  \code{\link{write.snps}}, \link[base]{typeof}
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
read.snps <- function(file, nlines=0, ncols=NULL, na=NA, what=integer(), extractIDs=TRUE, quiet=TRUE, ...) {
  if (is.null(ncols) & is.character(file)) {
    ncols <- get_ncols(file)
  }
  if (is.null(ncols)) stop('Cannot automagically detect number of columns to read as input file is a connection, not a character.')

  M <- matrix(scan(file, nlines=nlines, what=what, quiet=quiet, ...), ncol=ncols, byrow=TRUE)
  if (nrow(M) == 0) return(M)
  if (extractIDs) {
    rownames(M) <- M[,1]
    M <- M[,-1]
  }
  if (!is.na(na)) M[M==na] <- NA
  M
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
