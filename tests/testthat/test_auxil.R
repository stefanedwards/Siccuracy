library(testthat)
library(Siccuracy)

context('Testing auxillary functions')

test_that('`get_nlines` works and fails', {
  n <- sample.int(5000, 1)
  m <- 10
  r <- matrix(rnorm(n*m), ncol=m)
  fn <- tempfile(fileext='.txt')
  write.snps(r, fn)
  res <- get_nlines(fn)
  expect_equal(res, n)
  
  fn <- tempfile()
  expect_error(get_nlines(fn))
  
  writeLines('', fn)
  expect_equal(get_nlines(fn), 1) #because writeLines ends with new line.
  
  fn <- tempfile()
  writeLines('', fn, sep='')
  expect_equal(get_nlines(fn), 0)
  
  fn <- tempfile()
  s <- c('a','','b')
  writeLines(s, fn, sep='\n')
  expect_equal(get_nlines(fn), length(s))
  
  fn <- tempfile()
  s <- c('a','','','b')
  writeLines(s, fn, sep='\n')
  expect_equal(get_nlines(fn), length(s))
  
  fn <- tempfile()
  s <- c('a\n','\n','b')
  writeLines(s, fn, sep='')
  expect_equal(get_nlines(fn), length(s))
  
})

test_that('`get_ncolumns works', {
  n <- sample.int(5000, 1)
  m <- 10
  r <- matrix(rnorm(n*m), ncol=m)
  fn <- tempfile(fileext='.txt')
  write.snps(r, fn)
  res <- get_ncols(fn)
  expect_equal(res, m+1)

  # Also works on characters
  rownames(r) <- sample(combn(letters, 4, FUN=paste, collapse=''), n)
  fn <- tempfile(fileext='.txt')
  write.snps(r, fn, row.names=TRUE)
  expect_equal(get_ncols(fn), m+1, info='get_ncols works with characters')
  
  # Also with combination of tabs and spaces
  writeLines('a b\tc\t"d e f"\n12 3 ', con=fn) 
  expect_equal(get_ncols(fn), 4, info='get_ncols works with most')
})

test_that('`get_firstcolumn works', {
  df <- data.frame(a=1:10, b=letters[1:10], c=rnorm(10), d=c(1:5, letters[1:5]), stringsAsFactors = FALSE)
  fn <- tempfile(fileext='.txt')
  write.table(df, fn, append=FALSE, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
  
  expect_equal(get_firstcolumn(fn), df$a)
  expect_equal(get_firstcolumn(fn, 'numeric'), as.numeric(df$a))
  expect_equal(get_firstcolumn(fn, 'character'), as.character(df$a))
  
  expect_equal(get_firstcolumn(fn, class=c('NULL')), data.frame())
  expect_equal(get_firstcolumn(fn, class=c('NULL','character')), df$b)
  expect_equal(get_firstcolumn(fn, class=c('NULL','NULL','numeric')), df$c)
  expect_equal(get_firstcolumn(fn, class=c('NULL','NULL','NULL','character')), df$d)
  
  expect_error(get_firstcolumn(fn, class=c('NULL','NULL','NULL','numeric')))

  expect_equivalent(get_firstcolumn(fn, class=c('integer','character')), df[,1:2])
  
})

 test_that('Rcpp `get_nlines` and `get_ncols` (scan) can handle very, very long lines', {
   cols <- 20000
   rows <- 5
   m <- matrix(sample.int(9, cols*rows, replace=TRUE), ncol=cols)
   fn <- tempfile()
   write.snps(m, fn)
  
   expect_equal(get_nlines(fn), nrow(m))
   expect_equal(get_ncols(fn), ncol(m)+1)
})
