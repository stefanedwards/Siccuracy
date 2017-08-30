library(testthat)
library(Siccuracy)

context('SHAPEIT haps/sample file format')

test_that('Reading haps/sample format', {
  s <- system.file('extdata/testdata/simulated.haps', package='Siccuracy', mustWork=TRUE)
  
  h <- read.haps(s)
  
  # Basically contents of is.haps
  expect_true(all(c('ID_1','ID_2','missing') %in% names(h$samples)))
  expect_equal(2*nrow(h$samples), ncol(h$haps))
  expect_equal(nrow(h$haps), nrow(h$map))
  expect_gte(ncol(h$samples), 3)

  
})

test_that('Writing haps/sample format in haps/sample format', {
  s <- system.file('extdata/testdata/simulated.haps', package='Siccuracy', mustWork=TRUE)
  
  h <- read.haps(s)
  
  t1 <- tempfile()
  t2 <- tempfile()
  
  write.haps(h, t1, t2)
  
  s1 <- readLines(s)
  t1 <- readLines(t1)
  expect_equal(s1, t1)
  
  s <- system.file('extdata/testdata/simulated.sample', package='Siccuracy', mustWork=TRUE)
  s2 <- readLines(s)
  t2 <- readLines(t2)
  expect_equal(s2, t2)
})

test_that('Writing haps/sample objects to snp format', {
  s <- system.file('extdata/testdata/simulated.haps', package='Siccuracy', mustWork=TRUE)
  
  h <- read.haps(s)
  
  t1 <- tempfile()
  newID <- write.snps(h, t1)
  x <- read.snps(t1)
  
  expect_equal(nrow(x), nrow(h$samples))
  expect_equal(ncol(x), nrow(h$haps))
  
})

