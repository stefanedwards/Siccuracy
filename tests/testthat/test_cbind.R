library(testthat)
library(Siccuracy)

test_that('cbind_SNPs binds columns correctly', {
  nchr <- 5
  fns <- replicate(nchr, tempfile())
  # Different columns per SNP file.
  cols <- sample.int(10, size=nchr) + 10
  SNPs <- sapply(cols, Siccuracy:::make.true, n =8)
  null <- mapply(write.snps, fn=fns, x=SNPs)
  
  total <- do.call(cbind, SNPs)
  
  fnout <- tempfile()
  stat <- cbind_SNPs(fns, fnout, int=TRUE)
  res <- read.snps(fnout, what=integer())
  expect_equal(res, total)
  
  stat <- cbind_SNPs(fns, fnout, int=FALSE)
  res <- read.snps(fnout, what=numeric())
  expect_equal(res, total)  
})

test_that('Skiplines skip lines', {
  nchr <- 5
  fns <- replicate(nchr, tempfile())
  # Different columns per SNP file.
  cols <- sample.int(10, size=nchr) + 10
  SNPs <- sapply(cols, Siccuracy:::make.true, n =16)
  null <- mapply(write.snps, fn=fns, x=SNPs)
  
  total <- do.call(cbind, SNPs)
  
  fnout <- tempfile()
  stat <- cbind_SNPs(fns, fnout, int=TRUE, skiplines=5)
  res <- read.snps(fnout, what=integer())
  expect_equal(res, total[-c(1:5),])  
})

test_that('Excluded IDs are not outputted',{
  nchr <- 5
  fns <- replicate(nchr, tempfile())
  # Different columns per SNP file.
  cols <- sample.int(10, size=nchr) + 10
  SNPs <- sapply(cols, Siccuracy:::make.true, n =16)
  null <- mapply(write.snps, fn=fns, x=SNPs)
  
  total <- do.call(cbind, SNPs)
  
  exclude <- sample.int(nrow(total), 5)
  fnout <- tempfile()
  stat <- cbind_SNPs(fns, fnout, int=TRUE, excludeids = exclude)
  res <- read.snps(fnout, what=integer())
  expect_equal(res, total[-exclude,])  
})

test_that('Numerical formats are respected',{
  nchr <- 5
  fns <- replicate(nchr, tempfile())
  # Different columns per SNP file.
  cols <- sample.int(10, size=nchr) + 10
  SNPs <- sapply(cols, Siccuracy:::make.true, n =16)
  null <- mapply(write.snps, fn=fns, x=SNPs)
  
  total <- do.call(cbind, SNPs)
  
  fnout <- tempfile()
  n <- cbind_SNPs(fns, fnout, int=FALSE, format='8.4')
  s <- scan(fnout, what='character', quiet=TRUE, nlines=1)[-1]
  expect(all(nchar(s) == 6), 'Wrong number of characters in written output.')  # spaces are stripped; expected character length is 4+2 (decimals + period + leading digit).
  
  n <- cbind_SNPs(fns, fnout, int=FALSE, format='10.8') 
  s <- scan(fnout, what='character', quiet=TRUE)
  expect(length(s) == nrow(total), 'Space as separators are not missing.')
  
})