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
  stat <- cbind_SNPs(fns, fnout)
  res <- read.snps(fnout, what=numeric())
  expect_equal(res, total)
  
  stat <- cbind_SNPs(fns, fnout, format=integer())
  res <- read.snps(fnout, what=int())
})
