library(testthat)
library(Siccuracy)

test_that('cbind_snp_files binds columns correctly', {
  #cat('cbind_SNPs binds columns correctly\n')
  nchr <- 5
  fns <- replicate(nchr, tempfile())
  # Different columns per SNP file.
  cols <- sample.int(10, size=nchr) + 10
  SNPs <- sapply(cols, Siccuracy:::make.true, n =8)
  null <- mapply(write.snps, file=fns, x=SNPs)
  
  total <- do.call(cbind, SNPs)
  
  fnout <- tempfile()
  stat <- cbind_snp_files(fns, fnout)
  res <- read.snps(fnout, what=integer())
  expect_equal(res, total)
  
  stat <- cbind_snp_files(fns, fnout)
  res <- read.snps(fnout, what=numeric())
  expect_equal(res, total)  
})

test_that('Skiplines skip lines', {
  #cat('Skiplines skip lines\n')
  nchr <- 5
  fns <- replicate(nchr, tempfile())
  # Different columns per SNP file.
  cols <- sample.int(10, size=nchr) + 10
  SNPs <- sapply(cols, Siccuracy:::make.true, n =16)
  null <- mapply(write.snps, file=fns, x=SNPs)
  
  total <- do.call(cbind, SNPs)
  
  fnout <- tempfile()
  stat <- cbind_snp_files(fns, fnout, skiplines=5)
  res <- read.snps(fnout, what=integer())
  expect_equal(res, total[-c(1:5),])  
})

test_that('Excluded IDs are not outputted',{
  #cat('Excluded IDs are not outputted','\n')
  nchr <- 5
  fns <- replicate(nchr, tempfile())
  # Different columns per SNP file.
  cols <- sample.int(10, size=nchr) + 10
  SNPs <- sapply(cols, Siccuracy:::make.true, n =16)
  null <- mapply(write.snps, file=fns, x=SNPs)
  
  total <- do.call(cbind, SNPs)
  
  exclude <- sample.int(nrow(total), 5)
  fnout <- tempfile()
  stat <- cbind_snp_files(fns, fnout, excludeids = exclude)
  res <- read.snps(fnout, what=integer())
  expect_equal(res, total[-exclude,])  
})

test_that('Mismatching ID are warned and thrown out',{
  #cat('Excluded IDs are not outputted','\n')
  nchr <- 5
  fns <- replicate(nchr, tempfile())
  # Different columns per SNP file.
  cols <- sample.int(10, size=nchr) + 10
  SNPs <- sapply(cols, Siccuracy:::make.true, n =16)
  rownames(SNPs[[1]])[c(3,5)] <- c(100, 105)
  null <- mapply(write.snps, file=fns, x=SNPs)
  
  total <- do.call(cbind, SNPs)
  
  fnout <- tempfile()
  expect_warning( stat <- cbind_snp_files(fns, fnout) )
  expect_equal(stat, nrow(total) - 2)
  res <- read.snps(fnout, what=integer())
  expect_equal(res, total[-c(3,5),])  
})

