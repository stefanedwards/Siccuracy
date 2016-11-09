library(testthat)
library(Siccuracy)

test_that('Masking SNPs works...', {
  SNPs <- Siccuracy:::make.true(10, 10)
  snpfile <- tempfile()
  write.snps(SNPs, snpfile)
  
  maskIDs <- 5:9
  maskSNPs <- sample.int(ncol(SNPs), ncol(SNPs)*.3)
  
  fn <- tempfile()
  res <- mask_SNPs(snpfile, fn, maskIDs = maskIDs, maskSNPs = maskSNPs, na=9)
  m <- read.snps(fn)

  SNPs[maskIDs, maskSNPs] <- 9
  expect_equal(SNPs, m)
})

test_that('Dropping individuals works', {
  SNPs <- Siccuracy:::make.true(15, 10)
  snpfile <- tempfile()
  write.snps(SNPs, snpfile)
  
  maskIDs <- sample.int(nrow(SNPs), nrow(SNPs) * .5)
  dropIDs <- sample(maskIDs, length(maskIDs) * 0.5)
  
  fn <- tempfile()
  res <- mask_SNPs(snpfile, fn, maskIDs = maskIDs, maskSNPs=1, dropIDs=dropIDs, na=9)
  m <- read.snps(fn)
  
  expect_equal(nrow(m), nrow(SNPs) - length(dropIDs))
  expect_true(all(!dropIDs %in% rownames(m)))
  expect_true(all(m[as.character(setdiff(maskIDs, dropIDs)),1] == 9))  
  
})

test_that('Dropping SNPs works', {
  SNPs <- Siccuracy:::make.true(10, 20)
  snpfile <- tempfile()
  write.snps(SNPs, snpfile)
  
  maskIDs <- sample.int(nrow(SNPs), nrow(SNPs) * .3)
  maskSNPs <- sample.int(ncol(SNPs), ncol(SNPs) * .5)
  dropSNPs <- sample(maskSNPs, length(maskSNPs) * 0.5)  
  
  fn <- tempfile()
  res <- mask_SNPs(snpfile, fn, maskIDs=maskIDs, maskSNPs=maskSNPs, dropSNPs=dropSNPs, na=9)
  m <- read.snps(fn)
  
  SNPs[maskIDs,maskSNPs] <- 9
  SNPs <- SNPs[,-dropSNPs]
  
  expect_equal(ncol(m), ncol(SNPs) - length(dropSNPs))
  expect_equal(SNPs, m)
})