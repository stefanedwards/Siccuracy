library(testthat)
library(Siccuracy)

test_that('Masking SNPs works...', {
  SNPs <- Siccuracy:::make.true(10, 10)
  snpfile <- tempfile()
  write.snps(SNPs, snpfile)
  
  maskIDs <- 5:9
  maskSNPs <- sample.int(ncol(SNPs), ncol(SNPs)*.3)
  
  fn <- tempfile()
  res <- mask_snp_file(snpfile, fn, masking=list(maskIDs, maskSNPs))
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
  res <- mask_snp_file(snpfile, fn, masking=list(maskIDs, 1), dropIDs=dropIDs, na=9)
  m <- read.snps(fn)
  
  expect_equal(nrow(m), nrow(SNPs) - length(dropIDs))
  expect_true(all(!dropIDs %in% rownames(m)))
  expect_true(all(m[as.character(setdiff(maskIDs, dropIDs)),1] == 9))  
  
})

# Dropping SNPs works ----
test_that('Dropping SNPs works', {
  SNPs <- Siccuracy:::make.true(8, 9)
  snpfile <- tempfile()
  write.snps(SNPs, snpfile)
  
  maskIDs <- sample.int(nrow(SNPs), nrow(SNPs) * .3)
  maskSNPs <- sample.int(ncol(SNPs), ncol(SNPs) * .5)
  dropSNPs <- sample(maskSNPs, length(maskSNPs) * 0.5)  
  snps <- setdiff(1:ncol(SNPs), dropSNPs)
  
  
  fn <- tempfile()
  res <- mask_snp_file(snpfile, fn, snps=snps, masking=list(maskIDs, maskSNPs), na=9)
  m <- read.snps(fn)
  
  SNPs[maskIDs,maskSNPs] <- 9
  SNPs <- SNPs[,snps]
  
  expect_equal(SNPs, m)
})


# Multiple maskings ----
test_that('We can mask multiple patterns', {
  SNPs <- Siccuracy:::make.true(9, 12)
  snpfile <- tempfile()
  write.snps(SNPs, snpfile)
  
  masking = list(
    list(1:3, c(1,3,5,7,9,11)),
    list(4:6, c(2,4,6,8,10,12)),
    list(7:9, c(1,3,9,10,12))
  )
  
  fn <- tempfile()
  res <- mask_snp_file(snpfile, fn, masking=masking, na=9)
  m <- read.snps(fn)
  
  for (i in 1:length(masking)) {
    s <- masking[[i]]
    SNPs[s[[1]],s[[2]]] <- 9
  }
  
  expect_equal(SNPs, m)
})

# Using SnpsInNew ----
test_that('Using `snpsinnew`', {
  SNPs <- Siccuracy:::make.true(10, 14)
  snpfile <- tempfile()
  write.snps(SNPs, snpfile)
  
  snps <- sample.int(14, 8)
  
  masking = list(
    list(1:5, 1:4),
    list(6:10, 5:8)
  )
  
  fn <- tempfile()
  res <- mask_snp_file(snpfile, fn, snps=snps, masking=masking, snpsinnew = TRUE, na=9)
  m <- read.snps(fn)
  
  SNPs <- SNPs[,snps]
  for (i in 1:length(masking)) {
    s <- masking[[i]]
    SNPs[s[[1]],s[[2]]] <- 9
  }
  expect_equal(SNPs, m)
  
})


# Setting empty columns ----
test_that('Setting empty columns', {
  SNPs <- Siccuracy:::make.true(8, 5)
  snpfile <- tempfile()
  write.snps(SNPs, snpfile)
  
  snps <- c(1,0,2,0,3,0,4,0,5)
  
  fn <- tempfile()
  res <- mask_snp_file(snpfile, masking=5, fn, snps=snps, na=9)
  m <- read.snps(fn)
  
  expect_true(all(m[,c(2,4,6,8,9)] == 9))
  expect_true(all(m[,-c(2,4,6,8,9)] != 9))
  expect_equal(m[,c(1,3,5,7)], SNPs[,-5])
})
