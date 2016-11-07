library(Siccuracy)
library(testthat)

context('Testing merging SNP chips.')

test_that('Two, partially overlapping SNP chips can be merged.', {
  full.snps <- Siccuracy:::make.true(5, 18)
  rownames(full.snps) <- 1:nrow(full.snps)
  hd.snps <- c(1:3,5:9,11:12)
  ld.snps <- setdiff(1:12, c(2,3,5,8,10))
  
  hd.id <- c(1,3)
  ld.id <- c(2,4,5)
  
  hd <- tempfile()
  write.snps(full.snps[hd.id, hd.snps, drop=FALSE], hd)
  write.snps(full.snps[2, hd.snps, drop=FALSE], hd, append=TRUE)
  
  ld <- tempfile()
  full.snps2 <- full.snps
  full.snps2[5,1] <- 2
  write.snps(full.snps2[ld.id, ld.snps, drop=FALSE], ld)
  
  mergedfn <- tempfile()
  #mergedfn <- 'merged.txt'
  res <- rbind_SNPs(hdid=c(0, hd.id, 7), ldid=ld.id, hdpos=hd.snps, ldpos=ld.snps, hdfn=hd, ldfn=ld, fnout=mergedfn, na=9)
  
  expect_equivalent(res, length(hd.id)+length(ld.id))
  
  merged <- read.snps(mergedfn)
  
  # Test collected
  unique.cols <- sort(unique(c(hd.snps, ld.snps)))
  animals <- c(hd.id, ld.id)
  m <- full.snps[animals, 1:max(unique.cols), drop=FALSE]
  
  m[as.character(hd.id),setdiff(1:max(unique.cols),hd.snps)] <- 9
  m[as.character(ld.id),setdiff(1:max(unique.cols),ld.snps)] <- 9
  expect_equal(merged, m)
})

test_that('Order of ID\' doesn\'t matter.', {
  full.snps <- Siccuracy:::make.true(20, 30)
  hd.snps <- sample.int(ncol(full.snps), ncol(full.snps)*0.65)
  ld.snps <- sample.int(ncol(full.snps), ncol(full.snps)*0.20)
  
  hd.id <- sample.int(nrow(full.snps), nrow(full.snps)*0.5)
  ld.id <- setdiff(1:nrow(full.snps), hd.id)
  
  hd <- tempfile()
  ld <- tempfile()
  write.snps(full.snps[hd.id, hd.snps, drop=FALSE], hd)
  write.snps(full.snps[ld.id, ld.snps, drop=FALSE], ld)
  
  mergedfn <- tempfile()
  #mergedfn <- 'merged.txt'
  res <- rbind_SNPs(hdid=hd.id, ldid=ld.id, hdpos=hd.snps, ldpos=ld.snps, hdfn=hd, ldfn=ld, fnout=mergedfn, na=9)
  
  expect_equivalent(res, length(hd.id)+length(ld.id))
  
  merged <- read.snps(mergedfn)
  
  expect_equal(merged[1:length(hd.id),hd.snps],full.snps[hd.id, hd.snps])
  expect_equal(merged[(length(hd.id)+1):nrow(merged),ld.snps],full.snps[ld.id, ld.snps])

  expect_true(all(as.integer(merged[1:length(hd.id),setdiff(1:ncol(merged), hd.snps)]) == 9))
  expect_true(all(as.integer(merged[(length(hd.id)+1):nrow(merged),setdiff(1:ncol(merged), ld.snps)]) == 9))
})

test_that('Merging doesn\'t skip first row when starting in middle of ID list.', {
  hd <- 200:210
  ld <- 210:220
  both <- c(hd[3:7],ld[3:7])
  snps <- matrix(1, ncol=1, nrow=length(both), dimnames=list(both, 'a'))
  hd.fn <- tempfile()
  ld.fn <- tempfile()
  write.snps(snps, hd.fn)
  write.snps(snps, ld.fn)
  merged.fn <- tempfile()
  res <- rbind_SNPs(hdid=hd, ldid=ld, hdpos=1, ldpos=1, hdfn=hd.fn, ldfn=ld.fn, fnout=merged.fn)
  expect_equivalent(res, length(both))
  
  merged <- read.snps(merged.fn)
  expect_true(all(names(merged) %in% as.character(both)))
})

test_that('Searching exists when row not found in ID list.', {
  snps <- matrix(1, ncol=1, nrow=2)
  hd.fn <- tempfile()
  write.snps(snps, hd.fn)
  
  merged.fn <- tempfile()
  res <- rbind_SNPs(hdid=sample.int(10)+3, ldid=sample.int(10)+3, hdpos=1, ldpos=1, hdfn=hd.fn, ldfn=hd.fn, fnout=merged.fn)
  expect_equivalent(res, 0)
  merged <- read.snps(merged.fn)
  expect_equal(nrow(merged), 0)
  expect_equal(ncol(merged), 0)
})


