library(Siccuracy)
library(testthat)

context('Testing merging SNP chips.')

test_that('Two, partially overlapping SNP chips can be merged.', {
  full.snps <- Siccuracy:::make.true(5, 18)
  rownames(full.snps) <- 1:nrow(full.snps)
  hd.snps <- c(1:3,5:9,11:12)
  ld.snps <- setdiff(1:12, c(2,3,5,8,10))
  
  hd.id <- c(1,3,5)
  ld.id <- c(2,4,5)
  
  hd <- tempfile()
  write.snps(full.snps[hd.id, hd.snps, drop=FALSE], hd)
  write.snps(full.snps[2, hd.snps, drop=FALSE], hd, append=TRUE)
  
  ld <- tempfile()
  full.snps2 <- full.snps
  full.snps2[5,1] <- 2
  write.snps(full.snps[ld.id, ld.snps, drop=FALSE], ld)
  
  mergedfn <- tempfile()
  #mergedfn <- 'merged.txt'
  res <- rbind_SNPs(hdid=c(0, hd.id, 7), ldid=ld.id, hdpos=hd.snps, ldpos=ld.snps, hdfn=hd, ldfn=ld, fnout=mergedfn, na=9)
  
  unique.cols <- sort(unique(c(hd.snps, ld.snps)))
  animals <- sort(unique(c(hd.id, ld.id)))
  ld2 <- setdiff(ld.id, hd.id)
  ld3 <- setdiff(hd.snps, ld.snps)
  hd2 <- setdiff(ld.snps, hd.snps)
  m <- full.snps[animals, unique.cols, drop=FALSE]
  m[ld2,match(ld3, unique.cols)] <- 9
  m[hd.id,match(hd2, unique.cols)] <- 9
  
  merged <- read.snps(mergedfn)
  
  m <- m[match(rownames(merged), rownames(m)),]
  expect_equal(merged, m)
})

test_that('Order of ID\' doesn\'t matter.', {
  hd <- sample.int(500, 10)
  ld <- c(sample.int(500, 10), sample(hd, 3))
  both <- unique(c(hd,ld))
  both <- sample(both, length(both))
  
  snps <- matrix(1, ncol=1, nrow=length(both), dimnames=list(both, 'a'))
  hd.fn <- tempfile()
  write.snps(snps[as.character(hd),,drop=FALSE], hd.fn)
  ld.fn <- tempfile()
  write.snps(snps[as.character(ld),,drop=FALSE], ld.fn)
  merged.fn <- tempfile()
  res <- rbind_SNPs(hdid=hd, ldid=ld, hdpos=1, ldpos=1, hdfn=hd.fn, ldfn=ld.fn, fnout=merged.fn)
  
  merged <- read.snps(merged.fn, what=numeric(), quiet=TRUE)
  expect_true(all(names(merged) %in% as.character(both)))
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
  merged <- read.snps(merged.fn, what=numeric(), quiet=TRUE)
  expect_true(all(names(merged) %in% as.character(both)))
})

test_that('Searching exists when row not found in ID list.', {
  hd <- sample.int(500, 10)
  ld <- c(sample.int(500, 10), sample(hd, 3))
  both <- unique(c(hd,ld))
  both <- sample(both, length(both))
  
  snps <- matrix(1, ncol=1, nrow=length(both), dimnames=list(both, 'a'))
  hd.fn <- tempfile()
  write.snps(snps[as.character(setdiff(ld, hd)),,drop=FALSE], hd.fn)
  ld.fn <- tempfile()
  write.snps(snps[as.character(hd),,drop=FALSE], ld.fn)
  merged.fn <- tempfile()
  res <- rbind_SNPs(hdid=hd, ldid=ld, hdpos=1, ldpos=1, hdfn=hd.fn, ldfn=ld.fn, fnout=merged.fn)
  merged <- read.snps(merged.fn, what=numeric(), quiet=TRUE)
  expect_equal(nrow(merged), 0)
  #expect_true(all(names(merged) %in% as.character(both)))
  
  
  write.snps(snps[as.character(c(ld,hd)),,drop=FALSE], hd.fn)
  write.snps(snps[as.character(hd,ld),,drop=FALSE], ld.fn)
  res <- rbind_SNPs(hdid=hd, ldid=ld, hdpos=1, ldpos=1, hdfn=hd.fn, ldfn=ld.fn, fnout=merged.fn)
  merged <- read.snps(merged.fn, what=numeric(), quiet=TRUE)
  expect_true(all(names(merged) %in% as.character(both)))
  
})

test_that('Order of columns can be randomised.', {
  snps <- matrix(rep(1:8, each=12), ncol=8, nrow=12, dimnames = list(1:12, NULL))
  hdid <- c(4,9,3,1,2,10,11,12)
  ldid <- setdiff(1:12, hdid)
  hd.pos <- c(4, 5, 1, 8, 2, 3, 7, 6)  # sample.int(9,9)
  ld.pos <- c(6,3,4,8)
  
  hdfn <- tempfile()
  write.snps(snps[hdid, hd.pos, drop=FALSE], hdfn)
  ldfn <- tempfile()
  write.snps(snps[ldid, ld.pos, drop=FALSE], ldfn)
  
  merged <- tempfile()
  res <- rbind_SNPs(hdid=hdid, ldid=ldid, hdpos=hd.pos, ldpos=ld.pos, hdfn=hdfn, ldfn=ldfn, fnout=merged, na=9)
  merged <- read.snps(merged, what=numeric(), na = 9)
  
  for (i in 1:ncol(snps)) {
    expect_equal(merged[as.character(hdid),i], rep(i, length(hdid)), info=paste('hdid vs. pos.', i), check.attributes=FALSE)
    if (i %in% ld.pos)
      expect_equal(merged[as.character(ldid),i], rep(i, length(ldid)), info=paste('ldid vs. pos.', i), check.attributes=FALSE)
    if (!i %in% ld.pos)
      expect_true(all(is.na(merged[as.character(ldid),i])), info=paste('ld is NA on pos.', i))
  }
})
