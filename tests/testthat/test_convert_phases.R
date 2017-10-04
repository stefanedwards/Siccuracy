library(testthat)
library(Siccuracy)

context('Testing converting phase files to genotype files')

test_that('Converting phase-files to genotype-files works for integers', {
  phase <- Siccuracy:::make.phase(4,12)
  phasefn <- tempfile()
  write.snps(phase, phasefn)
  # slow summation
  true <- matrix(0, nrow(phase)/2, ncol=ncol(phase), dimnames = list(unique(rownames(phase)), NULL))
  for (i in seq.int(2,nrow(phase),2)) {
    true[i/2,] <- phase[i-1,] + phase[i,]
  }
  
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn)
  b <- read.snps(truefn)
  expect_equal(n, nrow(true))
  expect_equal(b, true)
})


test_that('Converting phase-files to genotype-files works for numeric values', {
  phase <- Siccuracy:::make.phase(4,12)
  phase[] <- phase[] + round(runif(prod(dim(phase)), max=0.5), 2)
  phasefn <- tempfile()
  write.snps(phase, phasefn)
  # slow summation
  true <- matrix(0, nrow(phase)/2, ncol=ncol(phase), dimnames = list(unique(rownames(phase)), NULL))
  for (i in seq.int(2,nrow(phase),2)) {
    true[i/2,] <- phase[i-1,] + phase[i,]
  }
  true[true > 2] <- NA
  true[true < 0] <- NA
  
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn)
  b <- read.snps(truefn, what = numeric(), na=9)
  expect_equal(n, nrow(true))
  expect_equal(b, true)
})

test_that('Converting phase-files with missing values (integers)', {
  phase <- Siccuracy:::make.phase(4,12)
  phasefn <- tempfile()
  write.snps(phase, phasefn)
  # slow summation
  true <- matrix(0, nrow(phase)/2, ncol=ncol(phase), dimnames = list(unique(rownames(phase)), NULL))
  for (i in seq.int(2,nrow(phase),2)) {
    true[i/2,] <- phase[i-1,] + phase[i,]
  }
  true[true > 2] <- NA
  true[true < 0] <- NA
  
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn, na=9)
  b <- read.snps(truefn, what = integer(), na=9)
  expect_equal(n, nrow(true))
  expect_equal(b, true)
})

test_that('Converting phase-files with missing values (numeric)', {
  phase <- Siccuracy:::make.phase(4,12)
  phase[] <- phase[] + round(runif(prod(dim(phase)), max=0.5), 2)
  phasefn <- tempfile()
  write.snps(phase, phasefn)
  # slow summation
  true <- matrix(0, nrow(phase)/2, ncol=ncol(phase), dimnames = list(unique(rownames(phase)), NULL))
  for (i in seq.int(2,nrow(phase),2)) {
    true[i/2,] <- phase[i-1,] + phase[i,]
  }
  true[true > 2] <- NA
  true[true < 0] <- NA
  
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn, na=9)
  b <- read.snps(truefn, what = numeric(), na=9)
  expect_equal(n, nrow(true))
  expect_equal(b, true)
})

test_that('Maximum number of nlines', {
  phase <- Siccuracy:::make.phase(8,12)
  phasefn <- tempfile()
  write.snps(phase, phasefn)
  # slow summation
  true <- matrix(0, nrow(phase)/2, ncol=ncol(phase), dimnames = list(unique(rownames(phase)), NULL))
  for (i in seq.int(2,nrow(phase),2)) {
    true[i/2,] <- phase[i-1,] + phase[i,]
  }
  
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn, nlines=3)
  b <- read.snps(truefn)
  expect_equal(n, 3)
  expect_equal(b, true[1:3,])
})
