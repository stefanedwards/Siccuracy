context('Testing converting phase files to genotype files')

test_that('Converting phase-files to genotype-files works for integers', {
  phase <- Siccuracy:::make.true(8,12)
  phase[,1] <- rep(1:(nrow(phase)/2), each=2)
  phasefn <- tempfile()
  write.snps(phase, phasefn, row.names = FALSE)
  # slow summation
  true <- matrix(0, nrow(phase)/2, ncol=ncol(phase)-1, dimnames = list(1:(nrow(phase)/2), NULL))
  for (i in seq.int(2,nrow(phase),2)) {
    true[i/2,] <- phase[i-1,-1] + phase[i,-1]
  }
  
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn, int=TRUE)
  b <- read.snps(truefn)
  expect_equal(n, nrow(true))
  expect_equal(b, true)
})


test_that('Converting phase-files to genotype-files works for numeric values', {
  phase <- Siccuracy:::make.true(8,12)
  phase[] <- phase[] + round(runif(prod(dim(phase))), 2)
  phase[,1] <- rep(1:(nrow(phase)/2), each=2)
  phasefn <- tempfile()
  write.snps(phase, phasefn, row.names = FALSE)
  # slow summation
  true <- matrix(0, nrow(phase)/2, ncol=ncol(phase)-1, dimnames = list(1:(nrow(phase)/2), NULL))
  for (i in seq.int(2,nrow(phase),2)) {
    true[i/2,] <- phase[i-1,-1] + phase[i,-1]
  }
  
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn, int=FALSE)
  b <- read.snps(truefn, what = numeric())
  expect_equal(n, nrow(true))
  expect_equal(b, true)
})

test_that('Converting phase-files with missing values (integers)', {
  phase <- Siccuracy:::make.imputed(Siccuracy:::make.true(8, 12))
  phase[,1] <- rep(1:(nrow(phase)/2), each=2)
  phasefn <- tempfile()
  write.snps(phase, phasefn, row.names = FALSE)
  # slow summation
  true <- matrix(0, nrow(phase)/2, ncol=ncol(phase)-1, dimnames = list(1:(nrow(phase)/2), NULL))
  for (i in seq.int(2,nrow(phase),2)) {
    true[i/2,] <- phase[i-1,-1] + phase[i,-1]
  }
  
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn, int=TRUE, naval=9)
  b <- read.snps(truefn, what = integer(), na=9)
  expect_equal(n, nrow(true))
  expect_equal(b, true)
})

test_that('Converting phase-files with missing values (numeric)', {
  phase <- Siccuracy:::make.imputed(Siccuracy:::make.true(8, 12))
  phase[] <- phase[] + round(runif(prod(dim(phase))), 2)  
  phase[,1] <- rep(1:(nrow(phase)/2), each=2)
  phasefn <- tempfile()
  write.snps(phase, phasefn, row.names = FALSE)
  # slow summation
  true <- matrix(0, nrow(phase)/2, ncol=ncol(phase)-1, dimnames = list(1:(nrow(phase)/2), NULL))
  for (i in seq.int(2,nrow(phase),2)) {
    true[i/2,] <- phase[i-1,-1] + phase[i,-1]
  }
  
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn, int=FALSE, naval=9)
  b <- read.snps(truefn, what = numeric(), na=9)
  expect_equal(n, nrow(true))
  expect_equal(b, true)
})