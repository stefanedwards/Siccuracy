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

test_that('Converting phase-files respects numerical format', {
  phase <- Siccuracy:::make.imputed(Siccuracy:::make.true(2, 3))
  phase[] <- phase[] + round(runif(prod(dim(phase))), 3)  
  phasefn <- tempfile()
  write.snps(phase, phasefn, row.names = TRUE)
  truefn <- tempfile()
  n <- convert_phases(phasefn, truefn, int=FALSE, naval=9, numeric.format='8.4')
  s <- scan(truefn, what='character', quiet=TRUE)[-1]
  expect(all(nchar(s) == 6), 'Wrong number of characters in written output.')  # spaces are stripped; expected character length is 4+2 (decimals + period + leading digit).
  n <- convert_phases(phasefn, truefn, int=FALSE, naval=9, numeric.format='10.6')
  s <- scan(truefn, what='character', quiet=TRUE)[-1]
  expect(all(nchar(s) == 8), 'Wrong number of characters in written output.')  # spaces are stripped; expected character length is 6+2 (decimals + period + leading digit).

  context('Test result when number of decimals exceed width')
  n <- convert_phases(phasefn, truefn, int=FALSE, naval=9, numeric.format='10.8') 
  s <- scan(truefn, what='character', quiet=TRUE)
  expect(length(s) == 1, 'Space as separators are not missing.')
})