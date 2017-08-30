library(testthat)
library(Siccuracy)

test_that('Oxford format is correctly parsed into 3D array', {
  
  # make test data
  tfn <- tempfile()
  gens <- matrix(as.integer(paste0(rep(1:5, times=6), rep(c(1,1,1,2,2,2),each=5), rep(c(1,2,3), each=5))), ncol=6)
  tst <- data.frame(1, letters[1:5], 1:5, 'A','C')
  write.table(cbind(tst, gens), tfn, row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  
  oxford <- read.oxford(tfn)
  
  cols <- (1:ncol(gens)) %% 3
  expect_true(all(gens[,cols == 1] %% 1 == 0))
  expect_true(all(gens[,cols == 2] %% 2 == 0))
  #expect_true(all(gens[,cols == 0] %% 3 == 0))

  expect_equivalent(gens[,cols == 1], oxford$probs[,,'hom1'])
  expect_equivalent(gens[,cols == 2], oxford$probs[,,'het'])
  expect_equivalent(gens[,cols == 0], oxford$probs[,,'hom2'])
  
})

