context('Basic imputation accuracies')

make.true <- function(n,m) {
  true <- matrix(sample(0:2, size = n*m, replace=TRUE), ncol=m)  # fill true with random 0, 1, or 2.
  # add non-segregating site
  i <- sample.int(m, 1)
  true[2:n,i] <- true[1,i]
  true
}
make.imputed <- function(true) {
  m <- ncol(true)
  n <- nrow(true)
  imputed <- true
  imputed[sample.int(n*m, floor(n*m*0.5))] <- sample(0:2, size=floor(n*m*0.5), replace=TRUE) # change half the elements
  imputed[sample.int(n*m, floor(n*m*0.1))] <- NA  # some elements are missing.
  imputed
}

test_that("Lengths of returned vectors matches expected lengths",{
  ts <- Siccuracy:::make.test(15, 21)
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized=TRUE, fast=TRUE)
  expect_equal(length(results$means), ncol(true))
  expect_equal(length(results$vars), ncol(true))
  expect_equal(length(results$rowcors), nrow(true))
  expect_equal(length(results$matcor), 1)
  expect_equal(length(results$colcors), ncol(true))
})

test_that("Results matches R's correlations (standardized=FALSE)",{
  true <- make.true(15,21); imputed <- make.imputed(true)
  truefn <- tempfile('true', fileext = '.txt')
  imputefn <- tempfile('imputed', fileext = '.txt')
  write.snps(true, truefn)
  write.snps(imputed, imputefn)
  results <- imputation_accuracy1(truefn=truefn, imputefn=imputefn, nSNPs=ncol(true), nAnimals=nrow(true), NAval=9, standardized=FALSE)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='complete.obs'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='complete.obs'))
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$rowcors, row1, tolerance=1e-9)
  expect_equal(results$colcors, col1, tolerance=1e-9)
})

test_that("Results matches R's correlations (standardized=TRUE)",{
  true <- make.true(21,15); imputed <- make.imputed(true)
  truefn <- tempfile('true', fileext = '.txt')
  imputefn <- tempfile('imputed', fileext = '.txt')
  write.snps(true, truefn)
  write.snps(imputed, imputefn)
  results <- imputation_accuracy1(truefn=truefn, imputefn=imputefn, nSNPs=ncol(true), nAnimals=nrow(true), NAval=9, standardized=TRUE)
  m <- apply(true, 2, mean)
  v <- apply(true, 2, var)
  #v[v==0] <- 1
  true <- scale(true, m, v)
  imputed <- scale(imputed, m, v)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='complete.obs'))
  col1 <- rep(NA, ncol(true))
  col1[v!=0] <- sapply(which(v!=0), function(i) cor(true[,i], imputed[,i], use='complete.obs'))
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$rowcors, row1, tolerance=1e-9)
  expect_equal(results$colcors, col1, tolerance=1e-9)
})
