context('Imputation accuracies')

# make.true <- function(n,m) {
#   true <- matrix(sample(0:2, size = n*m, replace=TRUE), ncol=m)  # fill true with random 0, 1, or 2.
#   # add non-segregating site
#   i <- sample.int(m, 1)
#   true[2:n,i] <- true[1,i]
#   true
# }
# make.imputed <- function(true) {
#   m <- ncol(true)
#   n <- nrow(true)
#   imputed <- true
#   imputed[sample.int(n*m, floor(n*m*0.5))] <- sample(0:2, size=floor(n*m*0.5), replace=TRUE) # change half the elements
#   imputed[sample.int(n*m, floor(n*m*0.1))] <- NA  # some elements are missing.
#   imputed
# }

test_that("Lengths of returned vectors matches expected lengths (fast)",{
  ts <- Siccuracy:::make.test(15, 21)
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized=TRUE, fast=TRUE)
  true <- ts$true
  expect_equal(length(results$means), ncol(true))
  expect_equal(length(results$vars), ncol(true))
  expect_equal(length(results$rowcors), nrow(true))
  expect_equal(length(results$matcor), 1)
  expect_equal(length(results$colcors), ncol(true))
  expect_equal(length(results$rowID), nrow(true))
  expect_equal(results$rowID, as.integer(rownames(true)))
})
test_that("Lengths of returned vectors matches expected lengths (fast)",{
  ts <- Siccuracy:::make.test(15, 21)
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized=TRUE, fast=FALSE)
  true <- ts$true
  expect_equal(length(results$means), ncol(true))
  expect_equal(length(results$vars), ncol(true))
  expect_equal(length(results$rowcors), nrow(true))
  expect_equal(length(results$matcor), 1)
  expect_equal(length(results$colcors), ncol(true))
  expect_equal(length(results$rowID), nrow(true))
  expect_equal(results$rowID, as.integer(rownames(true)))
})

test_that("Results matches R's correlations (standardized=FALSE)",{
  ts <- Siccuracy:::make.test(15, 21)
  true <- ts$true
  imputed <- ts$imputed
  
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))

  context('Fast')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, fast=TRUE)
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$rowcors, row1, tolerance=1e-9)
  expect_equal(results$colcors, col1, tolerance=1e-9)
  expect_equal(results$rowID, as.integer(rownames(true)))
  
  context('Adaptive')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, fast=FALSE)
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$rowcors, row1, tolerance=1e-9)
  expect_equal(results$colcors, col1, tolerance=1e-9)
  expect_equal(results$rowID, as.integer(rownames(true)))
})


test_that("Results matches R's correlations (standardized=TRUE)",{
  ts <- Siccuracy:::make.test(15, 21)
  true <- ts$true
  imputed <- ts$imputed
  
  m <- apply(true, 2, mean)
  v <- apply(true, 2, sd)
  true <- scale(true, m, v)
  imputed <- scale(imputed, m, v)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  context('Fast')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=TRUE, fast=TRUE)
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$rowcors, row1, tolerance=1e-9)
  expect_equal(results$colcors, col1, tolerance=1e-9)
  expect_equal(results$means, m, tolerance=1e-9)
  expect_equal(results$sds, v, tolerance=1e-9)
  
  
  context('Adaptive')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=TRUE, fast=FALSE)
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$rowcors, row1, tolerance=1e-9)
  expect_equal(results$colcors, col1, tolerance=1e-9)
  expect_equal(results$means, m, tolerance=1e-9)
  expect_equal(results$sds, v, tolerance=1e-9)
})

test_that("Adaptive works with more true rows, and imputed are shuffled",{
  ts <- Siccuracy:::make.test(15, 21)
  true <- ts$true
  # Shuffle and drop some rows of imputed.
  r <- sample.int(nrow(true),8)
  imputed <- ts$imputed[r,]
  write.snps(imputed, ts$imputedfn)
  
  mat2 <- cor(as.vector(true[r,]), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(r, function(i) cor(true[i,], ts$imputed[i,], use='na.or.complete'))[order(r)]
  col2 <- sapply(1:ncol(true), function(i) cor(true[r,i], imputed[,i], use='na.or.complete'))  
    
  # Cannot expect error from fast subroutine.
  #results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=TRUE, fast=TRUE)
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, fast=FALSE)
  expect_equal(results$matcor, mat2, tolerance=1e-9)
  expect_equal(results$rowcors, row2, tolerance=1e-9)
  expect_equal(results$colcors, col2, tolerance=1e-9)
  expect_equal(results$means, rep(0, length(results$means)), tolerance=1e-9)
  expect_equal(results$sds, rep(1, length(results$sds)), tolerance=1e-9)  
  expect_length(results$means, ncol(true))
  expect_length(results$sds, ncol(true))
})

test_that("Adaptive works with less true rows, and imputed are shuffled",{
  ts <- Siccuracy:::make.test(15, 21)
  imputed <- ts$imputed
  # Shuffle and drop some rows of imputed.
  r <- sample.int(nrow(imputed),8)
  true <- ts$true[r,]
  write.snps(true, ts$truefn)
  
  mat2 <- cor(as.vector(true), as.vector(imputed[r,]), use = 'complete.obs')
  row2 <- sapply(r, function(i) cor(ts$true[i,], ts$imputed[i,], use='na.or.complete'))
  col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[r,i], use='na.or.complete'))  
  
  # Cannot expect error from fast subroutine.
  #results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=TRUE, fast=TRUE)
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, fast=FALSE)
  expect_equal(results$matcor, mat2, tolerance=1e-9)
  expect_equal(results$rowcors, row2, tolerance=1e-9)
  expect_equal(results$colcors, col2, tolerance=1e-9)
  expect_equal(results$means, rep(0, length(results$means)), tolerance=1e-9)
  expect_equal(results$sds, rep(1, length(results$sds)), tolerance=1e-9)  
  expect_length(results$means, ncol(true))
  expect_length(results$sds, ncol(true))
})