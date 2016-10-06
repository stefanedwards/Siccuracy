context('Imputation accuracies')

test_that("Lengths of returned vectors matches expected lengths (non-adaptive)",{
  ts <- Siccuracy:::make.test(15, 21)
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized=TRUE, adaptive=FALSE)
  true <- ts$true
  expect_equal(length(results$means), ncol(true))
  expect_equal(length(results$sds), ncol(true))
  expect_equal(length(results$rowcors), nrow(true))
  expect_equal(length(results$matcor), 1)
  expect_equal(length(results$colcors), ncol(true))
  expect_equal(length(results$rowID), nrow(true))
  expect_equal(results$rowID, as.integer(rownames(true)))
})
test_that("Lengths of returned vectors matches expected lengths (adaptive)",{
  ts <- Siccuracy:::make.test(15, 21)
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized=TRUE, adaptive=TRUE)
  true <- ts$true
  expect_equal(length(results$means), ncol(true))
  expect_equal(length(results$sds), ncol(true))
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
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, adaptive=FALSE)
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$rowcors, row1, tolerance=1e-9)
  expect_equal(results$colcors, col1, tolerance=1e-9)
  expect_equal(results$rowID, as.integer(rownames(true)))
  
  context('Adaptive')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, adaptive=TRUE)
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
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=TRUE, adaptive=FALSE)
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$rowcors, row1, tolerance=1e-9)
  expect_equal(results$colcors, col1, tolerance=1e-9)
  expect_equal(results$means, m, tolerance=1e-9)
  expect_equal(results$sds, v, tolerance=1e-9)
  
  
  context('Adaptive')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=TRUE, adaptive=TRUE)
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
    
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, adaptive=TRUE)
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
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, adaptive=TRUE)
  expect_equal(results$matcor, mat2, tolerance=1e-9)
  expect_equal(results$rowcors, row2, tolerance=1e-9)
  expect_equal(results$colcors, col2, tolerance=1e-9)
  expect_equal(results$means, rep(0, length(results$means)), tolerance=1e-9)
  expect_equal(results$sds, rep(1, length(results$sds)), tolerance=1e-9)  
  expect_length(results$means, ncol(true))
  expect_length(results$sds, ncol(true))
})