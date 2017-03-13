library(testthat)
library(Siccuracy)

context('Imputation accuracies')

test_that("Lengths of returned vectors matches expected lengths (non-adaptive)",{
  ts <- Siccuracy:::make.test(15, 21)
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized=TRUE, adaptive=FALSE)
  true <- ts$true
  expect_equal(length(results$snps$means), ncol(true))
  expect_equal(length(results$snps$sds), ncol(true))
  expect_equal(length(results$animals$cors), nrow(true))
  expect_equal(length(results$matcor), 1)
  expect_equal(length(results$snps$cors), ncol(true))
  expect_equal(length(results$animals$rowID), nrow(true))
  expect_equal(results$animals$rowID, as.integer(rownames(true)))
})
test_that("Lengths of returned vectors matches expected lengths (adaptive)",{
  ts <- Siccuracy:::make.test(15, 21)
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized=TRUE, adaptive=TRUE)
  true <- ts$true
  expect_equal(length(results$snps$means), ncol(true))
  expect_equal(length(results$snps$sds), ncol(true))
  expect_equal(length(results$animals$cors), nrow(true))
  expect_equal(length(results$matcor), 1)
  expect_equal(length(results$snps$cors), ncol(true))
  expect_equal(length(results$animals$rowID), nrow(true))
  expect_equal(results$animals$rowID, as.integer(rownames(true)))
})

test_that("Results matches R's correlations (standardized=FALSE)",{
  ts <- Siccuracy:::make.test(15, 21)
  true <- ts$true
  imputed <- ts$imputed
  
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')

  #context('Fast')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, na=9, standardized=FALSE, adaptive=FALSE)
  expect_equal(results$matcor, mat1, tolerance=1e-15)
  expect_equal(results$animals$cors, row1, tolerance=1e-15)
  expect_equal(results$snps$cors, col1, tolerance=1e-15)
  expect_equal(results$animals$rowID, as.integer(rownames(true)))
  
  #context('Adaptive')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, na=9, standardized=FALSE, adaptive=TRUE)
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  expect_equal(results$animals$rowID, as.integer(rownames(true)))
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
  
  #context('Fast')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, na=9, standardized=TRUE, adaptive=FALSE)
  expect_equal(results$matcor, mat1, tolerance=1e-15)
  expect_equal(results$animals$cors, row1, tolerance=1e-15)
  expect_equal(results$snps$cors, col1, tolerance=1e-15)
  expect_equal(results$snps$means, m, tolerance=1e-9)
  expect_equal(results$snps$sds, v, tolerance=1e-9)
  
  
  #context('Adaptive')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, na=9, standardized=TRUE, adaptive=TRUE)
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  expect_equal(results$snps$means, m, tolerance=1e-9)
  expect_equal(results$snps$sds, v, tolerance=1e-9)
})

# Gene dosages ----

test_that('Imputation accuracies handles gene dosages, i.e. numeric values', {
  ts <- Siccuracy:::make.test(15, 21)
  r <- sample.int(prod(dim(ts$imputed)), prod(dim(ts$imputed))*0.5)
  imputed <- ts$imputed
  imputed[r] <- imputed[r] + round(rnorm(length(r), sd=0.3), 2)
  imputed[imputed < 0] <- 0
  imputed[imputed > 2] <- 2
  write.snps(imputed, ts$imputedfn)
  true <- ts$true
    
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  suppressWarnings(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')))
  
  results <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive=FALSE)
  expect_equal(results$matcor, mat1, tolerance=1e-8)
  expect_equal(results$animals$cors, row1, tolerance=1e-8)
  expect_equal(results$snps$cors, col1, tolerance=1e-8)

  results <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive=TRUE)
  expect_equal(results$matcor, mat1, tolerance=1e-8)
  expect_equal(results$animals$cors, row1, tolerance=1e-8)
  expect_equal(results$snps$cors, col1, tolerance=1e-8)
  
  m <- apply(true, 2, mean)
  v <- apply(true, 2, sd)
  true <- scale(true, m, v)
  imputed <- scale(imputed, m, v)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  results <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = TRUE, adaptive=FALSE)
  expect_equal(results$matcor, mat1, tolerance=1e-8)
  expect_equal(results$animals$cors, row1, tolerance=1e-8)
  expect_equal(results$snps$cors, col1, tolerance=1e-8)
  
  results <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = TRUE, adaptive=TRUE)
  expect_equal(results$matcor, mat1, tolerance=1e-8)
  expect_equal(results$animals$cors, row1, tolerance=1e-8)
  expect_equal(results$snps$cors, col1, tolerance=1e-8)
  
})

test_that('Non-adaptive handles missing SNPs in true files (exact match btw. true and genotyped)',{
  ts <- Siccuracy:::make.test(15, 21)
  imputed <- ts$imputed

  true <- ts$imputed
  true[,3] <- NA 
  true[1,4] <- NA
  true[1:2,5] <- NA
  write.snps(true, ts$truefn)
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = FALSE, adaptive = FALSE)
  
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)

  #context('Standardized')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = FALSE)
  
  m <- apply(true, 2, mean, na.rm=TRUE)
  v <- apply(true, 2, sd, na.rm=TRUE)
  true <- scale(true, m, v)
  imputed <- scale(imputed, m, v)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  v[is.na(v)] <- 0.0
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  expect_equal(results$snps$means, m)
  expect_equal(results$snps$sds, v)
  
})

test_that('Non-adaptive handles missing SNPs in true files',{
  ts <- Siccuracy:::make.test(15, 21)
  imputed <- ts$imputed
  
  true <- ts$true
  true[,3] <- NA 
  true[1,4] <- NA
  true[1:2,5] <- NA
  write.snps(true, ts$truefn)
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = FALSE, adaptive = FALSE)
  
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  
  #context('Standardized')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = FALSE)
  
  m <- apply(true, 2, mean, na.rm=TRUE)
  v <- apply(true, 2, sd, na.rm=TRUE)
  true <- scale(true, m, v)
  imputed <- scale(imputed, m, v)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  v[is.na(v)] <- 0.0
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  expect_equal(results$snps$means, m)
  expect_equal(results$snps$sds, v)
  
})

test_that('Adaptive handles missing SNPs in true files (exact match btw. true and genotyped)',{
  ts <- Siccuracy:::make.test(15, 21)
  imputed <- ts$imputed
  
  true <- ts$imputed
  true[,3] <- NA 
  true[1,4] <- NA
  true[1:2,5] <- NA
  write.snps(true, ts$truefn)
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = FALSE, adaptive = TRUE)
  
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  
  #context('Standardized')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = TRUE)
  
  m <- apply(true, 2, mean, na.rm=TRUE)
  v <- apply(true, 2, sd, na.rm=TRUE)
  true <- scale(true, m, v)
  imputed <- scale(imputed, m, v)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  v[is.na(v)] <- 0.0
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  expect_equal(results$snps$means, m)
  expect_equal(results$snps$sds, v)
  
})

# Adaptive handles missing SNPs in true files ----
test_that('Adaptive handles missing SNPs in true files',{
  ts <- Siccuracy:::make.test(15, 21)
  imputed <- ts$imputed
  
  true <- ts$true
  true[,3] <- NA 
  true[1,4] <- NA
  true[1:2,5] <- NA
  write.snps(true, ts$truefn)
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = FALSE, adaptive = TRUE)
  
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  
  #context('Standardized')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = TRUE)
  
  m <- apply(true, 2, mean, na.rm=TRUE)
  v <- apply(true, 2, sd, na.rm=TRUE)
  true <- scale(true, m, v)
  imputed <- scale(imputed, m, v)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  v[is.na(v)] <- 0.0
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  expect_equal(results$snps$means, m)
  expect_equal(results$snps$sds, v)
  
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
  expect_warning(col2 <- sapply(1:ncol(true), function(i) cor(true[r,i], imputed[,i], use='na.or.complete'))  , regexp = 'the standard deviation is zero')
    
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, na=9, standardized=FALSE, adaptive=TRUE)
  expect_equal(results$matcor, mat2, tolerance=1e-9)
  expect_equal(results$animals$cors, row2, tolerance=1e-9)
  expect_equal(results$snps$cors, col2, tolerance=1e-9)
  expect_equal(results$snps$means, rep(0, length(results$snps$means)), tolerance=1e-9)
  expect_equal(results$snps$sds, rep(1, length(results$snps$sds)), tolerance=1e-9)  
  expect_length(results$snps$means, ncol(true))
  expect_length(results$snps$sds, ncol(true))
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
  expect_warning(col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[r,i], use='na.or.complete'))  , regexp = 'the standard deviation is zero')
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, na=9, standardized=FALSE, adaptive=TRUE)
  expect_equal(results$matcor, mat2, tolerance=1e-9)
  expect_equal(results$animals$cors, row2, tolerance=1e-9)
  expect_equal(results$snps$cors, col2, tolerance=1e-9)
  expect_equal(results$snps$means, rep(0, length(results$snps$means)), tolerance=1e-9)
  expect_equal(results$snps$sds, rep(1, length(results$snps$sds)), tolerance=1e-9)  
  expect_length(results$snps$means, ncol(true))
  expect_length(results$snps$sds, ncol(true))
})

# Animal's genotype has no variance ---
test_that('Animal\'s genotype has no variance', {
  ts <- Siccuracy:::make.test(15, 21)
  true <- ts$true
  true[2,] <- 2
  write.snps(true, ts$truefn)
  
  # No standardization, as this changes each element of row 2 -- and it gets variance!
  imputed <- ts$imputed
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  suppressWarnings(row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete')))
  suppressWarnings(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')))
  
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive = TRUE)
  expect_equal(res$matcor, mat1, tolerance=1e-9)
  expect_equal(res$animals$cors, row1, tolerance=1e-9)
  expect_equal(res$snps$cors, col1, tolerance=1e-9)
  
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive = FALSE)
  expect_equal(res$matcor, mat1, tolerance=1e-9)
  expect_equal(res$animals$cors, row1, tolerance=1e-9)
  expect_equal(res$snps$cors, col1, tolerance=1e-9)
})

test_that('User-provided centering works',{
  ts <- Siccuracy:::make.test(31, 87)

  m <- runif(ncol(ts$true), 0, 2)
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = FALSE, center=m)

  true <- scale(ts$true, m, scale=FALSE)
  imputed <- scale(ts$imputed, m, scale=FALSE)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  
  #context('Adaptive')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = TRUE, center=m)
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  
})

test_that('User-provided scaling works',{
  ts <- Siccuracy:::make.test(31, 87)
  
  v <- runif(ncol(ts$true), 0.1, 2)
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = FALSE, scale=v)
  
  true <- scale(ts$true, center=FALSE, scale=v)
  imputed <- scale(ts$imputed, center=FALSE, scale=v)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  
  #context('Adaptive')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = TRUE, scale=v)
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  
})

# User-provided allele frequencies works: ----
test_that('User-provided allele frequencies works:',{
  ts <- Siccuracy:::make.test(10, 13)
  
  p <- seq(0.01, 0.05, length.out=ncol(ts$true))
  m <- 2*p
  v <- 2*p*(1-p)
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = FALSE, p=p)
  
  true <- scale(ts$true, center=m, scale=v)
  imputed <- scale(ts$imputed, center=m, scale=v)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  
  #context('Adaptive')
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, standardized = TRUE, adaptive = TRUE, p=p)
  
  expect_equal(results$matcor, mat1, tolerance=1e-9)
  expect_equal(results$animals$cors, row1, tolerance=1e-9)
  expect_equal(results$snps$cors, col1, tolerance=1e-9)
  
})

# Adaptive and non-adaptive returns same results ----
test_that('Adaptive and non-adaptive returns same results',{
  ts <- Siccuracy:::make.test(47, 108)
  ids <- sample.int(nrow(ts$true), nrow(ts$true) * 0.6)
  
  true <- ts$true
  true[sample.int(prod(dim(true)), prod(dim(true))*0.1)] <- NA
  t1 <- tempfile()
  t2 <- tempfile()
  write.snps(true, t1)
  write.snps(true[ids,], t2)
  
  imputed <- ts$imputed[ids,]
  write.snps(ts$imputed[ids,], ts$imputedfn)
  
  # standardized has to be disabled (or provided by other means), as the two true-files differs!
  resn <- imputation_accuracy(t2, ts$imputedfn, standardized=FALSE, adaptive=FALSE)
  resa <- imputation_accuracy(t1, ts$imputedfn, standardized=FALSE, adaptive=TRUE)
  expect_equal(resn$matcor, resa$matcor)
  expect_equal(resn$snps$cors, resa$snps$cors)
  expect_equal(resn$animals$cors, resa$animals$cors[match(resn$animals$rowID, resa$animals$rowID)])

  het <- heterozygosity(t1)  
  resn <- imputation_accuracy(t2, ts$imputedfn, standardized=TRUE, adaptive=FALSE, p=het$p)
  resa <- imputation_accuracy(t1, ts$imputedfn, standardized=TRUE, adaptive=TRUE, p=het$p)
  expect_equal(resn$matcor, resa$matcor)
  expect_equal(resn$snps$cors, resa$snps$cors)
  expect_equal(resn$animals$cors, resa$animals$cors[match(resn$animals$rowID, resa$animals$rowID)])
})
  
test_that('True-animals gets correlation of 1',{
  ts <- Siccuracy:::make.test(47, 108)
  tid <- 1:15
  imputed <- ts$imputed
  true <- ts$true
  imputed[tid,] <- true[tid,]
  write.snps(imputed, ts$imputedfn)
  
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')
  
  res1 <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive = FALSE)
  expect_equal(res1$snps$cors, col1)
  expect_equal(res1$animals$cors, row1)
  expect_equal(res1$matcor, mat1)
  expect_equal(res1$animals$cors[tid], rep(1, length(tid)))
  
  res1 <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive = TRUE)
  expect_equal(res1$snps$cors, col1)
  expect_equal(res1$animals$cors, row1)
  expect_equal(res1$matcor, mat1)
  expect_equal(res1$animals$cors[tid], rep(1, length(tid)))
  
  m <- apply(true, 2, mean)
  v <- apply(true, 2, sd)
  true <- scale(true, m, v)
  imputed <- scale(imputed, m, v)
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))

  res2 <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = TRUE, adaptive = FALSE)
  expect_equal(res2$snps$cors, col2)
  expect_equal(res2$animals$cors, row2)
  expect_equal(res2$matcor, mat2)
  expect_equal(res2$animals$cors[tid], rep(1, length(tid)))
  
  res2 <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = TRUE, adaptive = TRUE)
  expect_equal(res2$snps$cors, col2)
  expect_equal(res2$animals$cors, row2)
  expect_equal(res2$matcor, mat2)
  expect_equal(res2$animals$cors[tid], rep(1, length(tid)))
})

test_that('True-columns gets correlation of 1',{
  ts <- Siccuracy:::make.test(47, 108)
  imputed <- ts$imputed
  true <- ts$true
  cid <- sample.int(ncol(true),ncol(true)*0.15)
  imputed[,cid] <- true[,cid]
  write.snps(imputed, ts$imputedfn)
  
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')
  
  cid <- setdiff(cid, which(is.na(col1)))
  
  res1 <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive = FALSE)
  expect_equal(res1$snps$cors, col1)
  expect_equal(res1$animals$cors, row1)
  expect_equal(res1$matcor, mat1)
  expect_equal(res1$snps$cors[cid], rep(1, length(cid)))
  
  res1 <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive = TRUE)
  expect_equal(res1$snps$cors, col1)
  expect_equal(res1$animals$cors, row1)
  expect_equal(res1$matcor, mat1)
  expect_equal(res1$snps$cors[cid], rep(1, length(cid)))
  
  m <- apply(true, 2, mean)
  v <- apply(true, 2, sd)
  true <- scale(true, m, v)
  imputed <- scale(imputed, m, v)
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  res2 <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = TRUE, adaptive = FALSE)
  expect_equal(res2$snps$cors, col2)
  expect_equal(res2$animals$cors, row2)
  expect_equal(res2$matcor, mat2)
  expect_equal(res1$snps$cors[cid], rep(1, length(cid)))
  
  res2 <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = TRUE, adaptive = TRUE)
  expect_equal(res2$snps$cors, col2)
  expect_equal(res2$animals$cors, row2)
  expect_equal(res2$matcor, mat2)
  expect_equal(res1$snps$cors[cid], rep(1, length(cid)))
})

test_that('Excluding SNPs by given NA allele frequencies does or does not break',{
  ts <- Siccuracy:::make.test(15, 21)
  
  p <- rep(0.5, ncol(ts$true))
  p[4] <- 0
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, p=p)
  
  true <- scale(ts$true, 2*p, 2*p*(1-p))
  imputed <- scale(ts$imputed, 2*p, 2*p*(1-p))
  true <- true[,-4]
  imputed <- imputed[,-4]
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  suppressWarnings(col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')))

  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors[-4], col2)
  
  p[4] <- NA
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, p=p)
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors[-4], col2)
  
})

test_that('Excluding SNPs by given NA allele frequencies does or does not break, non-adaptive',{
  ts <- Siccuracy:::make.test(15, 21)
  
  p <- rep(0.5, ncol(ts$true))
  p[4] <- 0
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, p=p, adaptive = FALSE)
  
  true <- scale(ts$true, 2*p, 2*p*(1-p))
  imputed <- scale(ts$imputed, 2*p, 2*p*(1-p))
  true <- true[,-4]
  imputed <- imputed[,-4]
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors[-4], col2)
  
  p[4] <- NA
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, p=p, adaptive = FALSE)
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors[-4], col2)
  
})


test_that('Excluding individuals or SNPs from correations', {
  ts <- Siccuracy:::make.test(15, 21)
  noi <- c(3,8)
  nos <- c(2,9,10)
  
  
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=FALSE, excludeSNPs=nos)#, excludeIDs=noi)#, excludeSNPs=nos)
  
  true <- ts$true
  #true[noi,] <- NA
  true[,nos] <- NA
  imputed <- ts$imputed
  
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  suppressWarnings(col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')))
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors, col2)

  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=FALSE, excludeIDs=noi)#, excludeSNPs=nos)
  
  true <- ts$true
  true[noi,] <- NA
  #true[,nos] <- NA
  imputed <- ts$imputed
  
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors, col2)
  
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=FALSE, excludeIDs=noi, excludeSNPs=nos)
  
  true <- ts$true
  true[noi,] <- NA
  true[,nos] <- NA
  imputed <- ts$imputed
  
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  suppressWarnings(col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')))
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors, col2)
  
})

# Excluding individuals or SNPs from correations, adaptive ----
test_that('Excluding individuals or SNPs from correations, adaptive', {
  ts <- Siccuracy:::make.test(15, 21)
  noi <- c(3,8)
  nos <- c(2,9,10)
  
  write.snps(ts$imputed[sample.int(nrow(ts$imputed)),], ts$imputedfn)
  
  
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=TRUE, excludeSNPs=nos)#, excludeIDs=noi)#, excludeSNPs=nos)
  
  true <- ts$true
  #true[noi,] <- NA
  true[,nos] <- NA
  imputed <- ts$imputed
  
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  suppressWarnings(col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')))
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors, col2)
  
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=TRUE, excludeIDs=noi)#, excludeSNPs=nos)
  
  true <- ts$true
  true[noi,] <- NA
  #true[,nos] <- NA
  imputed <- ts$imputed
  #imputed[noi,] <- NA
  
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  expect_warning(col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')), regexp = 'the standard deviation is zero')
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors, col2)
    
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=TRUE, excludeIDs=noi, excludeSNPs=nos)
  
  true <- ts$true
  true[noi,] <- NA
  true[,nos] <- NA
  imputed <- ts$imputed
  
  mat2 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row2 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  suppressWarnings(col2 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete')))
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors, col2)
  
})


# Constant animal -----
test_that('Animal is constant', {
  ts <- Siccuracy:::make.test(15, 21)
  true <- ts$true
  
  # No standardization, as it induces variance in a row
  true[5,] <- 2
  true[,8] <- 2
  write.snps(true, ts$truefn)
  
  
  mat2 <- cor(as.vector(true), as.vector(ts$imputed), use = 'complete.obs')
  suppressWarnings(row2 <- sapply(1:nrow(true), function(i) cor(true[i,], ts$imputed[i,], use='na.or.complete')))
  suppressWarnings(col2 <- sapply(1:ncol(true), function(i) cor(true[,i], ts$imputed[,i], use='na.or.complete')))
  expect_true(is.na(row2[5]))
  expect_true(is.na(col2[8]))
  
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive=FALSE)
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors, col2)
  
  res <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized = FALSE, adaptive=TRUE)
  expect_equal(res$matcor, mat2)  
  expect_equal(res$animals$cors, row2)
  expect_equal(res$snps$cors, col2)
})

# Counting correct and incorrect ----
test_that('Counting correct and incorrect works, no dosages', {
  ts <- Siccuracy:::make.test(15, 21)
  true <- ts$true
  true[sample.int(prod(dim(true)), prod(dim(true))*0.4)] <- NA
  write.snps(true, ts$truefn)
  imputed <- ts$imputed
  
  rownames(imputed) <- rownames(true) <- NULL
  
  comp <- true == imputed
  true <- is.na(true)
  imputed <- is.na(imputed)
  both.na <- true & imputed
  only.tru <- true & !imputed
  only.imp <- !true & imputed
  
  row.correct <- rowSums(comp, na.rm = TRUE)
  row.na.imp <- rowSums(only.imp)
  row.na.tru <- rowSums(only.tru)
  row.na.both <- rowSums(both.na)
  
  col.correct <- colSums(comp, na.rm=TRUE)
  col.na.imp <- colSums(only.imp)
  col.na.tru <- colSums(only.tru)
  col.na.both <- colSums(both.na)

  results <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=FALSE)
  expect_equal(results$snps$correct, col.correct)
  expect_equal(results$snps$true.na, col.na.tru)
  expect_equal(results$snps$imp.na, col.na.imp)
  expect_equal(results$snps$both.na, col.na.both)

  expect_equal(results$animals$correct, row.correct)
  expect_equal(results$animals$true.na, row.na.tru)
  expect_equal(results$animals$imp.na, row.na.imp)
  expect_equal(results$animals$both.na, row.na.both)

  results <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=TRUE)
  expect_equal(results$snps$correct, col.correct)
  expect_equal(results$snps$true.na, col.na.tru)
  expect_equal(results$snps$imp.na, col.na.imp)
  expect_equal(results$snps$both.na, col.na.both)
  
  expect_equal(results$animals$correct, row.correct)
  expect_equal(results$animals$true.na, row.na.tru)
  expect_equal(results$animals$imp.na, row.na.imp)
  expect_equal(results$animals$both.na, row.na.both)
  
})

test_that('Counting correct and incorrect works, dosages', {
  ts <- Siccuracy:::make.test(15, 21)
  true <- ts$true
  true[sample.int(prod(dim(true)), prod(dim(true))*0.4)] <- NA
  write.snps(true, ts$truefn)
  imputed <- ts$imputed
  p <- sample.int(prod(dim(imputed)), prod(dim(imputed))*0.7)
  imputed[p] <- imputed[p] + rnorm(length(p), mean=0, sd=0.2)
  imputed[imputed > 2] <- 2
  imputed[imputed < 0] <- 0
  write.snps(imputed, ts$imputedfn)
  
  rownames(imputed) <- rownames(true) <- NULL
  
  comp <- abs(true - imputed) < 0.2
  true <- is.na(true)
  imputed <- is.na(imputed)
  both.na <- true & imputed
  only.tru <- true & !imputed
  only.imp <- !true & imputed

  row.correct <- rowSums(comp, na.rm = TRUE)
  row.na.imp <- rowSums(only.imp)
  row.na.tru <- rowSums(only.tru)
  row.na.both <- rowSums(both.na)
  
  col.correct <- colSums(comp, na.rm=TRUE)
  col.na.imp <- colSums(only.imp)
  col.na.tru <- colSums(only.tru)
  col.na.both <- colSums(both.na)
  
  results <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=FALSE, tol = 0.2)
  expect_equal(results$snps$correct, col.correct)
  expect_equal(results$snps$true.na, col.na.tru)
  expect_equal(results$snps$imp.na, col.na.imp)
  expect_equal(results$snps$both.na, col.na.both)
  
  expect_equal(results$animals$correct, row.correct)
  expect_equal(results$animals$true.na, row.na.tru)
  expect_equal(results$animals$imp.na, row.na.imp)
  expect_equal(results$animals$both.na, row.na.both)

  results <- imputation_accuracy(ts$truefn, ts$imputedfn, standardized=FALSE, adaptive=TRUE, tol = 0.2)
  expect_equal(results$snps$correct, col.correct)
  expect_equal(results$snps$true.na, col.na.tru)
  expect_equal(results$snps$imp.na, col.na.imp)
  expect_equal(results$snps$both.na, col.na.both)
  
  expect_equal(results$animals$correct, row.correct)
  expect_equal(results$animals$true.na, row.na.tru)
  expect_equal(results$animals$imp.na, row.na.imp)
  expect_equal(results$animals$both.na, row.na.both)
  
    
})

# What happens when an ID is repeated in the imputed file? -----
test_that('Repeated IDs are handled somehow?', {
  ts <- Siccuracy:::make.test(5, 13)
  i2 <- Siccuracy:::make.imputed(ts$true)
  
  write.snps(i2, ts$imputedfn, append=TRUE)
  
  
  i0 <- read.snps(ts$imputedfn)
  t1 <- ts$true
  i1 <- ts$imputed
  
  # Works as expected in non-adaptive because there are fewer true rows.
  res1 <- imputation_accuracy(ts$truefn, ts$imputedfn, adaptive = FALSE, standardized = FALSE)
  
  mat2 <- cor(as.vector(t1), as.vector(i1), use = 'complete.obs')
  row2 <- sapply(1:nrow(t1), function(i) cor(t1[i,], i1[i,], use='na.or.complete'))
  suppressWarnings(col2 <- sapply(1:ncol(t1), function(i) cor(t1[,i], i1[,i], use='na.or.complete')))
  
  expect_equal(res1$matcor, mat2)
  expect_equal(res1$snps$cors, col2)
  expect_equal(res1$animals$cors, row2)
  
  # Now, what happens with adaptive???
  
  expect_warning(res2 <- imputation_accuracy(ts$truefn, ts$imputedfn, adaptive = TRUE, standardized = FALSE))
  
  
  # If both files have repeated IDs
  tf2 <- tempfile()
  write.snps(t1, tf2)
  write.snps(t1, tf2, append=TRUE)
  res3 <- imputation_accuracy(tf2, ts$imputedfn, adaptive = FALSE, standardized = FALSE)
  
  res4 <- imputation_accuracy(tf2, ts$imputedfn, adaptive = TRUE, standardized = FALSE)
  expect_equal(res3, res4)
  
  
})

