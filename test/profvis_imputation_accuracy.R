# propfing with profvis
# Guide: https://support.rstudio.com/hc/en-us/articles/218221837-Profiling-with-RStudio

# We examine the effect of different real types in the Fortran code.
# This is not straightforward as the code must be duplicated and the real-types updated;
# next we have to make the R wrapper functions accommondate that we want to do proofing.
# So far, I have added an argument 'full=FALSE', where full=TRUE indicates using maximum real type
# (i.e. double precision).
# In imputation_accuracy function, the subroutine is chosen as:
# if (full) {
#   subroutine <- ifelse(fast, 'imp_acc_fast_full','imp_acc')
# } else {
#   subroutine <- ifelse(fast, 'imp_acc_fast','imp_acc')
# }

library(profvis)

test_that('Using less than full double precision still gives same accuracy', {
  ts <- Siccuracy:::make.test(15, 21)
  true <- ts$true
  imputed <- ts$imputed
  
  results <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, fast=TRUE, full=FALSE)
  mat1 <- cor(as.vector(true), as.vector(imputed), use = 'complete.obs')
  row1 <- sapply(1:nrow(true), function(i) cor(true[i,], imputed[i,], use='na.or.complete'))
  col1 <- sapply(1:ncol(true), function(i) cor(true[,i], imputed[,i], use='na.or.complete'))
  # Note how we lessen the tolerance...
  expect_equal(results$matcor, mat1, tolerance=1e-7)
  expect_equal(results$rowcors, row1, tolerance=1e-7)
  expect_equal(results$colcors, col1, tolerance=1e-7)
})

profvis({
  ts <- Siccuracy:::make.test(1000, 10000)  
  results1 <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, fast=TRUE, full=FALSE)
  Sys.sleep(30)
  results2 <- imputation_accuracy(truefn=ts$truefn, imputefn=ts$imputedfn, NAval=9, standardized=FALSE, fast=TRUE, full=TRUE)
  expect_equal(results1$matcor, results2$matcor, tolerance=1e-6)
  expect_equal(results1$rowcors, results2$rowcors, tolerance=1e-6)
  expect_equal(results1$colcors, results2$colcors, tolerance=1e-6)
})

# Ultimateley, we cannot say anything about the memory usage using different real types,
# at least not in Windows.
