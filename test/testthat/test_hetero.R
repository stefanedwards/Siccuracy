context('Testing calculation of heterozygosity')

test_that('Heterozygosity works with NA-values',{
  true <- Siccuracy:::make.true(9,18)
  true[sample(prod(dim(true)), 15)] <- NA
  truefn <- tempfile('true', fileext = '.txt')
  write.snps(true, truefn)
  res <- heterozygosity(truefn)

  p <- apply(true, 2, function(x) (sum(x==1, na.rm=TRUE) + sum(x==0, na.rm=TRUE)*2)/(2*sum(!is.na(x))))
  q <- apply(true, 2, function(x) (sum(x==1, na.rm=TRUE) + sum(x==2, na.rm=TRUE)*2)/(2*sum(!is.na(x))))

  Hobs <- apply(true, 2, function(x) sum(x==1, na.rm=TRUE)/sum(!is.na(x)))
  Hexp <- 2*p*q

  expect_equal(p, 1-q)
  expect_equal(res$p, p, tolerance=1e-7)
  expect_equal(1-res$p, q, tolerance=1e-7)
  expect_equal(res$n, apply(true, 2, function(x) sum(!is.na(x))))
  expect_equal(res$Hobs, Hobs, tolerance=1e-7)
  expect_equal(res$Hexp, Hexp, tolerance=1e-7)
})
