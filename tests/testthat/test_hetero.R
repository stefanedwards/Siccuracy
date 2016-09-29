context('Testing calculation of heterozygosity')

test_that('Heterozygosity works with a single population',{
  true <- Siccuracy:::make.true(9,18)
  truefn <- tempfile('true', fileext = '.txt')
  write.snps(true, truefn)

  p <- apply(true, 2, function(x) (sum(x==1, na.rm=TRUE) + sum(x==0, na.rm=TRUE)*2)/(2*sum(!is.na(x))))
  q <- apply(true, 2, function(x) (sum(x==1, na.rm=TRUE) + sum(x==2, na.rm=TRUE)*2)/(2*sum(!is.na(x))))
  
  Hobs <- apply(true, 2, function(x) sum(x==1, na.rm=TRUE)/sum(!is.na(x)))
  Hexp <- 2*p*q
  
  res <- heterozygosity(truefn)
  expect_equal(nrow(res), ncol(true))
  expect_equal(p, 1-q)
  expect_equal(res$p, p, tolerance=1e-7)
  expect_equal(1-res$p, q, tolerance=1e-7)
  expect_equal(res$n, apply(true, 2, function(x) sum(!is.na(x))))
  expect_equal(res$Hobs, Hobs, tolerance=1e-7)
  expect_equal(res$Hexp, Hexp, tolerance=1e-7)
  expect_equal(res$n, rep(nrow(true),ncol(true)))
})

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


test_that('Heterozygosity works with a two populations and missing values',{
  true <- Siccuracy:::make.true(9,4)
  true[sample(prod(dim(true)), 8)] <- NA
  truefn <- tempfile('true', fileext = '.txt')
  write.snps(true, truefn)
  
  populations <- rep(letters[1:2], each=5)[1:nrow(true)]
  
  pa <- apply(true[populations=='a',], 2, function(x) (sum(x==1, na.rm=TRUE) + sum(x==0, na.rm=TRUE)*2)/(2*sum(!is.na(x))))
  qa <- apply(true[populations=='a',], 2, function(x) (sum(x==1, na.rm=TRUE) + sum(x==2, na.rm=TRUE)*2)/(2*sum(!is.na(x))))
  pb <- apply(true[populations=='b',], 2, function(x) (sum(x==1, na.rm=TRUE) + sum(x==0, na.rm=TRUE)*2)/(2*sum(!is.na(x))))
  qb <- apply(true[populations=='b',], 2, function(x) (sum(x==1, na.rm=TRUE) + sum(x==2, na.rm=TRUE)*2)/(2*sum(!is.na(x))))
  na <- apply(true[populations=='a',], 2, function(x) sum(!is.na(x)))
  nb <- apply(true[populations=='b',], 2, function(x) sum(!is.na(x)))

  expect_equal(pa, 1-qa)
  expect_equal(pb, 1-qb)
  
  Hobsa <- apply(true[populations=='a',], 2, function(x) sum(x==1, na.rm=TRUE)/sum(!is.na(x)))
  Hobsb <- apply(true[populations=='b',], 2, function(x) sum(x==1, na.rm=TRUE)/sum(!is.na(x)))
  Hexpa <- 2*pa*qa
  Hexpb <- 2*pb*qb
  
  res <- heterozygosity(truefn, population = populations)
  expect_equal(nrow(res), 2*ncol(true)) # two populations in this
  expect_equal(with(res, p[populations=='a']), pa, tolerance=1e-7)
  expect_equal(with(res, p[populations=='b']), pb, tolerance=1e-7)
  expect_equal(with(res, n[populations=='a']), na, tolerance=1e-7)
  expect_equal(with(res, n[populations=='b']), nb, tolerance=1e-7)
  expect_equal(with(res, Hobs[populations=='a']), Hobsa, tolerance=1e-7)
  expect_equal(with(res, Hobs[populations=='b']), Hobsb, tolerance=1e-7)
  expect_equal(with(res, Hexp[populations=='a']), Hexpa, tolerance=1e-7)
  expect_equal(with(res, Hexp[populations=='b']), Hexpb, tolerance=1e-7)
})

test_that('What happens with numeric inputs?', {
  true <- Siccuracy:::make.true(9,4)
  true[] <- true[] + round(runif(prod(dim(true))), 2)
  true[true > 2] <- 2
  truefn <- tempfile('true', fileext = '.txt')
  write.snps(true, truefn)
  
  res <- heterozygosity(truefn)  
  
})