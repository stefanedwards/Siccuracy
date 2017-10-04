library(testthat)
library(Siccuracy)

test_that('convert_plinkA correctly reads and writes files', {
  context('No calling to PLINK is performed.')
  
  raw <- Siccuracy:::make.true(8, 18)
  raw[,1] <- letters[1:nrow(raw)]
  raw[,2] <- raw[,1]
  raw[,3:4] <- 0
  raw[,5] <- sample.int(2, size=nrow(raw), replace=TRUE)
  raw[,6] <- round(rnorm(nrow(raw)), 5)
  raw[5,9] <- NA
  raw[6,10] <- NA
  rawfn <- tempfile()
  write.snps(raw, rawfn, row.names = FALSE, na='NA')  

  newIDs <- 1:nrow(raw) + 30
  true <- apply(raw[,-c(1:6)], 2, as.integer)
  rownames(true) <- newIDs

  newfn <- tempfile()
  
  res <- convert_plinkA(rawfn, newfn, newID = 31)
  s <- read.snps(newfn, na=9)
  expect_equal(s, true)
})

test_that('convert_plinkA correctly uses a scalar as newID', {
  raw <- Siccuracy:::make.true(8, 18)
  raw[,1] <- letters[1:nrow(raw)]
  raw[,2] <- raw[,1]
  raw[,3:4] <- 0
  raw[,5] <- sample.int(2, size=nrow(raw), replace=TRUE)
  raw[,6] <- round(rnorm(nrow(raw)), 5)
  raw[5,9] <- NA
  raw[6,10] <- NA
  rawfn <- tempfile()
  write.snps(raw, rawfn, row.names = FALSE, na='NA') 
  
  true <- apply(raw[,-c(1:6)], 2, as.integer)
  newIDoffset <- 1000
  rownames(true) <- newIDoffset + 1:nrow(true) - 1
  
  newfn <- tempfile()
  
  res <- convert_plinkA(rawfn, newfn, newID = newIDoffset)
  s <- read.snps(newfn, na=9)
  expect_equal(s, true)
  
  context('Limiting number of rows...')
  res <- convert_plinkA(rawfn, newfn, newID = newIDoffset, nlines = 5)
  s <- read.snps(newfn, na=9)
  expect_equal(s, true[1:5,])
})
