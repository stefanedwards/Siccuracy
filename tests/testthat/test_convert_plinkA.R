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
  
  res <- convert_plinkA(rawfn, newfn, newID = newIDs)
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
  rownames(true) <- newIDoffset + 1:nrow(true)
  
  newfn <- tempfile()
  
  res <- convert_plinkA(rawfn, newfn, newID = newIDoffset)
  s <- read.snps(newfn, na=9)
  expect_equal(s, true)
  
  context('Limiting number of rows...')
  res <- convert_plinkA(rawfn, newfn, newID = newIDoffset, nlines = 5)
  s <- read.snps(newfn, na=9)
  expect_equal(s, true[1:5,])
})

test_that('convert_plinkA correctly uses custom data.frame for newIDs', {
  raw <- Siccuracy:::make.true(8, 18)
  raw[,1] <- rep(LETTERS[1:3], each=3)[1:nrow(raw)]
  raw[,2] <- rep(letters[1:4], times=3)[1:nrow(raw)]
  raw[,3:4] <- 0
  raw[,5] <- sample.int(2, size=nrow(raw), replace=TRUE)
  raw[,6] <- round(rnorm(nrow(raw)), 5)
  raw[5,9] <- NA
  raw[6,10] <- NA
  raw[,18] <- sample.int(nrow(raw)*5, nrow(raw), replace=FALSE)
  raw <- raw[sample.int(nrow(raw)),]
  rawfn <- tempfile()
  write.snps(raw[,-18], rawfn, row.names = FALSE, na='NA') 
  
  true <- apply(raw[,-c(1:6,18)], 2, as.integer)
  rownames(true) <- as.character(raw[,18])
  
  newIDs <- data.frame(famID=raw[,1], sampID=raw[,2], newID=raw[,18])
  newIDs <- newIDs[sample.int(nrow(newIDs)),]
  
  newfn <- tempfile()
  res <- convert_plinkA(rawfn, newfn, newID = newIDs)
  s <- read.snps(newfn, na=9)
  expect_equal(s, true)
  
  context('Limiting number of rows...')
  res <- convert_plinkA(rawfn, newfn, newID = newIDs, nlines = 5)
  s <- read.snps(newfn, na=9)
  expect_equal(s, true[1:5,])
})