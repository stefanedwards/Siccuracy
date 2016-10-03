test_that('convert_plinkA correctly reads and writes files', {
  context('No calling to PLINK is performed.')
  
  raw <- Siccuracy:::make.true(8, 14)
  raw[,1] <- letters[1:nrow(raw)]
  raw[,2] <- raw[,1]
  raw[,3:4] <- 0
  raw[,5] <- sample.int(2, size=nrow(raw), replace=TRUE)
  raw[,6] <- round(rnorm(nrow(raw)), 5)
  raw[5,9] <- NA
  rawfn <- tempfile()
  write.snps(raw, rawfn, row.names = FALSE, na='NA')  

  newIDs <- 1:nrow(raw) + 30

  newfn <- tempfile()
  
  .Fortran('convertplinka', rawfn=as.character(rawfn), outputfn=as.character(newfn), newID=as.integer(newIDs),
                  ncol=as.integer(14-6), nrow=as.integer(8), 
                  stat=integer(1), PACKAGE='Siccuracy', NAOK=TRUE)
  
  convert_plinkA(rawfn, newfn, newID = newIDs)
})