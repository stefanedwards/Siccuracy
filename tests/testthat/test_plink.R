library(testthat)
library(Siccuracy)

context('Converting plink binary format to AlphaImpute format')

# Requires test files in ../..(/inst)/extdata/testdata, which can be constructed using ../../tools/make_plink_test_files.R

.datadir <- function(...) {
  # Find basedirectory
  BASEDIR <- getwd()
  while (basename(BASEDIR) != 'Siccuracy') BASEDIR <- dirname(BASEDIR)
  if (file.exists(file.path(BASEDIR, 'inst'))) BASEDIR <- file.path(BASEDIR, 'inst')
  file.path(BASEDIR, 'extdata', 'testdata', c(...))
}


expect_equal(file.exists(.datadir('')), TRUE, info=paste0('Run tools/make_plink_test_files.R to create `', .datadir(''), '`.'))
             
# Basic plink conversion works (`simple1`) ####
# Prepared in advance from make_plink_test_files.R:
# * binary files
# * --recode A file ('raw') file, although this file counts major alleles
# * 
test_that('Basic plink conversion works (`simple1`)', {
  fn1 <- tempfile()
  res <- convert_plinkA(.datadir('simple1.raw'), fn1)
  snp1 <- read.snps(fn1)
  
  fn2 <- tempfile()
  res <- convert_plink(.datadir('simple1'), fn2, countminor=FALSE)
  snp2 <- read.snps(fn2)
  expect_equal(snp1, snp2)
})

test_that('Basic plink routine correctly takes newID specifiers', {
  newID <- 100
  fn2 <- tempfile()
  res <- convert_plink(.datadir('simple1'), fn2, newID=newID)
  expect_equal(res$newID, newID + 1:6)
  expect_equal(res$newID, get_firstcolumn(fn2))
  
  newID <- c(55, 99, 10)
  res <- convert_plink(.datadir('simple1'), fn2, newID=newID)
  expect_equal(res$newID, newID)
  expect_equal(res$newID, get_firstcolumn(fn2))
  
  # Shuffled newID
  newID <- get_firstcolumn(.datadir('simple1.fam'), class=rep('character', 2), col.names=c('famID', 'sampID'))
  newID$newID <- 1:nrow(newID)
  newID <- newID[sample.int(nrow(newID)),]
  res <- convert_plink(.datadir('simple1'), fn2, newID=newID)
  newID <- newID[order(newID$newID),]
  expect_equal(res, newID)
  expect_equal(res$newID, get_firstcolumn(fn2))
  
  # To many newID's
  newID <- get_firstcolumn(.datadir('simple1.fam'), class=rep('character', 2), col.names=c('famID', 'sampID'))
  newID <- rbind(newID, data.frame(famID=c('a','a','b','b','b','d','d'), sampID=c('p','q','a','b','c','d','e')))
  newID <- newID[sample.int(nrow(newID)),]  
  newID$newID <- sample.int(nrow(newID))
  res <- convert_plink(.datadir('simple1'), fn2, newID=newID)
  snp1 <- read.snps(fn2)  
})

test_that('Duplicate IDs lead to issues', {
  fn2 <- tempfile()
  newID <- get_firstcolumn(.datadir('simple1.fam'), class=rep('character', 2), col.names=c('famID', 'sampID'))
  newID$order <- 1:nrow(newID)
  newID <- rbind(newID, data.frame(famID=c('a','a','b','b','b','d','d'), sampID=c('x','d','o','p','x','a','b'), order=NA))
  newID <- newID[sample.int(nrow(newID)),]  
  newID$newID <- sample.int(nrow(newID))
  res <- convert_plink(.datadir('simple1'), fn2, newID=newID)
  snp1 <- read.snps(fn2)  
  expect_true(any(duplicated(res[,1:2])))
  expect_equal(nrow(res), 1+length(get_firstcolumn(.datadir('simple1.fam'), class='character')))
  
  ni <- newID[!is.na(newID$order), ]
  ni <- ni[order(ni$order),]
  cols <- c('famID','sampID','newID')
  expect_equivalent(ni[,cols], res[!is.na(res$order),cols])
  
  expect_equal(res$newID, get_firstcolumn(fn2), info='Same length, i.e. including duplicate label.')
  
})

test_that('Count minor/major allele works', {
  fn1 <- tempfile()
  res <- convert_plinkA(.datadir('simple1.raw'), fn1)
  snp1 <- read.snps(fn1)
  
  fn2 <- tempfile()
  res <- convert_plink(.datadir('simple1'), fn2, countminor=FALSE)
  expect_equal(snp1, read.snps(fn2))
  snp.major <- read.snps(fn2, na=9)
  
  
  res <- convert_plink(.datadir('simple1'), fn2, countminor=TRUE)
  snp.minor <- read.snps(fn2, na = 9)
  expect_equal(snp.major, 2-snp.minor)
})


test_that('Simple converter accepts minimal allele frequencies', {
  fn1 <- tempfile()
  res <- convert_plinkA(.datadir('simple2_30.raw'), fn1)  # maf 0.30
  snp1 <- read.snps(fn1)  
  
  fn2 <- tempfile()
  res <- convert_plink(.datadir('simple2'), fn2, countminor=FALSE, maf=0.30)
  snp2 <- read.snps(fn2)
  expect_equal(snp1, snp2)
})