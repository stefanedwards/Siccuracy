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
       
test_that('Lowmem plink conversion does not blow up',{
  outfiles <- c(paste0('simple2_chr',1:3,'.txt'), paste0('bed',1:3))
  res <- convert_plink(.datadir('simple2'), outfn='eh', method='lowmem', fragment='chr', fragmentfns=outfiles)
  res
})

      
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
  expect_equal(res$newID$newID, newID + 1:6)
  expect_equal(res$newID$newID, get_firstcolumn(fn2))
  
  newID <- c(55, 99, 10)
  res <- convert_plink(.datadir('simple1'), fn2, newID=newID)
  expect_equal(res$newID$newID, newID)
  expect_equal(res$newID$newID, get_firstcolumn(fn2))
  
  # Shuffled newID
  newID <- get_firstcolumn(.datadir('simple1.fam'), class=rep('character', 2), col.names=c('famID', 'sampID'))
  newID$newID <- 1:nrow(newID)
  newID <- newID[sample.int(nrow(newID)),]
  res <- convert_plink(.datadir('simple1'), fn2, newID=newID)
  newID <- newID[order(newID$newID),]
  expect_equal(res$newID, newID)
  expect_equal(res$newID$newID, get_firstcolumn(fn2))
  
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
  expect_true(any(duplicated(res$newID[,1:2])))
  expect_equal(nrow(res$newID), 1+length(get_firstcolumn(.datadir('simple1.fam'), class='character')))
  
  ni <- newID[!is.na(newID$order), ]
  ni <- ni[order(ni$order),]
  cols <- c('famID','sampID','newID')
  expect_equivalent(ni[,cols], res$newID[!is.na(res$newID$order),cols])
  
  expect_equal(res$newID$newID, get_firstcolumn(fn2), info='Same length, i.e. including duplicate label.')
  
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

test_that('Simple maf filtering work', {
  fn1 <- tempfile()
  res <- convert_plinkA(.datadir('simple2.raw'), fn1)  # maf 0.30
  snp1 <- read.snps(fn1)  
  
  m <- ncol(snp1)
  extract <- seq.int(1, m, by=3)
  res <- convert_plink(.datadir('simple2'), fn1, extract=extract, countminor=FALSE)
  snp2 <- read.snps(fn1)  
  expect_equal(snp2, snp1[,extract])
})

test_that('Filtering works', {
  inp <- .datadir('simple2')
  fn <- tempfile()
  fam <- get_firstcolumn(Siccuracy:::sub_ext(inp, 'fam'), class=c('character','character'), col.names=c('famID','sampID'))
  rs <- get_firstcolumn(Siccuracy:::sub_ext(inp, 'bim'), class=c('NULL','character'))
  raw <- read.table(Siccuracy:::sub_ext(inp, 'raw'), header=TRUE, as.is = TRUE)
  names(raw)[1:2] <- names(fam)
  #simple2.raw is recoded by plink.
  
  res <- convert_plink(inp, outfn='', method='dryrun')
  expect_equal(res$ncol, length(rs))
  expect_equal(res$nlines, nrow(fam))
  expect_equal(sum(res$keep), nrow(fam))
  expect_equal(sum(res$extract), length(rs))
  
  res <- convert_plink(inp, outfn='', method='dryrun', remove=list(famID='b'))
  expect_equal(as.logical(res$keep), fam$famID != 'b')
  
  remove=list(sampID=letters[4:6])
  keep='a'
  rows <- (fam$famID %in% keep | fam$sampID %in% keep) & !fam$sampID %in% remove$sampID
  res <- convert_plink(inp, outfn=fn, method='dryrun', remove=remove, keep=keep)
  expect_equal(as.logical(res$keep), rows )
  res <- convert_plink(inp, outfn=fn, method='simple', remove=remove, keep=keep, countminor=FALSE)
  expect_equal(get_firstcolumn(fn), res$newID$newID[as.logical(res$keep)])
  expect_equal(res$newID[as.logical(res$keep),1:2], fam[rows,])
  expect_equivalent(read.snps(fn, na=9), as.matrix(raw[(raw$famID %in% keep | raw$sampID %in% keep) & !raw$sampID %in% remove$sampID, -c(1:6)]))
  
  keep <- logical(nrow(fam))
  keep[sample.int(nrow(fam), 5)] <- TRUE
  remove <- fam$sampID[sample.int(nrow(fam), 4)]
  res <- convert_plink(inp, outfn=fn, method='dryrun', remove=remove, keep=keep)
  expect_equal(as.logical(res$keep), (keep) & !(fam$famID %in% remove | fam$sampID %in% remove))
  res <- convert_plink(inp, outfn=fn, method='simple', remove=remove, keep=keep, countminor=FALSE)
  expect_equal(get_firstcolumn(fn), res$newID$newID[as.logical(res$keep)])
  expect_equal(res$newID[as.logical(res$keep),1:2], fam[(keep) & !(fam$famID %in% remove | fam$sampID %in% remove),])
  expect_equivalent(read.snps(fn, na=9), as.matrix(raw[(keep) & !(raw$famID %in% remove | raw$sampID %in% remove), -c(1:6)]))
  
  remove <- keep <- logical(nrow(fam))
  keep[sample.int(nrow(fam), 5)] <- TRUE
  remove[sample.int(nrow(fam), 5)] <- TRUE
  res <- convert_plink(inp, outfn=fn, method='dryrun', remove=remove, keep=keep)
  expect_equal(sum(res$keep), sum(keep & !remove))
  res <- convert_plink(inp, outfn=fn, method='simple', remove=remove, keep=keep, countminor=FALSE)
  expect_equal(get_firstcolumn(fn), res$newID$newID[as.logical(res$keep)])
  expect_equal(res$newID[as.logical(res$keep),1:2], fam[keep & !remove,])
  expect_equivalent(read.snps(fn, na=9), as.matrix(raw[keep & !remove, -c(1:6)]))
  
  
  context('Loci')
  extract <- sample(rs, max(length(rs), 20))
  exclude <- sample(rs, 5)
  res <- convert_plink(inp, outfn=fn, method='simple', countminor=FALSE, extract=extract, exclude=exclude)
  expect_equal(sum(res$extract), length(setdiff(extract, exclude)))  
  expect_equivalent(read.snps(fn, na=9), as.matrix(raw[,which(rs %in% setdiff(extract, exclude))+6]))
  
  
  extract <- logical(length(rs))
  extract[] <-  TRUE
  extract[sample.int(length(rs), 8)] <- FALSE
  res <- convert_plink(inp, outfn=fn, method='simple', countminor=FALSE, extract=extract, exclude=exclude)
  expect_equal(as.logical(res$extract), extract & !rs %in% exclude)  
  expect_equivalent(read.snps(fn, na=9), as.matrix(raw[,which( extract & !(rs %in% exclude) )+6]))
  
})