library(testthat)
library(Siccuracy)
library(vcfR)

# test files created with tools/make_vcf_test_files.R

data(vcfR_example) # loads object vcf

context('vcf files/objects')

# write.snps ------
test_that('write.snps works on vcf objects', {
  tfn <- tempfile()
  newIDs <- write.snps(vcf, file = tfn, na='9')
  raw_vcf <- readLines(tfn)  
  pfn <- system.file('extdata/testdata/vcfR.txt', package='Siccuracy')
  vcf_by_plinkA <- readLines(pfn)

  expect_equal(get_firstcolumn(tfn), get_firstcolumn(pfn))
  
  newID_by_plink <- read.table(system.file('extdata/testdata/vcfR.newids.txt', package='Siccuracy'), as.is=TRUE, header=TRUE)
  # Plink split sample names into family and sample by '_'; but these were only present a couple of places
  r <- with(newID_by_plink, famID != sampID)
  newID_by_plink$sampID[r] <- with(newID_by_plink[r,], paste(famID, sampID, sep='_'))
  newID_by_plink$famID <- NULL
  
  expect_equal(newID_by_plink, newIDs)
})

test_that('write.snps fails when given bad NewID', {
  newIDs <- data.frame(ID=colnames(vcf@gt), newID=letters[1:(dim(vcf)['gt_cols'])])
  expect_error(write.snps(vcf, file=tempfile(), newID = newIDs),
               regexp = 'all\\(c\\("sampID", "newID"\\) %in% names\\(newID\\)\\) is not TRUE')
  
  newIDs2 <- with(newIDs, data.frame(sampID=ID, newID=1:length(ID)))[1:5,]
  write.snps(vcf, file=tempfile(), newID = newIDs2)
})

test_that('write.snps handles some bad newIDs', {
  # write.snps fills in missing newIDs:
  newIDs <- data.frame(ID=colnames(vcf@gt), newID=letters[1:(dim(vcf)['gt_cols'])])

  newIDs2 <- with(newIDs, data.frame(sampID=ID, newID=1:length(ID)))[1:5,]
  v <- write.snps(vcf, file=tempfile(), newID = newIDs2)
  expect_type(v, 'list')
  expect_equivalent(nrow(v), as.integer(dim(vcf)['gt_cols'] - 1)) ## bloody names...
})

## write.snps returns all columns of newID -----
test_that('write.snps returns all columns of newID', {
  newIDs <- data.frame(sampID=colnames(vcf@gt), newID=1:(dim(vcf)['gt_cols']) +10, stringsAsFactors = FALSE)
  newIDs$nice <- letters[sample.int(length(letters), nrow(newIDs), replace=TRUE)]

  tfn <- tempfile()
  v <- write.snps(vcf, file=tfn, newID = newIDs)
  newIDs <- newIDs[-1,] # removes 'FORMAT' entry
  rownames(newIDs) <- NULL
  expect_equal(v, newIDs)
  snps <- read.snps(tfn)
  expect_equal(rownames(snps), as.character(newIDs$newID))
})

## write.snps fails when given vcfR without gt -----
test_that('write.snps fails when given vcfR without gt', {
  # object here is given by:
  # data(vcfR_test)  # return object named vcfR_test
  # dput(vcfR_test)
  # pasting back without gt slot
  
  v <- new("vcfR"
      , meta = c("##fileformat=VCFv4.3", "##fileDate=20090805", "##source=myImputationProgramV3.1", 
                 "##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta", 
                 "##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>", 
                 "##phasing=partial", "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">", 
                 "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">", 
                 "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">", 
                 "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">", 
                 "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">", 
                 "##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">", 
                 "##FILTER=<ID=q10,Description=\"Quality below 10\">", "##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">", 
                 "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", 
                 "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">", 
                 "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">", 
                 "##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">"
      )
      , fix = structure(c("20", "20", "20", "20", "20", "14370", "17330", "1110696", 
                          "1230237", "1234567", "rs6054257", NA, "rs6040355", NA, "microsat1", 
                          "G", "T", "A", "T", "GTC", "A", "A", "G,T", NA, "G,GTCT", "29", 
                          "3", "67", "47", "50", "PASS", "q10", "PASS", "PASS", "PASS", 
                          "NS=3;DP=14;AF=0.5;DB;H2", "NS=3;DP=11;AF=0.017", "NS=2;DP=10;AF=0.333,0.667;AA=T;DB", 
                          "NS=3;DP=13;AA=T", "NS=3;DP=9;AA=G"), .Dim = c(5L, 8L), .Dimnames = list(
                            NULL, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", 
                                    "INFO")))
      )
  expect_error(write.snps(v, file=tempfile()))
})


# Imputation accuracy -----------------

data(vcfR_example) # loads object vcf

## imputation_accuracy on VCF objects are entirely true (object: vcf) ------
test_that('imputation_accuracy on VCF objects are entirely true (object: vcf)', {
  m <- vcfR::maf(vcf)
  res <- imputation_accuracy(vcf, vcf)
  
  expect_equal(res$matcor, 1.0)
  
  snps.0 <- snps.1 <- rep(1, dim(vcf)['variants'])
  NAs <- m[,'Count'] <= 1 # subset of rows with no variance
  snps.0[res$snps$sds == 0 | is.na(res$snps$sd)] <- NA
  expect_equal(nrow(res$snps), length(snps.1))
  expect_equal(res$snps$cors, snps.0)
  expect_equal(res$snps$correct.pct, snps.1)
  
  
  samples.1 <- rep(1, dim(vcf)['gt_cols'] - 1)
  expect_equal(nrow(res$animals), length(samples.1))
  expect_equal(res$animals$cors, samples.1)
  expect_equal(res$animals$correct.pct, samples.1)
  
  expect_equal(res$snps$both.na, as.integer(m[,'NA']))
})