# Check against Andrew's results in secret directory, that totally doesn't make sense unless you have the files here.

setwd('C:\\Users\\shojedw\\Documents\\Projects\\Siccuracy\\1-snp_1000-ld_200-hd_0-multi_100')

# Run Andrew's script ---------------

if (!file.exists('andrews_combined_results.Rdata')) {
  source('analyzeResults.r')
  save(combined, true, genotypes, targetSnps, maf, stdev, file='andrews_combined_results.Rdata')
} else {
  load('andrews_combined_results.Rdata')
}


# Redo with Siccuracy ------
library(Siccuracy)
library(testthat)

test_that('Input and true genotypes corresponds', {
  my.true <- read.snps('trueGenotype.txt')
  ttrue <- true  # From Andrew's script
  dimnames(ttrue) <- list(NULL, NULL)
  
  dimnames(my.true) <- list(NULL, NULL)
  expect_equal(my.true, ttrue)
  
  my.genotypes <- read.snps('genotypes.txt')
  dimnames(my.genotypes) <- list(NULL, NULL)
  
  tgeno <- genotypes   # From Andrew's script
  dimnames(tgeno) <- list(NULL, NULL)
  expect_equal(tgeno, my.genotypes)
})

test_that('My imputation_accuracy gives same statistics as Andrew\'s maf and targetSNPs', {
  res <- imputation_accuracy('trueGenotype.txt', 'genotypes.txt', standardize=TRUE)
  expect_equivalent(res$snps$means, maf)
  
})

test_that('FastPhase gives similar results:', {
  # Copy-paste from Andrew's script
  tmp = readLines("postPhase_hapguess_switch.out")
  
  tmp = tmp[nchar(tmp)> 500]
  tmp = gsub("A", 0, tmp, fixed=TRUE)
  tmp = gsub("C", 1, tmp, fixed=TRUE)
    
  rawValues = as.matrix(read.table(textConnection(tmp)))
  even = (1:nrow(rawValues))%%2 == 0
  odd = (1:nrow(rawValues))%%2 == 1
    
  fastPhase = rawValues[even,] + rawValues[odd,]
  fastPhase.values = getMetrics(fastPhase, true, "fastPhase")
  # End of copy-paste
  
  values <- structure(fastPhase.values[,3], .Names=as.character(fastPhase.values[,2]))
    
  # Note that Andrew used the median to summarise correlations
  
  res <- imputation_accuracy(true, fastPhase, standardize=FALSE)
  expect_true(abs(median(res$animals$cors, na.rm=TRUE) - values['breedingValueCorrelation']) < 1e-9)
  expect_true(abs(median(res$snps$cors, na.rm=TRUE) - values['gwasCorrelation']) < 1e-9)  
})

test_that('eagleOutput has similar results',{
  # Copy-paste
  eagleRaw = read.table("eagleOutput.haps")
  reference = eagleRaw[,4]
  eagleRaw = as.matrix(eagleRaw[,-(1:5)])
  even = (1:ncol(eagleRaw))%%2 == 0
  odd = (1:ncol(eagleRaw))%%2 == 1
  
  eagle = eagleRaw[,even] + eagleRaw[,odd]
  eagle[reference == "C",] = 2- eagle[reference == "C",]
  eagle = t(eagle)
  
  eagle.values = getMetrics(eagle, true, "eagle")
  ###
  
  values <- structure(eagle.values[,3], .Names=as.character(eagle.values[,2]))
  
  ##
  
  haps <- read.haps('eagleOutput')
  expect_equal(reference == 'C', haps$map$A1 == 'C')
  
  h2 <- extract.snps(haps)
  h2[, haps$map$A1 == 'C'] <- 2 - h2[, haps$map$A1 == 'C'] ## switching coding
  
  dimnames(eagle) <- dimnames(h2)
  expect_equal(eagle, h2)
  
  res <- imputation_accuracy(true, h2, standardized=FALSE)
  expect_true(abs(median(res$animals$cors, na.rm=TRUE) - values['breedingValueCorrelation']) < 1e-9)
  expect_true(abs(median(res$snps$cors, na.rm=TRUE) - values['gwasCorrelation']) < 1e-9)  
  
})
