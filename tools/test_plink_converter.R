library(Siccuracy)

outfn <- "C:/Users/shojedw/Documents/Projects/Siccuracy.extra/stuff/cache_canine_imputation/cornell_canine.txt"
inp <- "C:/Users/shojedw/Documents/Projects/Siccuracy.extra/stuff/cache_canine_imputation/cornell_canine"
chrfn <- function(i) sprintf('C:/Users/shojedw/Documents/Projects/Siccuracy.extra/stuff/cache_canine_imputation/cornell_canine_chr%02i.txt', i)

res <- convert_plink(inp, outfn, method='lowmem', fragmentfns = chrfn, remerge=FALSE)
print(res$fragmentfns)