# running time of get_nlines vs. get_firstcol:
library(Siccuracy)

fn <- tempfile()
print('Creating file:')
system.time(write.snps(Siccuracy:::make.true(6000, 5000), fn))

print('get_firstcolumn:')
system.time(length(get_firstcolumn(fn)))
print('get_nlines:')
system.time(get_nlines(fn))

t <- Siccuracy:::make.true(20, 30)
rownames(t) <- replicate(nrow(t), {paste0(letters[sample.int(26, runif(1,3,8))], collapse='')}, simplify=TRUE)
fn <- tempfile()
write.snps(t, fn)
print('get_firstcolumn:')
system.time(n1 <- length(get_firstcolumn(fn, 'character')))
print('get_nlines:')
system.time(n2 <- get_nlines(fn))
