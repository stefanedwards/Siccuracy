# Make PLINK test files, deposits in data-folder.
# Requires PLINK is downloaded in tools-directory.

if (basename(getwd()) != 'tools') {
  warning('Setting working directory to `tools` subdirectory!')
  try(setwd('tools'))
  try(setwd('../tools'))
  stopifnot(basename(getwd()) == 'tools')
  warning('Working directory is now ', getwd())
}

plink <- function(bfile, ..., out=NA, stdout='') {
  args <- c(...)
  if (!is.null(bfile)) args <- c(args, '-bfile', bfile)
  if (!is.na(out)) args <- c(args, '--out', out)
  if (.Platform$OS.type == 'unix') {
    args=c(args,'--threads', Sys.getenv('NSLOTS', '1'))
    bin='/exports/igmm/software/pkg/el7/apps/plink/1.90b1g/bin/plink'  # specific for running on Eddie3.
    if (!file.exists(bin)) bin <- 'plink'  # use local located
    stopifnot(file.exists(bin))
    s <- system2(bin, args=args, stdout=stdout,stderr=stdout)
  } else {
    bin='plink.exe' # Got as https://www.cog-genomics.org/static/bin/plink160607/plink_win64.zip
    stopifnot(file.exists(bin))
    s <- system2(bin, args=args, stdout=stdout,stderr=stdout)
  }
  invisible(s)
}

.datadir <-  function(...) file.path('..', 'inst', 'extdata', 'testdata', c(...))
null <- dir.create(.datadir(''), recursive = TRUE)

# Make simple example from website ####
# see examples https://www.cog-genomics.org/plink2/formats#bed and http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
fam <- data.frame(
  famID='a',
  sampleID=letters[1:6],
  pedp=0,
  pedm=0,
  sex=1,
  y=1,
  snp1a=c('G','A','0','A','A','A'),
  snp1b=c('G','A','0','A','A','A'),
  snp2a=c(2,0,1,2,2,2),
  snp2b=c(2,0,2,2,2,2),
  snp3a=c('C','A','A','0','0','A'),
  snp3b=c('C','C','C','0','0','A'),
  stringsAsFactors = FALSE
)
# see examples https://www.cog-genomics.org/plink2/formats#bed and http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
write.table(fam, 'test.ped', row.names=FALSE, col.names=FALSE, quote=FALSE)
map <- data.frame(chr=1, rs=c('snp1','snp2','snp3'), ps=0, id=1:3)
write.table(map, 'test.map', row.names=FALSE, col.names=FALSE, quote=FALSE)
plink(bfile=NULL, '--file test', '--make-bed',out=.datadir('simple1'))
plink(bfile=.datadir('simple1'), '--recode A', out=.datadir('simple1'))


# A bigger example with multiple chromosomes ####
# Limit chromosomes to 23 as plink should otherwise be called with --nchrom.


set.seed(140111)
n <- 10
m <- 25
fam <- data.frame(
  famID=sample(letters[1:2], n, replace=TRUE, prob=c(0.7, 0.3)),
  sampleID=letters[sample.int(15, n, replace=FALSE)],
  pedp=0,
  pedm=0,
  sex=1,
  y=2
)
# Populate random chromosomes
random.allele <- function(n) {
  al <- sample(c('A','C','G','T'), 2, replace=FALSE)
  p <- runif(1, 0.1, 0.5)
  ge <- matrix(sample(al, size=2*n, replace=TRUE, prob=c(p,1-p)), ncol=2)
  list(ge[,1], ge[,2])
}
snps <- t(do.call(cbind, replicate(m, random.allele(n))))
miss <- ceiling(sample.int(prod(dim(snps))-2, n*m*0.1)/2)*2
snps[c(miss,miss-1)] <- NA
snps <- t(snps)
majmin <- do.call(rbind, lapply(seq(1,m*2,by=2), function(i) {names(sort(table(na.omit(snps[,(i):(i+1)])), decreasing = TRUE))}))
chrs <- sample.int(10, 2)
chrs <- c(chrs, m-sum(chrs))
map <- cbind(data.frame(chr=rep.int(1:3, chrs), rs=paste0('snp', 1:m), cm=0, pos=1:m), majmin)

write.table(cbind(fam, snps), 'test.ped', row.names=FALSE, col.names=FALSE, quote=FALSE, na='0')
write.table(map, 'test.map', row.names=FALSE, col.names=FALSE, quote=FALSE)

plink(bfile=NULL, '--file test', '--make-bed', '--maf 0.0001' ,out=.datadir('simple2'))
plink(bfile=.datadir('simple2'), '--recode A', out=.datadir('simple2'))
plink(bfile=.datadir('simple2'), '--recode A', '--maf 0.30', out=.datadir('simple2_30'))
