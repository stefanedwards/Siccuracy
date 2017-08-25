# Uses vcf data from library vcfR package;
# writes it to a file,
# uses > plink < to recode to plinkA format,
# which is finally deposited as a AlphaImpute format.

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


require(vcfR)
require(Siccuracy)

data(vcfR_example)
# vcf

vcf.fn <- 'vcfR.vcf.gz'
vcfR::write.vcf(vcf, file = vcf.fn)
plink(bfile=NULL, '--vcf', vcf.fn,  '--recode A', '--allow-extra-chr', out='vcfR')

newIds <- convert_plinkA('vcfR.raw', outfn = .datadir('vcfR.txt'))
write.table(newIDs, .datadir('vcfR.newids.txt'), row.names=FALSE, col.names = TRUE, quote=FALSE)


