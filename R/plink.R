#' Converts PLINK binary format to flat format.
#' 
#' @references 
#' \itemize{
#'  \item PLINK v. 1.07 BED file format: \url{http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml}
#'  \item Shaun Purvell and Christopher Chang. \emph{PLINK v. 1.90} \url{https://www.cog-genomics.org/plink2}
#'  \item Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) \href{http://www.gigasciencejournal.com/content/4/1/7}{Second-generation PLINK: rising to the challenge of larger and richer datasets.} \emph{GigaScience}, 4.
#' }
#

make.bed <- function() {
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
  system2("C:/Users/shojedw/Documents/Projects/Siccuracy/Siccuracy/tools/plink.exe", args=c('--file test','--make-bed'))
}