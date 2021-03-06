#' Stefan's imputation accuracy package
#' 
#' Calculation of imputation accuracies and similar allele stats.
#' 
#' This package was developed for working with SNP files in the format used in \href{http://www.alphagenes.roslin.ed.ac.uk/}{AlphaSuite}. 
#' The format of these consists of first column with \emph{integer} sample ID (Individual or Animal ID),
#' and subsequent column counting minor alleles at each loci:
#' 
#' \tabular{lcccccc}{
#'  1001  \tab    0   \tab   0    \tab 1    \tab  2   \tab  1   \tab   2 \cr
#'  1002    \tab 0.05   \tab 0.50   \tab 1.50   \tab 1.80   \tab 0.95     \tab 9 \cr
#' }
#' 
#' In this example, the file consists of two individuals (\code{1001} and \code{1002}); 
#' first individual has genotypes (\code{0}, \code{0}, \code{1}, \code{2}, \code{1}, \code{2}) for six loci. 
#' The second individual has been imputed, and the genotypes are given as real values of best guess. 
#' The last genotype was not sufficiently imputed and therefore set as missing with value \code{9}.
#'
#' Requirements for format: Space or tab separated; number of digits is unimportant.
#' 
#' @section Overview:
#' 
#' \itemize{
#'   \item \strong{Converters} \code{\link{convert_phases}}, \code{\link{convert_plink}}, \code{\link{convert_plinkA}}.
#'   \item \strong{File info} \code{\link{get_nlines}}, \code{\link{get_firstcolumn}}, \code{\link{get_ncols}}.
#'   \item \strong{Imputation accuracies} \code{\link{imputation_accuracy}} correlations between matrices, and tabulating correctly imputed genotypes.
#'   \item \strong{Heterozygosity} \code{\link{heterozygosity}} counting alleles.
#'   \item \strong{Mutators} \code{\link{cbind_snp_files}}, \code{\link{rbind_snp_files}} for combining multiple files, \code{\link{mask_SNPs}} for masking SNPs for simulating different density SNP chips.
#'   \item \strong{Reading and writing} \code{\link{write.snps}}, \code{\link{read.snps}},
#'                                      \code{\link{read.haps}}, \code{\link{read.oxford}},
#'                                      \code{\link[=VCF format]{read.vcf}}.
#' }
#' 
#' @section IDs:
#' 
#' The AlphaImpute format \strong{must} have integer IDs. 
#' In this software package, we have choosen an arbitrary upper limit of 20 digits for IDs.
#' Exceed this limit and the ID will be lumped together with the following genotype.
#' 
#' Regarding \emph{repeated} IDs, it is best to avoid these.
#' 
#' @section Output format:
#' 
#' Functions that write new files may have a pair of arguments (\code{int} and \code{format}) 
#' that specifies whether the outputted format are integers (\code{int=TRUE}) or have decimals.
#' For more information on how to specify the format, see \link{parseformat}.
#' 
#' @section Other formats:
#' 
#' \strong{Oxford} format is covered with \code{\link{read.oxford}}
#' and a method dispatch for imputation accuracy exists for this too.
#' 
#' \strong{SHAPEIT}'s haps/sample format is covered with \code{\link{read.haps}} 
#' and a method dispatch for imputation accuracy exists for this too.
#' 
#' Currently, \strong{VCF} objects by package \code{vcfR} are supported for functions
#' \code{\link[=write.snps.vcfR]{write.snps}} and \code{\link[=imputation_accuracy.vcfR]{imputation_accuracy}}.
#' For more information, see \link{VCF format}.
#' 
#' 
#' @docType package
#' @name Siccuracy
#' @aliases AlphaImpute-format
#' @useDynLib Siccuracy
#' @author Stefan McKinnon Edwards <sme@@iysik.com>
NULL
