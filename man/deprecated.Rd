% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{Siccuracy-deprecated}
\alias{Siccuracy-deprecated}
\alias{cbind_SNPs}
\alias{rowconcatenate}
\alias{imputation_accuracy3}
\alias{imputation_accuracy1}
\alias{mergeChips}
\alias{phasotogeno}
\alias{phasotogeno_int}
\alias{rbind_SNPs}
\title{Deprecated functions}
\usage{
cbind_SNPs(fns, fnout, nlines, ncols, skiplines, excludeids, int)

rowconcatenate(fns, fnout, nlines = NULL, ncols = NULL, skiplines = 0,
  excludeids = integer(0))

imputation_accuracy3(truefn, imputefn, nSNPs = NULL, nAnimals = NULL,
  NAval = 9, standardized = TRUE)

imputation_accuracy1(truefn, imputefn, nSNPs = NULL, nAnimals = NULL,
  NAval = 9, standardized = TRUE)

mergeChips(hdid, ldid, hdpos, ldpos, hdfn, ldfn, fnout, outpos = NULL,
  missing = 9)

phasotogeno(phasefn, genofn, ncol = NULL, nrow = NULL)

phasotogeno_int(phasefn, genofn, ncol = NULL, nrow = NULL)

rbind_SNPs(hdid, ldid, hdpos, ldpos, hdfn, ldfn, fnout, outcol = NULL,
  na = 9, format = NULL, int = TRUE)
}
\arguments{
\item{fns}{Character vector of filenames to concatenate. Max 255 characters per filename due to R-Fortran.}

\item{fnout}{Filename of resulting file}

\item{nlines}{Number of lines to read from each file. When \code{NULL} (default), use all lines.}

\item{ncols}{Number of columns to read from each file. When \code{NULL} (default), use all columns.}

\item{skiplines}{Integer, number of lines to skip before outputting. \code{nlines} counts towards total number of lines, i.e. both skipped and outputted.}

\item{excludeids}{Integer vector of first column id's to \emph{exclude} from the output. \code{nlines} also counts towards excluded lines.}

\item{int}{Logical (default \code{TRUE}), read and write integers.}

\item{truefn, imputefn}{Deprected, was \emph{filename} to files with true and imputed genotype matrix.}

\item{nSNPs}{Deprecated, use \code{ncol}.}

\item{nAnimals}{Deprecated, use \code{nlines}.}

\item{NAval}{Deprecated, use \code{na}.}

\item{standardized}{Logical, whether to center and scale genotypes by dataset in \code{true}-matrix.
Currently by subtracting column mean and dividing by column standard deviation.}

\item{hdid}{IDs of HD genotyped individuals. Corresponds to (subset of) first column of \code{hdfn}.}

\item{ldid}{IDs of LD genotyped individuals. Corresponds to (subset of) first column of \code{ldfn}.}

\item{hdpos}{Integer vector of where columns in \code{hdfn} are positioned in output file. Coerced to same length as columns in \code{hdfn}.}

\item{ldpos}{Integer vector of where columns in \code{ldfn} are positioned in output file. Coerced to same length as columns in \code{ldfn}.}

\item{hdfn}{Filename of HD genotypes.}

\item{ldfn}{Filename of LD genotypes.}

\item{outpos}{Integer vector of collective SNP positions. Default to sorted, union of \code{hdpos} and \code{ldpos}. Make it anything else and you get?}

\item{missing}{Missing value.}

\item{phasefn}{Filename of input file, every two rows are for same animal.}

\item{genofn}{Filename of intended output.}

\item{ncol}{Number of columns to read, if \code{NULL} (default), number is estimated from file.}

\item{nrow}{Number of rows to maximally read from \code{phasefn}. If \code{NULL}, no limit is used.}

\item{outcol}{Integer, number of columns in output file. When \code{NULL} (default), uses max value of \code{hdpos} and \code{ldpos}.}

\item{na}{Value of missing genotypes.}

\item{format}{Character, Fortran edit descriptors for output. See \link{parseformat}.}
}
\value{
Number of rows written.
}
\description{
Mostly renamed function names.

\code{cbind_SNPs} has been renamed to \code{\link{cbind_snp_files}}.

\code{rowconcatenate} has been renamed to \code{\link{cbind_snp_files}}.

\code{imputation_accuracy1} and \code{imputation_accuracy3} has been replaced by \code{\link{imputation_accuracy}}.
The difference between the two former functions is now covered by the \code{adaptive}-argument of the latter.

\code{mergeChips} has been renamed to \code{\link{rbind_snp_files}}.

\code{phasotogeno} and \code{phasotogeno_int} have been replaced by 
\code{\link{convert_phases}}. 
The difference between the two former functions is covered by the \code{int} argument in the new function.

\code{rbind_SNPs} has been renamed to \code{\link{rbind_snp_files}}.
}
