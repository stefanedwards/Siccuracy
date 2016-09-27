% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputation_accuracy.R
\name{imputation_accuracy}
\alias{imputation_accuracy}
\title{Imputation accuracy, aka. correlations}
\usage{
imputation_accuracy(truefn, imputefn, nSNPs = NULL, nAnimals = NULL,
  NAval = 9, standardized = TRUE, fast = FALSE)
}
\arguments{
\item{truefn}{Filename to true matrix. NB! Max. 255 characters!}

\item{imputefn}{Filename to imputed matrix.}

\item{nSNPs}{Integer, number of snps (i.e. number of columns minus 1).
If \code{NULL} (default), detected from \code{truefn}.}

\item{nAnimals}{Integer, number of rows.
If \code{NULL} (default), detected from \code{truefn}.}

\item{NAval}{Integer, value of missing genotype.}

\item{standardized}{Logical, whether to center and scale genotypes by dataset in \code{true}-matrix.
Currently by subtracting column mean and dividing by column variance.}

\item{fast}{Use method that does not check row IDs.}
}
\value{
List with following elements:
\describe{
  \item{\code{means}}{Column means of true matrix.}
  \item{\code{vars}}{Column variances of true matrix.}
  \item{\code{rowcors}}{Row-wise (animal-wise) correlation between true and imputed matrix.}
  \item{\code{matcor}}{Matrix-wise correlation between true and imputed matrix.}
  \item{\code{colcors}}{Column-wise (locus-wise) correlation between true and imputed matrix.}
}
}
\description{
Calculation of column-wise, row-wise, and matrix-wise correlations between
matrix in file \code{truefn} and \code{imputefn}.
Assumes first column in both files is an integer ID column and thus excluded.
Standardization (subtract mean, divide by variance) is done column-wise based
on means and variances of \code{truefn}.
Correlations are only performed on those rows that are found in \emph{both} files,
based on the first column (ID column).
}
\details{
This method stores the "true" matrix in memory with a low-precision real type,
and rows in the "imputed" matrix are read and matched by ID.
If there are no extra rows in either matrix and order of IDs is the same,
consider using \code{fast}-option that does not store either.
}
\seealso{
\code{\link{write.snps}} for writing SNPs to a file.
}
