% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heterozygosity.R
\name{heterozygosity}
\alias{heterozygosity}
\title{Calculates hetereozygosity stats for a single population}
\usage{
heterozygosity(fn, population = NULL, ncol = NULL, nlines = NULL)
}
\arguments{
\item{fn}{Filename of genotype matrix (0,1,2) with first column denoting ID.}

\item{population}{Vector of same length as rows in \code{fn}; defaults to \code{1}, coerced from factor to integer.}

\item{ncol}{Integer, number of SNP columns in \code{fn} When \code{NULL}, automagically detected with \code{get_ncols(fn)-1}.}

\item{nlines}{Integer, number of lines in \code{fn}. When \code{NULL}, automagically detected with \code{gen_nlines(fn)}.}
}
\value{
Data frame with columns
\describe{
  \item{\code{population}}{Population, as specified by argument \code{population}.}
  \item{\code{p}}{Vector of allele frequencies of alleles coded as '0'.}
  \item{\code{Hobs}}{Vector of observed heterozygosity (\eqn{H_{obs}}{Hobs}).}
  \item{\code{Hexp}}{Vector of expected heterozygosity (\eqn{H_{exp}}{Hexp}).}
  \item{\code{n}}{Vector of observed genotypes for each locus (ignoring NA-values). Alleles are twice this number.}
}
}
\description{
Calculate gene frequencies and observed and expected heterozygosites for locus.
These routines only work on a genotype file (e.g. as written by \code{\link{write.snps}}) with first column as ID column,
and subsequent columns are genotypes coded as \code{0}, \code{1}, or \code{2}.
If \code{population} is given, all calculations are performed on each population separately.
}
\details{
The gene frequencies, \eqn{p} and \eqn{q}, refers to the gene frequencies of genotypes \code{0} and \code{2}, respectively.
They are calculated as such, for each column:
\deqn{
  p = \frac{\# 0 \cdot 2 + \# 1 }{2 \cdot \# rows }
}{p = (2 * # homozygotes + # heterozygotes) / (2 * # rows)}
\deqn{
  q = 1 - p
}

Observed (\eqn{H_{obs}}{Hobs}) and expected (\eqn{H_{exp}}{Hexp}) heterozygosity are calculated as
\deqn{
  H_{obs} = \frac{\# 1}{\# rows}
  }{Hobs = #1/n}
\deqn{
  H_{exp} = 2 \cdot p \cdot q
  }{Hexp = 2*p*q}

\strong{Missing values:}
In the above, \code{NA} elements are ignored and thus do not count toward either genotype nor number of rows.
Genotypes with values smaller than 0 or greater than 2 are considered missing.

Finally, an inbreeding coeffiecent would be calculated as:
\deqn{
  F = \frac{H_{exp} - H_{obs}}{H_{exp}}
}{F = (Hexp - Hobs)/Hexp}
For which \eqn{F_{st}}{Fst} estimator to use, and how to combine across multiple
snps, see Bhatia et al. (2013).
}
\references{
Equations based on \url{http://www.uwyo.edu/dbmcd/molmark/practica/fst.html}.

Bhatia, Patterson, Sankararaman, and Price. Estimating and interpreting FST: The impact of rare variants.
\emph{Genome Research} (2013) 23: 1514-1521. \href{http://genome.cshlp.org/content/early/2013/07/16/gr.154831.113.full.pdf}{Preprint}
\href{http://dx.doi.org/10.1101/gr.154831.113}{doi: 10.1101/gr.154831.113}.
}
