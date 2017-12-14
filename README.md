# Siccuracy #

Stefan's imputation accuracy package. And with helper functions for working with AlphaSuite's SNP files.
Implemented in Fortran for speed.

## What is this repository for? ###

This package was developed for working with SNP files in the format used in [AlphaSuite](https://bitbucket.org/tutorials/markdowndemo). The format of these consists of first column with sample ID (Individual or Animal ID), and subsequent column counting minor alleles at each loci: 

```
1001      0    0    1    2    1    2
1002   0.05 0.50 1.50 1.80 0.95    9
```

In this example, the file consists of two individuals (1001 and 1002); first individual is has genotypes (`0`,`0`,`1`,`2`,`1`,`2`) for six loci. The second individual has been imputed, and the genotypes are given as real values of best guess. The last genotype was not sufficiently imputed and therefore set as missing with value ``9``.

Requirements for format: Space separated; number of digits is unimportant.

## To use ##

In R, load package with `library(Siccuracy)`. There exists a number of functions:

* Calculate imputation accuracy (i.e. correlation between two genotype matrices): `imputation_accuracy`
* Calculate hetereozygosity: `heterozygosity`
* Helper functions: `get_ncol`, `get_firstcolumn`, `get_nlines`
* Combining files: `rowconcatenate` and `mergeChips`
* Reading and writing files: `write.snps`, `read.snps`
* Convert phased files to genotype files: `phasotogeno`

**Most functions operate directly with files, not matrices loaded into memory!**


## How do I get set up? ###

Visiti https://github.com/stefanedwards/Siccuracy/releases to download 
the package for installation in R.
Zip-files are available for Windows users.

Alternatively, 
use `devtools` to install directly from github:

```
#!r

library(devtools)
install_github(repo='stefanedwards/Siccuracy', ref = "master")
```

## I found a bug ##

Lucky you.

Please report it using the issues tool.

## Can' some of these function be done in bash?

Presumable.

The `rowconcatenate` can be done by using `paste` in bash -- except it requires that you first strip off
first column of all files except the first.
