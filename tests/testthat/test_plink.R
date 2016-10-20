library(testthat)
library(Siccuracy)


system2('plink.exe', c('-bfile plink', '--recode A'))
convert_plinkA('plink.raw', 'plink.txt')


# subroutine readplink(ped, bim, bed, fnout, ncol, nlines, na)
.Fortran('readplink', ped=as.character('plink.bed'), bim=as.character('plink.bim'), bed=as.character('plink.bed'), fnout=as.character('new.txt'),
         ncol=as.integer(3), nlines=as.integer(6), na=as.integer(9), PACKAGE='Siccuracy')
