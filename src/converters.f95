! Converts phase-file to genotype file.
!numfmt should default to '5.2'
subroutine convert_phase_num(phasefn, genofn, ncol, nrow, naval, numfmt, lennumfmt)

  implicit none

  !! Arguments
  character(255), intent(in) :: phasefn, genofn
  character(len=lennumfmt) :: numfmt
  integer, intent(in) :: ncol, naval, lennumfmt
  integer, intent(inout) :: nrow

  !! Local variables
  integer :: stat, i, animalid
  real, dimension(ncol) :: linea, lineb, summ
  character(100) :: nChar, fmt

  write(nChar,*) ncol
  fmt='(i20,'//trim(adjustl(nChar))//'F'//trim(numfmt)//')'

  open(97, file=phasefn, status='OLD')
  open(98, file=genofn, status='UNKNOWN')
  i = 0
  do while (.TRUE.)
    read(97, *, iostat=stat) animalid, linea
    if (stat /= 0) exit
    read(97, *, iostat=stat) animalid, lineb
    summ = linea + lineb
    where (summ > 2.000 .or. summ < 0.000) summ = naval
    write(98, fmt) animalid, summ
    i = i + 1
    if (i == nrow) exit
  end do

  close(97)
  close(98)

  nrow = i

end subroutine

subroutine convert_phase_int(phasefn, genofn, ncol, nrow, naval, numfmt, lennumfmt)

  implicit none

  !! Arguments
  character(255), intent(in) :: phasefn, genofn
  character(len=lennumfmt) :: numfmt
  integer, intent(in) :: ncol, naval, lennumfmt
  integer, intent(inout) :: nrow

  !! Local variables
  integer :: stat, i, animalid
  integer, dimension(ncol) :: linea, lineb, summ
  character(100) :: nChar, fmt

  write(nChar,*) ncol
  fmt='(i20,'//trim(adjustl(nChar))//'i2)'

  open(97, file=phasefn, status='OLD')
  open(98, file=genofn, status='UNKNOWN')
  i = 0
  do while (.TRUE.)
    read(97, *, iostat=stat) animalid, linea
    if (stat /= 0) exit
    read(97, *, iostat=stat) animalid, lineb
    summ = linea + lineb
    where (summ > 2 .or. summ < 0) summ = naval
    write(98, fmt) animalid, summ
    i = i + 1
    if (i == nrow) exit
  end do

  close(97)
  close(98)

  nrow = i

end subroutine

!! Converting plink --recode A

subroutine convertplinka(rawfn, outputfn, newID, ncol, nrow, naval, stat) 

  implicit none
  
  !! Arguments
  character(255), intent(in) :: rawfn, outputfn
  integer, intent(in) :: ncol, nrow, naval
  integer, dimension(nrow) :: newID
  integer, intent(out) :: stat
  
  !! Local variables
  integer :: i
  integer, dimension(ncol) :: iSNPs
  character(3) :: chrNA
  character(100) :: nChar, fmt
  character(3), dimension(ncol) :: SNPs
  character(5), dimension(6) :: sixcolumns

  !! Set output format
  write(nChar,*) ncol
  fmt='(i20,'//trim(adjustl(nChar))//'i2)'
  
  write(chrNA, '(I2)') naval
  
  i = 0
  open(90, file=rawfn, status='OLD')
  open(91, file=outputfn, status='UNKNOWN')
  do while (.TRUE.)
    read(90, *, iostat=stat) sixcolumns, SNPs
    if (stat /= 0) exit
    i = i + 1
    where(SNPs == 'NA') SNPs = chrNA
    read(SNPs, '(I2)') iSNPs
    write(91, fmt) newID(i), iSNPs
    if (i == nrow) exit
  end do  
  close(90)
  close(91)
  
  stat = i

end subroutine convertplinka
