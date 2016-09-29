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
  print *,'I did something:',fmt

  open(97, file=phasefn, status='OLD')
  open(98, file=genofn, status='UNKNOWN')
  i = 0
  do while (.TRUE.)
    read(97, *, iostat=stat) animalid, linea
    if (stat /= 0) exit
    read(97, *, iostat=stat) animalid, lineb
    summ = linea + lineb
    where (summ >= naval) summ = naval
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
    where (summ >= naval)  summ = naval
    write(98, fmt) animalid, summ
    i = i + 1
    if (i == nrow) exit
  end do

  close(97)
  close(98)

  nrow = i

end subroutine