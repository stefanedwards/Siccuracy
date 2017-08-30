!! Auxiliary functions
!!
!! Get number of lines (Fortran style, but you might just do cat fn | wc -l).
!! Read first column of a file (i.e. animals)
!! ###Read first row. (Do it by R...).


!! Get number of lines
subroutine get_nlines(fn, nlines, stat)
  implicit none

  !! Arguments
  character(255), intent(in) :: fn
  integer, intent(out) :: nlines, stat

  !! Local variables
  character(len=1) :: one

  nlines = 0
  open(40, file=fn, status='OLD')
  do
    read(40, *, iostat=stat) one
    if (stat /= 0) exit
    nlines = nlines + 1
  end do
  close(40)

end subroutine

!! Get integers from first column
!! Date: September 3rd 2015
subroutine get_firstcolumn(fn, nlines, column, stat)
  !USE ISO_FORTRAN_ENV
  implicit none
  
  integer, parameter :: i16_kind = selected_int_kind(16)

  !! Arguments
  character(255), intent(in) :: fn
  integer, intent(in) :: nlines
  integer, dimension(nlines), intent(out) :: column
  integer, intent(out) :: stat

  !! Private variables
  integer(i16_kind) :: i, dummy

  open(30, file=fn, status='OLD')
  do i=1,nlines
    read(30, *, iostat=stat) dummy
    if (stat /= 0) exit
    column(i) = dummy
  end do
  close(30)
end subroutine

