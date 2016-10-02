! Cannot pass character vector from R, so we use this wrapper.
subroutine cbind_SNPs_Rwrap(files, fnin, fnout, nlines, ncols, skiplines, idlength, excludeids, result, lenfmt, userfmt)

  integer, intent(in) :: files, nlines, skiplines, idlength, lenfmt
  character(lenfmt), intent(in) :: userfmt
  character(255), intent(in) :: fnin, fnout
  integer, dimension(files), intent(in) :: ncols
  integer, intent(out) :: result
  integer, dimension(idlength), intent(in) :: excludeids

  integer :: i, stat
  character(255), dimension(files) :: fns

  open(64, file=fnin, status='OLD')
  do i=1,files
    read(64, *, iostat=stat) fns(i)
  end do
  close(64)

  call cbind_SNPs(files, fns, fnout, nlines, ncols, skiplines, idlength, excludeids, result, lenfmt, userfmt)

end subroutine cbind_SNPs_Rwrap


! Concatenates genotype matrices (from e.g. multiple chromosomes) into one.
! No row ID checking.
subroutine cbind_SNPs(files, fns, fnout, nlines, ncols, skiplines, idlength, excludeids, result, lenfmt, userfmt)
  implicit none

  integer, intent(in) :: files, nlines, skiplines, idlength, lenfmt
  character(lenfmt), intent(in) :: userfmt
  character(255), intent(in) :: fnout
  character(255), dimension(files), intent(in) :: fns
  integer, dimension(files), intent(in) :: ncols
  integer, dimension(idlength), intent(in) :: excludeids
  integer, intent(out) :: result

  integer :: i, j, id, stat, k
  integer, dimension(files) :: units
  real, dimension(:,:), allocatable :: row
  character(50) :: fmt0
  character(50), dimension(files) :: fmt
  character(4), dimension(files) :: advance

  allocate(row(files,maxval(ncols, 1)))

  advance(:) = 'no'
  advance(files) = 'yes'

  do i=1,files
    units(i) = 200 + i
    open(units(i), file=fns(i), status='OLD')
    write(fmt0, '(i5)') ncols(i)
    !fmt(i)='('//trim(adjustl(fmt0))//userfmt//')'
    fmt(i)='('//trim(adjustl(fmt0))//'F3)'
  end do
  open(55, file=fnout, status='UNKNOWN')

  do j=1,nlines
    fileloop: do i=1,files
      read(units(i), *, iostat=stat) id, row(i,1:ncols(i))
      if (stat /= 0) exit
      if (j <= skiplines) cycle
      if (idlength > 0) then
        do k=1, idlength
          if (excludeids(k) .eq. id) cycle fileloop
        enddo
      endif
      if (i == 1) write(55, '(i20)', advance='no') id
      write(55, fmt(i), advance=advance(i)) row(i,1:ncols(i))
    end do fileloop
    if (stat /= 0) exit
  end do

  deallocate(row)

  close(55)
  do i=1,files
    close(units(i))
  end do

  result=stat

end subroutine cbind_SNPs