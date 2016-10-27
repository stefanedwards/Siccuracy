
! Cannot pass character vector from R, so we use this wrapper.
subroutine cbindsnpsrwrapper(files, fnin, fnout, nlines, ncols, skiplines, idlength, excludeids, result, lenfmt, userfmt, asint)
  implicit none

  integer, intent(in) :: files, nlines, skiplines, idlength, lenfmt, asint
  character(lenfmt), intent(in) :: userfmt
  character(255), intent(in) :: fnin, fnout
  integer, dimension(files), intent(in) :: ncols
  integer, intent(out) :: result
  integer, dimension(idlength), intent(in) :: excludeids

  logical :: alsoasint
  integer :: i, stat
  character(255), dimension(files) :: fns

  open(64, file=fnin, status='OLD')
  do i=1,files
    read(64, '(A255)', iostat=stat) fns(i)
  end do
  close(64)

  alsoasint = asint==1
  
  call cbindsnps(files, fns, fnout, nlines, ncols, skiplines, idlength, excludeids, result, lenfmt, userfmt, alsoasint)

end subroutine cbindsnpsrwrapper


! Concatenates genotype matrices (from e.g. multiple chromosomes) into one.
! No row ID checking.
subroutine cbindsnps(files, fns, fnout, nlines, ncols, skiplines, idlength, excludeids, result, lenfmt, userfmt, asint)
  implicit none

  logical, intent(in) :: asint
  integer, intent(in) :: files, nlines, skiplines, idlength, lenfmt
  character(lenfmt), intent(in) :: userfmt
  character(255), intent(in) :: fnout
  character(255), dimension(files), intent(in) :: fns
  integer, dimension(files), intent(in) :: ncols
  integer, dimension(idlength), intent(in) :: excludeids
  integer, intent(out) :: result

  integer :: i, j, id, stat, k
  integer, dimension(files) :: units
  integer, dimension(:), allocatable :: rowint
  real, dimension(:), allocatable :: rowreal
  character(50) :: fmt0
  character(50), dimension(files) :: fmt
  character(4), dimension(files) :: advance, formatprefix

  if (asint .eqv. .true.) then
    allocate(rowint(maxval(ncols, 1)))
  end if
  allocate(rowreal(maxval(ncols, 1)))

  advance(:) = 'no'
  advance(files) = 'yes'
  
  do i=1,files
    units(i) = 200 + i
    open(units(i), file=fns(i), status='OLD')
    write(fmt0, '(i5)') ncols(i)
    fmt(i)='('//trim(adjustl(fmt0))//trim(adjustl(userfmt))//')'
  end do
  open(55, file=fnout, status='UNKNOWN')

  do j=1,nlines
    fileloop: do i=1,files
      read(units(i), *, iostat=stat) id, rowreal(1:ncols(i))
      if (stat /= 0) exit
      if (j <= skiplines) cycle
      if (idlength > 0) then
        do k=1, idlength
          if (excludeids(k) .eq. id) cycle fileloop
        enddo
      endif
      if (i == 1) write(55, '(i20)', advance='no') id
      if (asint .eqv. .true.) then
        rowint(1:ncols(i)) = NINT(rowreal(1:ncols(i)))
        write(55, fmt(i), advance=advance(i)) rowint(1:ncols(i))
      else
        write(55, fmt(i), advance=advance(i)) rowreal(1:ncols(i))
      end if
    end do fileloop
    if (stat /= 0) exit
  end do

  if (asint .eqv. .true.) then
    deallocate(rowint)
  end if
  deallocate(rowreal)

  close(55)
  do i=1,files
    close(units(i))
  end do

  result=stat

end subroutine cbindsnps



!! Merges two SNP chips

subroutine rbindsnps(fnhd, fnld, fnout, hdcols, ldcols, outcols, &
  nhd, hdid, nld, ldid, hdpos, ldpos,  missing, lenfmt, userfmt, stat, asint)

  integer, intent(in) :: hdcols, ldcols, outcols, nhd, nld, missing, lenfmt, asint
  character(lenfmt), intent(in) :: userfmt
  character(255), intent(in) :: fnhd, fnld, fnout
  integer, intent(in), dimension(nhd) :: hdid
  integer, intent(in), dimension(nld) :: ldid
  integer, intent(in), dimension(hdcols) :: hdpos
  integer, intent(in), dimension(ldcols) :: ldpos
  integer, intent(out) :: stat

  logical :: foundID, isint
  integer :: i, oldi, animalID
  
  character(50)  :: fmt
  integer, allocatable :: intoutput(:)
  real, allocatable :: line(:), realoutput(:)
  
  isint = asint == 1
  
  write(fmt, '(i5)') outcols
  fmt='(i20,'//trim(adjustl(fmt))//userfmt//')'
  
  if (isint .eqv. .true.) then
    allocate(intoutput(outcols))
    intoutput(:) = missing
  else
    allocate(realoutput(outcols))
    realoutput(:) = missing
  end if


  ! Open output file.
  open(71, file=fnout, status='UNKNOWN', iostat=stat)
  
  ! Read through HD file.
  allocate(line(hdcols))
  i = 0
  oldi = 0
  open(72, file=fnhd, status='OLD')
  do while (.TRUE.)
    read(72, *, iostat=stat) animalID, line
    if (stat /= 0) exit

    ! Find animal in HD IDs.
    foundID = .FALSE.
    !oldi = i
    !if (oldi == 0) oldi=1
    do while (.TRUE.)
      if (i == nhd) i = 0
      i = i + 1
      if (hdid(i) == animalID) then
        foundID = .TRUE.
        exit
      end if
      if (i == oldi) exit
      if (i == nhd .and. oldi == 0) exit
    enddo
    oldi = i
    if (foundID) then
      if (isint .eqv. .true.) then
        intoutput(hdpos) = NINT(line(1:hdcols))
        write(71, fmt) animalID, intoutput
      else
        realoutput(hdpos) = line(1:hdcols)
        write(71, fmt) animalID, realoutput
      end if
    end if

  end do
  close(72)
  deallocate(line)

  ! Read through LD file.
  allocate(line(ldcols))
  if (isint .eqv. .true.) then
    intoutput(:) = missing
  else
    realoutput(:) = missing
  endif
  
  i = 0
  oldi=0
  open(72, file=fnld, status='OLD')
  do while (.TRUE.)
    read(72, *, iostat=stat) animalID, line
    if (stat /= 0) exit


    ! Find animal in LD IDs.
    foundID = .FALSE.
    oldi = i
    if (oldi == 0) oldi=1
    do while (.TRUE.)
      if (i == nld) i = 0
      i = i + 1
      if (ldid(i) == animalID) then
        foundID = .TRUE.
        exit
      end if
      if (i == oldi) exit
      if (i == nld .and. oldi == 0) exit
    enddo
    oldi = i
    if (foundID) then
      if (isint .eqv. .true.) then
        intoutput(ldpos) = NINT(line(1:ldcols))
        write(71, fmt) animalID, intoutput
      else
        realoutput(ldpos) = line(1:ldcols)
        write(71, fmt) animalID, realoutput
      end if
    end if

  end do
  close(72)
  close(71)
  
  deallocate(line)

  if (isint .eqv. .true.) then
    deallocate(intoutput)
  else
    deallocate(realoutput)
  endif


end subroutine rbindsnps
