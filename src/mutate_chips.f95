
! Cannot pass character vector from R, so we use this wrapper.
subroutine cbindsnpsrwrapper(files, fnin, fnout, nlines, ncols, skiplines, & 
    idlength, excludeids, result, lenfmt, userfmt, asint)
  implicit none

  integer, intent(in) :: files, nlines, skiplines, idlength, lenfmt, asint
  character(lenfmt), intent(in) :: userfmt
  character(255), intent(in) :: fnin, fnout
  integer, dimension(files), intent(in) :: ncols
  integer, intent(out) :: result
  integer, dimension(idlength), intent(in) :: excludeids

  logical :: alsoasint
  integer :: i, stat
  !character(1000) :: ctag
  character(1000), dimension(files) :: fns
  
  open(64, file=fnin, status='OLD')
  do i=1,files
    read(64, '(A255)', iostat=stat) fns(i)
  end do
  close(64)


  !open(27, file='C:\Users\shojedw\Documents\Projects\test1.txt')
  !write(27, *) 'tag:', ctag
  !write(27, '(A7,A255)') 'input:', trim(adjustl(fnin))
  !write(27, *) 'files:', files
  !do i=1,files
  !  write(27, *) 'file:', i, '>>'//trim(adjustl(fns(i)))//'<<'
  !enddo
  !close(27)
  
  alsoasint = asint==1
  
  !ctag="Hello, is it me you're looking for?"
  call cbindsnpscore(files, fns, fnout, nlines, ncols, skiplines, &
      idlength, excludeids, result, lenfmt, userfmt, alsoasint)
  
  
end subroutine cbindsnpsrwrapper


! Concatenates genotype matrices (from e.g. multiple chromosomes) into one.
! No row ID checking.
! The argument ctag isn't really used for anything. It just somehow
! fixes a bug where fns is not populated.
subroutine cbindsnpscore(files, fns, fnout, nlines, ncols, skiplines, &
    idlength, excludeids, result, lenfmt, userfmt, asint)
  implicit none

  logical, intent(in) :: asint
  integer, intent(in) :: files, nlines, skiplines, idlength, lenfmt
  character(lenfmt), intent(in) :: userfmt
  !character(1000), intent(in) :: ctag
  character(255), intent(in) :: fnout
  character(1000), dimension(files), intent(in) :: fns
  integer, dimension(files), intent(in) :: ncols
  integer, dimension(idlength), intent(in) :: excludeids
  integer, intent(out) :: result

  integer :: i, j, id, stat, k
  integer, dimension(files) :: units
  integer, dimension(:), allocatable :: rowint
  real, dimension(:), allocatable :: rowreal
  character(50) :: fmt0
  character(50), dimension(files) :: fmt
  character(4), dimension(files) :: advance
  
  if (asint .eqv. .true.) then
    allocate(rowint(maxval(ncols, 1)))
  end if
  allocate(rowreal(maxval(ncols, 1)))

  advance(:) = 'no'
  advance(files) = 'yes'

  !open(27, file='C:\Users\shojedw\Documents\Projects\test2.txt')  
  !write(27, *) 'tag:', ctag  
  !write(27, *) 'from cbindsnps:'
  !write(27, *) 'files:', files
  !do i=1,files
  !  write(27, '(A6,I2,A70,A70)') 'file:', i, '>>'//trim(adjustl(fns(i)))//'<< ', fns(i)
  !enddo
  !close(27)
  
  do i=1,files
    units(i) = 200 + i
    open(units(i), file=fns(i), status='UNKNOWN')
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

end subroutine cbindsnpscore



!! Merges two SNP chips

subroutine rbindsnps(fnhd, fnld, fnout, hdcols, ldcols, outcols, &
  nhd, hdid, nld, ldid, hdpos, ldpos,  missing, lenfmt, userfmt, asint, stat, totallines)

  implicit none
  
  integer, parameter :: i16_kind = selected_int_kind(16)

  integer, intent(in) :: hdcols, ldcols, outcols, nhd, nld, missing, lenfmt, asint
  character(lenfmt), intent(in) :: userfmt
  character(255), intent(in) :: fnhd, fnld, fnout
  integer, intent(in), dimension(nhd) :: hdid
  integer, intent(in), dimension(nld) :: ldid
  integer, intent(in), dimension(hdcols) :: hdpos
  integer, intent(in), dimension(ldcols) :: ldpos
  integer, intent(out) :: stat, totallines

  logical :: foundID, isint
  integer :: i, oldi
  integer(i16_kind) :: animalID
  
  character(50)  :: fmt
  integer, allocatable :: intoutput(:)
  real, allocatable :: line(:), realoutput(:)
  logical, dimension(:), allocatable :: mask
  
  isint = asint == 1
  
  write(fmt, '(i20)') outcols
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
  allocate(mask(hdcols))
  mask(:) = hdpos(:) /= 0
  
  i = 0
  oldi = 0
  totallines=0
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
        intoutput(hdpos) = NINT(PACK(line, mask))
        write(71, fmt) animalID, intoutput
      else
        realoutput(hdpos) = PACK(line, mask)
        write(71, fmt) animalID, realoutput
      end if
      totallines = totallines+1
    end if

  end do
  close(72)
  deallocate(line)
  deallocate(mask)

  ! Read through LD file.
  allocate(line(ldcols))
  allocate(mask(ldcols))
  mask(:) = ldpos(:) /= 0
  
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
        intoutput(ldpos) = NINT(PACK(line, mask))
        write(71, fmt) animalID, intoutput
      else
        realoutput(ldpos) = PACK(line, mask)
        write(71, fmt) animalID, realoutput
      end if
      totallines = totallines+1
    end if

  end do
  close(72)
  close(71)
  
  deallocate(line)
  deallocate(mask)

  if (isint .eqv. .true.) then
    deallocate(intoutput)
  else
    deallocate(realoutput)
  endif

  if (stat == -1) stat = 0

end subroutine rbindsnps


!! Mask SNPs

subroutine masksnps2(fn, outfn, ncols, nlines, na, userfmt, lenuserfmt, asint, stat, &
    dropIDs, snpslength, snps, maskIDs, maps, maskstart, maskend, masklength, maskSNPs)

  implicit none

  integer, parameter :: i16_kind = selected_int_kind(16)

  ! Arguments
  character(255), intent(in) :: fn, outfn
  integer, intent(in) :: ncols, nlines, lenuserfmt, asint, snpslength, maps, masklength
  character(lenuserfmt), intent(in) :: userfmt
  double precision, intent(in) :: na
  integer, intent(out) :: stat
  integer, dimension(snpslength), intent(in) :: snps
  integer, dimension(nlines), intent(in) :: maskIDs, dropIDs
  integer, dimension(maps), intent(in) :: maskstart, maskend
  integer, dimension(masklength), intent(in) :: maskSNPs

    
  ! Local variables
  logical :: isint
  integer :: i, ndropSNPs, snpstocopy
  integer(i16_kind) :: animalID 
  character(50) :: fmt
  logical, dimension(:), allocatable :: isnotzero
  integer, dimension(:), allocatable :: intoutput, copysnps, copysnpsdest
  real, dimension(:), allocatable :: realoutput, inputline, outputline

  ! Construct input and output holders
  allocate(inputline(1:ncols)) 
  allocate(isnotzero(1:snpslength))

  isnotzero =  snps /= 0
  !print *, isnotzero
  allocate(copysnps(1:COUNT(isnotzero)))
  allocate(copysnpsdest(1:COUNT(isnotzero)))
  copysnps = PACK(snps, isnotzero)
  copysnpsdest = PACK( (/ (i, i=1, snpslength) /), isnotzero)
  
  !print *, 'snps:', snps
  !print *, 'copysnps:', copysnps
  !print *, 'copysnp dest:', copysnpsdest
  
  isint = asint == 1
  if (isint) then
    allocate(intoutput(snpslength))  
  endif
  
  allocate(realoutput(snpslength))
  realoutput(:) = na
  
  
  !print *, 'maskIDs:', maskIDs
  !print *, 'maskSNPs:', maskSNPs
  !print *, 'dropIDs:', dropIDs
  
  write(fmt, '(i10)') snpslength
  fmt='(i20,'//trim(adjustl(fmt))//userfmt//')'
  
  ! Loop: read, write, loop.
  open(119, file=fn, status='OLD')
  open(120, file=outfn, status='UNKNOWN')
  do i=1,nlines
    if (dropIDs(i) == 1) then
      read(119, *, iostat=stat) animalID
      cycle
    endif
    
    read(119, *, iostat=stat) animalID, inputline
    if (stat /= 0) exit
    !print *, animalID, inputline

    if (maskIDs(i) /= 0) then 
      inputline(maskSNPs(maskstart(maskIDs(i)):maskend(maskIDs(i)))) = na
    endif

    
    realoutput(copysnpsdest) = inputline(copysnps)
    !print *, animalID, realoutput
    
    
    if (isint) then
      intoutput = NINT(realoutput)
      write(120, fmt) animalID, intoutput
    else
      write(120, fmt) animalID, realoutput
    endif
    
    
  enddo
  close(119)
  close(120)
  
  
  deallocate(inputline)
  deallocate(realoutput)
  if (isint) then
    deallocate(intoutput)
  endif
  deallocate(copysnps)
  deallocate(copysnpsdest)
  
  if (stat == -1) stat=0
end subroutine    masksnps2
    
    