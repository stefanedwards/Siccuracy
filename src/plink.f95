! Stores entire genotype file in memory
subroutine readplinksimple(bed, fnout, ncol, nlines, na, newID, minor, maf, extract, keep, status)
  
  implicit none
  
  character(255), intent(in) :: bed, fnout
  integer, intent(in) :: ncol, nlines, na, minor
  integer, dimension(nlines), intent(in) :: newID, keep
  integer, dimension(ncol), intent(in) :: extract
  integer, intent(out) :: status
  double precision, intent(in) :: maf
  
  !! Local arguments
  byte :: readplinkmode, element, plinkmode
  byte, dimension(2) :: readmagicnumber, magicnumber
  logical :: checkmaf
  logical, dimension(ncol) :: masksnps
  integer :: stat, i, j, pos, k, snpcount, majorcount
  integer, dimension(4) :: codes
  integer, dimension(:), allocatable :: domasksnps
  integer, dimension(:,:), allocatable :: snps
  real :: allelefreq
  character(100) :: nChar, fmt
  
  ! Supported formats as per plink 1.9.
  data magicnumber/X'6c',X'1b'/,  plinkmode/X'01'/
  
  allocate(snps(nlines,ncol))
  snps(:,:) = -1
  
  masksnps=extract==1
  
  if (minor == 1) then
    codes = (/ 0, 1, 2, na /)
  else
    codes = (/ 2, 1, 0, na /)
  endif
  
  open(15, file=bed, status='OLD', ACCESS='STREAM', FORM='UNFORMATTED')
  read(15) readmagicnumber, readplinkmode
  if (all(readmagicnumber /= magicnumber) ) then
    status=-1
    close(15)
    return
  endif
  if (readplinkmode /= plinkmode) then
    status=-2
    close(15)
    return
  endif


  j=0  ! Sample-index
  k=1  ! SNP-index
  snpcount = 0
  majorcount = 0
  outer: do 
    read(15, iostat=stat) element
    if (stat /= 0) exit
    !print '(B10.8,B10.8,I3)', element, not(element), BIT_SIZE(element)
    !print *, IBITS(element, 0, 2), IBITS(element, 2, 2), IBITS(element, 4, 2), IBITS(element, 6, 2)
    inner: do i=0,6,2
      j = j + 1
      snpcount = snpcount + 1
      select case(IBITS(element, i, 2))
        case (0) ! homozygote
          snps(j,k) = codes(1)
        case (1) ! missing
          snps(j,k) = codes(4)
          snpcount = snpcount - 1
        case (2) ! heterozygote
          snps(j,k) = codes(2)
          majorcount = majorcount + 1
        case (3) ! homozygote, minor
          snps(j,k) = codes(3)
          majorcount = majorcount + 2
      endselect
      if (j == nlines) then
        if (snpcount /= 0) then
          allelefreq = majorcount / (snpcount*2.)
          masksnps(k) = masksnps(k) .and. allelefreq .ge. maf .and. allelefreq .le. (1-maf)
        endif
        j = 0
        snpcount = 0
        majorcount = 0
        k = k + 1
        cycle outer
      endif
    enddo inner
  enddo outer
  close(15)
  
  if (stat == -1) stat=0
  
  ! Write output
  write(nChar,*) count(masksnps)
  fmt='(i20,'//trim(adjustl(nChar))//'I2)'
  
  open(16, file=fnout)
  do i=1,nlines
    if (keep(i) == 1) write(16, fmt) newID(i), pack(snps(i,:), masksnps)
  enddo
  close(16)

  deallocate(snps)
  
end subroutine readplinksimple


!! Fragmented approach.
! First an R wrapper for retrieving filenames,
! then the actual subroutine (that calls readplinksimple)

subroutine convertplinkrwrapper(listfn, n, remerge, fragments, &
    bed, fnout, ncol, nlines, na, newID, minor, maf, extract, keep, status)
    ! 2nd line are readplinksimple arguments.
   
  implicit none
  
  ! Arguments
  character(255), intent(in) :: bed, fnout, listfn
  integer, intent(in) :: ncol, nlines, na, minor, remerge, n
  integer, dimension(nlines), intent(in) :: newID, keep
  integer, dimension(ncol), intent(in) :: extract, fragments
  integer, intent(out) :: status
  double precision, intent(in) :: maf
  
  ! Local variables
  character(255), dimension(n) :: filelist
  integer :: i
  
  open(77, file=listfn)
  read(77, *) (filelist(i), i=1, n)
  close(77)
  
  call convertplinkfragment(filelist((n/2+1):n), filelist(1:(n/2)), n/2, remerge==1, fragments,   &
      bed, fnout, ncol, nlines, na, newID, minor, maf, extract, keep, status)

end subroutine convertplinkrwrapper

subroutine convertplinkfragment(bedfilenames, flatfilenames, n, remerge, fragments, &
    bed, fnout, ncol, nlines, na, newID, minor, maf, extract, keep, status)

  implicit none
  ! Arguments
  logical, intent(in) :: remerge
  character(255), intent(in) :: bed, fnout
  character(255), dimension(n), intent(in) :: bedfilenames, flatfilenames
  integer, intent(in) :: ncol, nlines, na, minor, n
  integer, dimension(nlines), intent(in) :: newID, keep
  integer, dimension(ncol), intent(in) :: extract, fragments
  integer, intent(out) :: status
  double precision, intent(in) :: maf
  
  ! Local variables
  byte :: readplinkmode, element, plinkmode
  byte, dimension(2) :: readmagicnumber, magicnumber
  integer :: i,j,k,stat,ncoli
  integer, dimension(n) :: bedcons, flatcons
  integer, dimension(ncol) :: allpos
  integer, dimension(:), allocatable :: subset
  logical, dimension(ncol) :: mask
  

  ! Supported formats as per plink 1.9.
  data magicnumber/X'6c',X'1b'/,  plinkmode/X'01'/
  
  open(15, file=bed, status='OLD', ACCESS='STREAM', FORM='UNFORMATTED')
  read(15) readmagicnumber, readplinkmode
  if (all(readmagicnumber /= magicnumber) ) then
    status=-1
    close(15)
    return
  endif
  if (readplinkmode /= plinkmode) then
    status=-2
    close(15)
    return
  endif
  
  do i=1,n
    ! Delete file first.
    ! Opening for binary output does not truncate the file, so any previous data will be left trailing...
    open(100+i, file=bedfilenames(i), status='old', iostat=stat)
    if (stat == 0) close(100+i, status='delete')
    ! Then open for writing.
    open(100+i, file=bedfilenames(i), status='UNKNOWN', ACCESS='STREAM', form='UNFORMATTED')
    write(100+i) magicnumber, plinkmode
  enddo

  i = 1
  j = 1
  outer: do
    read(15, iostat=stat) element
    if (stat /= 0) exit
    write(100+fragments(i)) element
    j = j + 4
    if (j .ge. nlines) then
      j = 1
      i = i+1
    endif
  enddo outer

  close(15)
  do i=1,n
    close(100+i)
  enddo
  
  forall(i=1:ncol) allpos(i)=i
  do i=1,n
    mask=fragments==i
    ncoli=count(mask)
    allocate(subset(ncoli))
    subset = extract(pack(allpos,mask))
        !readplinksimple(bed, fnout, ncol, nlines, na, newID, minor, maf, extract, keep, status)
    call readplinksimple(bedfilenames(i), flatfilenames(i), ncoli, nlines, na, newID, minor, maf, subset, keep, stat)
    deallocate(subset)
  enddo
end subroutine 