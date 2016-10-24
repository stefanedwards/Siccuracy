! Stores entire genotype file in memory
subroutine readplinksimple(bed, fnout, ncol, nlines, na, newID, minor, status)
  
  implicit none
  
  character(255), intent(in) :: bed, fnout
  integer, intent(in) :: ncol, nlines, na, minor
  integer, dimension(nlines), intent(in) :: newID
  integer, intent(out) :: status
  
  !! Local arguments
  byte :: readplinkmode, element, plinkmode
  byte, dimension(2) :: readmagicnumber, magicnumber
  integer :: stat, i, j, pos, k
  integer, dimension(4) :: codes
  integer, dimension(:,:), allocatable :: snps
  character(100) :: nChar, fmt
  
  ! Supported formats as per plink 1.9.
  data magicnumber/X'6c',X'1b'/,  plinkmode/X'01'/
  
  allocate(snps(nlines,ncol))
  snps(:,:) = -1
  
  if (minor == 1) then
    codes = (/ 0, 1, 2, na /)
  else
    codes = (/ 2, 1, 0, na /)
  endif
  
  open(15, file=bed, status='OLD', ACCESS='STREAM', FORM='UNFORMATTED')
  read(15) readmagicnumber, readplinkmode
  if (all(readmagicnumber /= magicnumber) ) then
    status=-1
    return
  endif
  if (readplinkmode /= plinkmode) then
    status=-2
    return
  endif


  j=0  ! SNP-index
  k=1  ! Sample-index
  outer: do 
    read(15, iostat=stat) element
    if (stat /= 0) exit
    !print '(B10.8,B10.8,I3)', element, not(element), BIT_SIZE(element)
    !print *, IBITS(element, 0, 2), IBITS(element, 2, 2), IBITS(element, 4, 2), IBITS(element, 6, 2)
    inner: do i=0,6,2
      j = j + 1
      select case(IBITS(element, i, 2))
        case (0) ! homozygote
          snps(j,k) = codes(1)
        case (1) ! missing
          snps(j,k) = codes(4)
        case (2) ! heterozygote
          snps(j,k) = codes(2)
        case (3) ! homozygote, minor
          snps(j,k) = codes(3)
      endselect
      if (j == nlines) then
        j = 0
        k = k + 1
        cycle outer
      endif
    enddo inner
  enddo outer
  close(15)
  
  if (stat == -1) stat=0
  
  ! Write output
  write(nChar,*) ncol
  fmt='(i20,'//trim(adjustl(nChar))//'I2)'
  
  open(16, file=fnout)
  do i=1,nlines
    write(16, fmt) newID(i), snps(i,:)
  enddo
  close(16)
  
  deallocate(snps)
  
end subroutine readplinksimple