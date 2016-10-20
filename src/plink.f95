subroutine readplink(ped, bim, bed, fnout, ncol, nlines, na)
  
  implicit none
  
  character(255), intent(in) :: ped, bim, bed, fnout
  integer, intent(in) :: ncol, nlines, na
  
  !! Local arguments
  byte :: plinkmode, element
  byte, dimension(2) :: magicnumber
  integer :: stat, i, j, pos
  integer, dimension(:), allocatable :: snps
  
  allocate(snps(ncol*nlines))
  
  snps(:) = -1
  
  open(15, file=bed, status='OLD', ACCESS='STREAM', FORM='UNFORMATTED')
  read(15) magicnumber, plinkmode
  if (plinkmode == 1) then
    print *, 'SNP major-mode.'
  endif
  j=1
  do 
    read(15, iostat=stat) element
    if (stat /= 0) exit
    print '(B10.8,B10.8,I3)', element, not(element), BIT_SIZE(element)
    print *, IBITS(element, 0, 2), IBITS(element, 2, 2), IBITS(element, 4, 2), IBITS(element, 6, 2)
    do i=0,6,2
      select case(IBITS(element, i, 2))
        case (0) ! homozygote
          snps(j+i/2) = 0
        case (1) ! missing
          snps(j+i/2) = na
        case (2) ! heterozygote
          snps(j+i/2) = 1
        case (3) ! homozygote, minor
          snps(j+i/2) = 2
      endselect
      print '(I2, I2, B3.2, I2)', j+i/2, ibits(element, i, 2), ibits(element, i, 2), snps(j+i/2)
    enddo
    j=j+i/2
  enddo
  close(15)
  
  print '(I2)', snps
  print *, 'farewell'
  
  deallocate(snps)
  
end subroutine