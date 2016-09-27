!! Calculation of heterozygosity
!! Sauce: http://www.uwyo.edu/dbmcd/molmark/practica/fst.html
subroutine heterozygosity(fn, ncols, NAval, p, Hobs, Hexp, n)
  implicit none

  integer, parameter :: r8_kind = selected_real_kind(15, 307)

  !! Arguments
  character(255), intent(in) :: fn
  integer, intent(in) :: ncols, NAval
  integer, dimension(ncols), intent(out) :: n
  real(r8_kind), dimension(ncols), intent(out) :: p, Hobs, Hexp

  !! Local variables
  integer :: i, stat, animalid, j
  integer, dimension(ncols) :: genotype
  integer, dimension(ncols, 0:2) :: allelecount

  allelecount(:,:) = 0
  n(:) = 0

  open(45, file=fn, status='OLD')
  i = 0
  do while (.TRUE.)
    read(45, *, iostat=stat) animalid, genotype
    if (stat /= 0) exit
    do j=1,ncols
      if (genotype(j) == NAval) cycle
      allelecount(j,genotype(j)) = allelecount(j,genotype(j)) + 1
      n(j) = n(j) + 1
    end do
    i = i + 1
  end do
  close(45)

  p(:) = (allelecount(:,1) + allelecount(:,0)*2) / (2. * n(:) )
  q(:) = 1 - p

  Hobs(:) = allelecount(:,1) / (1. * n(:))
  Hexp(:) = 2.*p*q


end subroutine



