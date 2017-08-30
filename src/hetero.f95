!! Calculation of heterozygosity
!! Sauce: http://www.uwyo.edu/dbmcd/molmark/practica/fst.html
subroutine heterozygosity(fn, ncols, nlines, populations, npop, &
           p, Hobs, Hexp, n)
  implicit none

  integer, parameter :: r8_kind = selected_real_kind(15, 307)
  integer, parameter :: i16_kind = selected_int_kind(16)

  !! Arguments
  character(255), intent(in) :: fn
  integer, intent(in) :: ncols, nlines, npop
  integer, dimension(ncols,npop), intent(out) :: n
  integer, dimension(nlines), intent(inout) :: populations
  real(r8_kind), dimension( ncols, npop), intent(out) :: p, Hobs, Hexp

  !! Local variables
  integer :: i, stat, j
  integer(i16_kind) :: animalid
  integer, dimension(ncols) :: genotype
  integer, dimension(0:2,ncols, npop) :: allelecount

  allelecount(:,:,:) = 0
  n(:,:) = 0

  open(45, file=fn, status='OLD')
  i = 0
  do while (.TRUE.)
    read(45, *, iostat=stat) animalid, genotype
    if (stat /= 0) exit
    i = i + 1
    do j=1,ncols
      if (genotype(j) < 0 .or. genotype(j) > 2) cycle
      allelecount(genotype(j),j,populations(i)) = allelecount(genotype(j),j,populations(i)) + 1
      n(j,populations(i)) = n(j,populations(i)) + 1
    end do
  end do
  close(45)

  p(:,:) = (allelecount(1,:,:) + allelecount(0,:,:)*2) / (2. * n(:,:) )

  Hobs(:,:) = allelecount(1,:,:) / (1. * n(:,:))
  Hexp(:,:) = 2.*p*(1-p)

end subroutine




