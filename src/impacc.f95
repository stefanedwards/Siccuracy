!! imputation accuracy

!! Calculates imputation accuracy
!! truefn and imputefn are filepaths of respectively true genotypes and imputed genotypes
!! m, number of SNPs
!! NAval is an integer value of missing genotypes, which are ignored.
!! standardized, boolean, whether to standardize genotypes based on entire true genotypes.
!! rowcor, matcor, colcor, vector of results.
subroutine imp_acc_fast(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, means, vars, 
    rowcors, matcor, colcors, rowcorID)
  implicit none
  
  integer, parameter :: r8_kind = selected_real_kind(6, 37) ! 32-bit 

  !! Arguments
  character(255), intent(in) :: truefn, imputefn
  integer, intent(in) :: nSNPs, NAval, standardized, nAnimals
  !double precision, dimension(nSNPs), intent(out) :: colcor
  !double precision, dimension(nAnimals), intent(out) :: rowcor
  !double precision, intent(out) :: matcor
  real(r8_kind), dimension(nSnps), intent(out) :: means, vars, colcors
  real(r8_kind), dimension(nAnimals), intent(out) :: rowcors
  real(r8_kind), intent(out) :: matcor

  !! Private variables
  integer :: stat, animalID, nLines, i, j, nn
  real(r8_kind) :: t, t2, imp, imp2, tim, nan
  real(r8_kind), dimension(nSNPs) :: M, S, Mold, Sold
  real, dimension(nSNPs) :: genoin, true, imputed
  !! For running correlation on columns
  integer, dimension(nSNPs) :: cNA, nnn
  real(r8_kind), dimension(nSNPs) :: cx, cy, cx2, cxy, cy2
  !! For matrix
  real(r8_kind) :: mx, my, mx2, mxy, my2
  !! For row
  integer :: rNA
  real(r8_kind) :: rx, ry, rx2, rxy, ry2

  !! NAN
  !! Quiet NAN, double precision.
  !! Source, Steve Lionel (Intel) https://software.intel.com/en-us/forums/topic/294680
  nan = 0.

  !! Read through true genotype file and get column-wise mean and variance
  if (standardized == 1) then
    open(10, file=truefn, status='OLD')
    read(10, *, iostat=stat) animalID, genoin
    nLines = 1
    M(:) = genoin(:)
    S(:) = 0
    do
      read(10, *, iostat=stat) animalID, genoin
      if (stat /= 0) exit
 
      nLines = nLines + 1
      Mold(:) = M(:)
      Sold(:) = S(:)
      M = M + (genoin - Mold)/nLines
      S = S + (genoin - Mold) * (genoin - M)

    enddo
    close(10)
    means = M
    vars = S/(nLines - 1)
    !where (vars == 0.0) vars = 1.0
  else
    means(:) = 0
    vars(:) = 1
  end if

  !! Go through both files
  open(10, file=truefn, status='OLD')
  open(20, file=imputefn, status='OLD')
  i = 0
  mx=0; my=0; mx2=0; mxy=0; my2=0
  cx(:)=0; cy(:)=0; cx2(:)=0; cxy(:)=0; cy2(:)=0; cNA(:) = 0
  do
    i = i + 1
    read(10, *, iostat=stat) animalID, true
    if (stat /= 0) then
      exit
    endif
    read(20, *, iostat=stat) animalID, imputed
    if (stat /= 0) then
      exit
    endif

    rx = 0
    ry = 0
    rx2 = 0
    rxy = 0
    ry2 = 0
    rNA = 0

    do j=1,nSnps
      if (imputed(j) == NAval .or. true(j) == NAval .or. vars(j) == 0.) then
        rNA = rNA + 1
        cNA(j) = cNA(j) + 1
        !print *, animalID, j
        cycle
      end if

      t = (true(j)-means(j))/vars(j)
      imp = (imputed(j)-means(j))/vars(j)
      t2 = t*t
      imp2 = imp*imp
      tim = t*imp

      cx(j) = cx(j) + t
      mx = mx + t
      rx = rx + t

      cy(j) = cy(j) + imp
      my = my + imp
      ry = ry + imp

      cx2(j) = cx2(j) + t2
      mx2 = mx2 + t2
      rx2 = rx2 + t2

      cy2(j) = cy2(j) + imp2
      my2 = my2 + imp2
      ry2 = ry2 + imp2

      cxy(j) = cxy(j) + tim
      mxy = mxy + tim
      rxy = rxy + tim
    enddo
    nn = j - 1 - rNA
    !print *, animalID, nn, rx, ry, rx2, ry2, rxy, dummy
    rowcors(i) = (nn * rxy - rx * ry) / sqrt( (nn * rx2 - rx**2) * (nn*ry2 - ry**2 ) )
  enddo
  close(10)
  close(20)

  nn = (i-1) * nSNPs - sum(cNA)
  matcor = (nn * mxy - mx * my) / sqrt( (nn * mx2 - mx**2) * (nn*my2 - my**2) )
  nnn = i - 1  - cNA
  colcors = (nnn * cxy - cx * cy) / sqrt( (nnn * cx2 - cx**2) * (nnn*cy2 - cy**2) )
  if (standardized == 1) where (vars == 0) colcors = 1/nan

end subroutine

!! ============================================================================

!! Imputation accuracy, round 2:
!!
!! In this exciting development, we will be handling uneven true and imputed matrices, meaning that the rows in imputed may be different from
!! those in true. Columns should still be the same!
!! This implementation (unlike imp_acc2) solves this by storing the true matrix, as it shouldn't be too large.
!! If it is, well, go develop imp_acc2.
!!
!! param nAnimals is the number of animals (rows) in the true file.
subroutine imp_acc(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, means, sds, rowcors, matcor, colcors,rowcorID)
  implicit none
  
  integer, parameter :: r8_kind = selected_real_kind(6, 37) ! 32-bit 

  !! Arguments
  character(255), intent(in) :: truefn, imputefn
  integer, intent(in) :: nSNPs, NAval, standardized, nAnimals
  real(r8_kind), dimension(nSnps), intent(out) :: means, sds, colcors
  real(r8_kind), dimension(nAnimals), intent(out) :: rowcors
  integer, dimension(nAnimals), intent(out) :: rowcorID
  real(r8_kind), intent(out) :: matcor

  !! Private variables
  integer :: stat, start, animalID, commonrows, i, j, k, nn, maxanimal, minanimal, ianimalID
  real(r8_kind) :: t, t2, imp, imp2, tim, nan
  real(r8_kind), dimension(nSNPs) :: M, S, Mold, Sold
  real(r8_kind), dimension(nSNPs) :: genoin, imputed
  real(r8_kind), dimension(nAnimals, nSnps) :: trueMat
  !! For running correlation on columns
  integer, dimension(nSNPs) :: cNA, nnn
  integer, dimension(nAnimals) :: animalIndex
  real(r8_kind), dimension(nSNPs) :: cx, cy, cx2, cxy, cy2
  !! For matrix
  real(r8_kind) :: mx, my, mx2, mxy, my2
  !! For row
  integer :: rNA
  real(r8_kind) :: rx, ry, rx2, rxy, ry2

  nan = 0.0

  rowcorID(:) = 0
  rowcors(:) = 0.0
  colcors(:) = 0.0
  
  !! Read in true genotype fil
  open(10, file=truefn, status='OLD')
  i = 1
  do while (.TRUE.)
    read(10, *, iostat=stat) animalID, genoin
    if (stat /= 0) exit
    !print *, 'Read...', animalID
    animalIndex(i) = animalID
    trueMat(i, :) = genoin
    i = i + 1
  end do
  close(10)
  maxanimal = maxval(animalIndex)
  minanimal = minval(animalIndex)
  !print *, 'Min, max:', minanimal, maxanimal
  !print *, 'Index size', size(animalIndex)
  !print *, 'Index:', animalIndex

  !! Calculate scaling parameters (mean and var)
  if (standardized == 1) then
    M(:) = trueMat(1,:)
    S(:) = 0
    do i=2,nAnimals
      Mold(:) = M(:)
      Sold(:) = S(:)
      M = M + (trueMat(i, :) - Mold)/i
      S = S + (trueMat(i, :) - Mold) * (trueMat(i,:) - M)
    end do
    means = M
    sds = sqrt(S/(nAnimals - 1))
  else
    means(:) = 0
    sds(:) = 1
  end if

  !! Sneak peak into imputed
  ianimalID = 0
  open(20, file=imputefn, status='OLD')
  do while (.TRUE.)
    read(20, *, iostat=stat) animalID
    if (stat /= 0) return
    if (animalID > minanimal) then  ! imputed animals has later animals than true.
      rewind(20)
      exit
    end if
    if (ianimalID > animalID) then  ! imputed animals are not ordered
      rewind(20)
      exit
    end if
    if (animalID == minanimal) then
      backspace(20)
      exit
    end if
    ianimalID = animalID
  end do
  !print *, 'Min, max, imp.a.ID', minanimal, maxanimal, animalID

  !! Go through imputed
  mx=0; my=0; mx2=0; mxy=0; my2=0
  cx(:)=0; cy(:)=0; cx2(:)=0; cxy(:)=0; cy2(:)=0; cNA(:) = 0
  i = 1
  k = 1
  commonrows = 0
  do
    read(20, *, iostat=stat) ianimalID, imputed
    if (stat /= 0) then
      exit
    endif
    !print *, 'Got', ianimalID, 'as imputed...'

    ! Find true genotype
    start = i
    k = 0
    do while (animalIndex(i) /= ianimalID)
      !print *, 'Found', animalIndex(i),'...'
      i = i + 1
      k = k + 1
      if (i > nAnimals) i = 1
      if (i == start) exit
      if (k > nAnimals + 3) then
        !print *, 'Oh oh.'
        exit
      end if
    end do
    if (animalIndex(i) /= ianimalID) then
      i = start
      cycle
    end if
    !print *, 'Settled on', animalIndex(i)

    commonrows = commonrows + 1
    rowcorID(i) = ianimalID
    !print *, 'True:', trueMat(i,:)
    !print *, 'Impu;', imputed(:)

    rx = 0
    ry = 0
    rx2 = 0
    rxy = 0
    ry2 = 0
    rNA = 0

    do j=1,nSnps
      if (imputed(j) == NAval .or. trueMat(i,j) == Naval .or. sds(j) == 0.) then
        rNA = rNA + 1
        cNA(j) = cNA(j) + 1
        cycle
      end if

      t = (trueMat(i,j)-means(j))/sds(j)
      imp = (imputed(j)-means(j))/sds(j)
      t2 = t*t
      imp2 = imp*imp
      tim = t*imp

      cx(j) = cx(j) + t
      mx = mx + t
      rx = rx + t

      cy(j) = cy(j) + imp
      my = my + imp
      ry = ry + imp

      cx2(j) = cx2(j) + t2
      mx2 = mx2 + t2
      rx2 = rx2 + t2

      cy2(j) = cy2(j) + imp2
      my2 = my2 + imp2
      ry2 = ry2 + imp2

      cxy(j) = cxy(j) + tim
      mxy = mxy + tim
      rxy = rxy + tim

    end do ! end of j=1,nSNPs

    nn = j - 1 - rNA
    !print *, animalID, nn, rx, ry, rx2, ry2, rxy, dummy
    rowcors(i) = (nn * rxy - rx * ry) / sqrt( (nn * rx2 - rx**2) * (nn*ry2 - ry**2 ) )
    i = i + 1
  enddo
  close(10)
  close(20)

  nn = commonrows * nSNPs - sum(cNA)
  matcor = (nn * mxy - mx * my) / sqrt( (nn * mx2 - mx**2) * (nn*my2 - my**2) )
  nnn = commonrows  - cNA
  colcors = (nnn * cxy - cx * cy) / sqrt( (nnn * cx2 - cx**2) * (nnn*cy2 - cy**2) )
  if (standardized == 1)  where (sds == 0.) colcors = 1/nan

end subroutine
