!! imputation accuracy

!! Calculates imputation accuracy
!! truefn and imputefn are filepaths of respectively true genotypes and imputed genotypes
!! m, number of SNPs
!! NAval is an integer value of missing genotypes, which are ignored.
!! standardized, boolean, whether to standardize genotypes based on entire true genotypes.
!! rowcor, matcor, colcor, vector of results.
subroutine imp_acc_fast(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, means, sds, &
    usermeans, rowcors, matcor, colcors, rowID, iexids, iexsnps, &
    tol, colcorrect, coltruena, colimpna, colbothna, &
         rowcorrect, rowtruena, rowimpna, rowbothna)
  implicit none
  
  integer, parameter :: r8_kind = selected_real_kind(15, 307) ! double precision, 64-bit like, required for transferring to and fro R.

  !! Arguments
  character(255), intent(in) :: truefn, imputefn
  integer, intent(in) :: nSNPs, NAval, standardized, nAnimals, usermeans
  integer, dimension(nAnimals), intent(in) :: iexids
  integer, dimension(nSNPs), intent(in) :: iexsnps
  integer, dimension(nAnimals), intent(out) :: rowID, rowcorrect, rowtruena, rowimpna, rowbothna
  integer, dimension(nSNPs), intent(out) :: colcorrect, coltruena, colimpna, colbothna
  real(r8_kind), intent(in) :: tol
  real(r8_kind), dimension(nSnps), intent(inout) :: means, sds, colcors
  real(r8_kind), dimension(nAnimals), intent(out) :: rowcors
  real(r8_kind), intent(out) :: matcor

  !! Private variables
  integer :: stat, animalID, i, j, mn
  real(r8_kind) :: tru, imp, nan
  logical, dimension(nAnimals) :: exids
  real(r8_kind), dimension(nSNPs) :: M, S, Mold, Sold
  real, dimension(nSNPs) :: genoin, true, imputed
  !! For running correlation on columns
  integer, dimension(nSNPs) :: cNA, nLines
  real(r8_kind), dimension(nSNPs) :: cmp, cmq, cmt, cmi, csi, cst, csb
  !! For matrix
  real(r8_kind) :: mmp, mmq, mmt, mmi, msi, mst, msb
  !! For row
  integer :: rNA
  real(r8_kind) :: rmp, rmq, rmt, rmi, rsi, rst, rsb

  !! NAN
  !! Quiet NAN, double precision.
  !! Source, Steve Lionel (Intel) https://software.intel.com/en-us/forums/topic/294680
  nan = 0.

  exids = iexids == 1
  

  !! Read through true genotype file and get column-wise mean and variance
  if (standardized == 1 .and. usermeans == 0) then
    open(10, file=truefn, status='OLD')
    nLines(:) = 0
    S(:) = 0
    M(:) = -9
    Mold(:) = 0
    i = 0
    do
      read(10, *, iostat=stat) animalID, genoin
      if (stat /= 0) exit
  
      WHERE (genoin /= NAval)
        nLines = nLines + 1
        WHERE (nLines .eq. 1)
          M = genoin
        ENDWHERE
        WHERE (nLines .gt. 1 .and. M .gt. -9)
          Mold = M
          Sold = S
          M = M + (genoin - Mold)/nLines
          S = S + (genoin - Mold) * (genoin - M)
        ENDWHERE
      ENDWHERE
    enddo
    close(10)
    means = M
    sds = sqrt(S/(nLines - 1))
  elseif (usermeans == 0) then
    means(:) = 0
    sds(:) = 1
  end if

  where(iexsnps == 1) sds = 0.
  
  colcorrect(:) = 0
  coltruena(:) = 0
  colimpna(:) = 0
  colbothna(:) = 0
  rowcorrect(:) = 0
  rowtruena(:) = 0
  rowimpna(:) = 0
  rowbothna(:) = 0

  !! Go through both files
  open(10, file=truefn, status='OLD')
  open(20, file=imputefn, status='OLD')
  i = 0
  cmt(:)=0; cmp(:)=0; cst(:)=0; csi(:)=0; csb(:)=0; cNA(:)=0
  mst=0; msi=0; msb=0; mn = 0
  rowID(:) = 0
  do
    i = i + 1
    read(10, *, iostat=stat) animalID, true
    if (stat /= 0) then
      exit
    endif
    rowID(i) = animalID
    read(20, *, iostat=stat) animalID, imputed
    if (stat /= 0) then
      exit
    endif
    
    if (exids(i) .eqv. .true.) then
      cNA = cNA + 1
      cycle
    endif
    
    rNA = 0
    
    rst = 0
    rsi = 0
    rsb = 0
    
    !print *, 'Animal', animalID

    do j=1,nSnps
      if (imputed(j) == NAval .and. true(j) == NAval) then
        colbothna(j) = colbothna(j) + 1
        rowbothna(i) = rowbothna(i) + 1
      elseif (imputed(j) == NAval) then
        colimpna(j) = colimpna(j) + 1
        rowimpna(i) = rowimpna(i) + 1
      elseif (true(j) == NAval) then
        coltruena(j) = coltruena(j) + 1
        rowtruena(i) = rowtruena(i) + 1
      endif
      if (imputed(j) == NAval .or. true(j) == NAval .or. sds(j) == 0.) then
        rNA = rNA + 1
        cNA(j) = cNA(j) + 1
        cycle
      end if
      mn = mn + 1

      if (abs(true(j) - imputed(j)) .le. tol) then
        colcorrect(j) = colcorrect(j) + 1
        rowcorrect(i) = rowcorrect(i) + 1
      endif

      if (standardized == 1) then
        tru = (true(j)-means(j))/sds(j)
        imp = (imputed(j)-means(j))/sds(j)
      else
        tru = true(j)
        imp = imputed(j)
      endif
      
      ! rowcorrelation
      if (j - rNA .eq. 1) then
        rmt = tru
        rmi = imp
      elseif (j - rNA .gt. 1) then
        rmp = rmt  ! for X, m_n and m_{n_1}
        rmt = rmt + (tru - rmt) / (j - rNA)
        rst = rst + (tru - rmp) * (tru - rmt)
        rmq = rmi  ! for Y, m_n and m_{n_1}
        rmi = rmi + (imp - rmi) / (j - rNA)
        rsi = rsi + (imp - rmq) * (imp - rmi)
        rsb = rsb + (imp - rmi) * (tru - rmp)      
      endif
      
      ! column correlation
      if (i - cNA(j) .eq. 1) then
        cmt(j) = tru
        cmi(j) = imp
      elseif (i - cNA(j) .gt. 1) then
        cmp(j) = cmt(j)  ! for X, m_n and m_{n_1}
        cmt(j) = cmt(j) + (tru - cmt(j)) / (i - cNA(j))
        cst(j) = cst(j) + (tru - cmp(j)) * (tru - cmt(j))
        cmq(j) = cmi(j)  ! for Y, m_n and m_{n_1}
        cmi(j) = cmi(j) + (imp - cmi(j)) / (i - cNA(j))
        csi(j) = csi(j) + (imp - cmq(j)) * (imp - cmi(j))
        csb(j) = csb(j) + (imp - cmi(j)) * (tru - cmp(j))      
      endif 
      
      !print *, j, cmt(j), cmp(j), cst(j)
      
      ! matrix correlation
      if (mn .eq. 1) then
        mmt = tru
        mmi = imp
      elseif (mn .gt. 1) then
        mmp = mmt  ! for X, m_n and m_{n_1}
        mmt = mmt + (tru - mmt) / (mn)
        mst = mst + (tru - mmp) * (tru - mmt)
        mmq = mmi  ! for Y, m_n and m_{n_1}
        mmi = mmi + (imp - mmi) / (mn)
        msi = msi + (imp - mmq) * (imp - mmi)
        msb = msb + (imp - mmi) * (tru - mmp)         
      endif
    enddo
    if (abs(rst) .lt. 1E-24) rst=0
    if (abs(rsi) .lt. 1E-24) rsi=0
    rowcors(i) = rsb / (sqrt(rst) * sqrt(rsi))
  enddo
  close(10)
  close(20)

  ! Evil 32-bit hacks to ensure "no" variance really is perceived as no variance.
  if (abs(mst) .lt. 1E-24) mst=0
  if (abs(msi) .lt. 1E-24) msi=0
  matcor = msb / (sqrt(mst) * sqrt(msi))

  where(abs(cst) .lt. 1E-24) cst=0
  where(abs(csi) .lt. 1E-24) csi=0
  colcors = csb / (sqrt(cst) * sqrt(csi))
  if (standardized == 1) where (sds == 0) colcors = 1/nan

  !flush(6)

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
subroutine imp_acc(truefn, imputefn, nSNPs, nAnimals, NAval, standardized, means, sds, &
  usermeans, rowcors, matcor, colcors, rowID, iexids, iexsnps, &
    tol, colcorrect, coltruena, colimpna, colbothna, &
         rowcorrect, rowtruena, rowimpna, rowbothna)
  implicit none
  
  integer, parameter :: r8_kind = selected_real_kind(15, 150) ! double precision, 64-bit like, required for transferring to and fro R.

  character(255), intent(in) :: truefn, imputefn
  integer, intent(in) :: nSNPs, NAval, standardized, nAnimals, usermeans
  integer, dimension(nAnimals), intent(in) :: iexids
  integer, dimension(nSNPs), intent(in) :: iexsnps
  integer, dimension(nAnimals), intent(out) :: rowID, rowcorrect, rowtruena, rowimpna, rowbothna
  integer, dimension(nSNPs), intent(out) :: colcorrect, coltruena, colimpna, colbothna
  real(r8_kind), intent(in) :: tol
  real(r8_kind), dimension(nSnps), intent(inout) :: means, sds, colcors
  real(r8_kind), dimension(nAnimals), intent(out) :: rowcors
  real(r8_kind), intent(out) :: matcor
  !! Private variables
  logical, dimension(nAnimals) :: foundID
  integer :: stat, start, animalID, commonrows, i, j, k, l, maxanimal, minanimal, ianimalID
  real(r8_kind) :: tru, imp, nan
  real(r8_kind), dimension(nSNPs) :: M, S, Mold, Sold
  real(r8_kind), dimension(nSNPs) :: genoin, imputed
  real(r8_kind), dimension(nAnimals, nSnps) :: trueMat
  !! For running correlation on columns
  integer, dimension(nSNPs) :: cNA, nLines
  integer, dimension(nAnimals) :: animalIndex
  real(r8_kind), dimension(nSNPs) :: cmp, cmq, cmt, cmi, csi, cst, csb
  !! For matrix
  integer :: mn
  real(r8_kind) :: mmp, mmq, mmt, mmi, msi, mst, msb
  !! For row
  integer :: rNA
  real(r8_kind) :: rmp, rmq, rmt, rmi, rsi, rst, rsb  

  nan = 0.0

  rowID(:) = 0
  rowcors(:) = 0.0
  colcors(:) = 0.0
  foundID(:) = .false.
  
  colcorrect(:) = 0
  coltruena(:) = 0
  colimpna(:) = 0
  colbothna(:) = 0
  rowcorrect(:) = 0
  rowtruena(:) = 0
  rowimpna(:) = 0
  rowbothna(:) = 0
  
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

  !! Calculate scaling parameters (mean and var)
  if (standardized == 1 .and. usermeans == 0) then
    nLines(:) = 0
    S(:) = 0
    M(:) = -9
    Mold(:) = 0
    i = 0
    do i=1,nAnimals
      genoin = trueMat(i,:)

      WHERE (genoin /= NAval)
        nLines = nLines + 1
        WHERE (nLines .eq. 1)
          M = genoin
        ENDWHERE
        WHERE (nLines .gt. 1 .and. M .gt. -9)
          Mold = M
          Sold = S
          M = M + (genoin - Mold)/nLines
          S = S + (genoin - Mold) * (genoin - M)
        ENDWHERE
      ENDWHERE
    enddo
    close(10)
    means = M
    sds = sqrt(S/(nLines - 1))
  elseif (usermeans == 0) then
    means(:) = 0
    sds(:) = 1
  end if

  where(iexsnps == 1) sds = 0.

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
  cst(:)=0; csi(:)=0; csb(:)=0; cNA(:)=0
  mst=0; msi=0; msb=0; mn = 0
  i = 1
  k = 1
  l = 0
  commonrows = 0
  do
    l = l + 1
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

    
    rowID(i) = ianimalID
    foundID(i) = .true.

    if (iexids(i) == 1) cycle

    commonrows = commonrows + 1
    
    rNA = 0
    rst = 0
    rsi = 0
    rsb = 0

    do j=1,nSnps
      if (imputed(j) == NAval .or. trueMat(i,j) == Naval .or. sds(j) == 0.) then
        rNA = rNA + 1
        cNA(j) = cNA(j) + 1
        cycle
      end if
      mn = mn + 1

      if (standardized == 1) then
        tru = (trueMat(i,j)-means(j))/sds(j)
        imp = (imputed(j)-means(j))/sds(j)
      else
        tru = trueMat(i,j)
        imp = imputed(j)
      endif

      ! rowcorrelation
      if (j - rNA .eq. 1) then
        rmt = tru
        rmi = imp
      elseif (j - rNA .gt. 1) then
        rmp = rmt  ! for X, m_n and m_{n_1}
        rmt = rmt + (tru - rmt) / (j - rNA)
        rst = rst + (tru - rmp) * (tru - rmt)
        rmq = rmi  ! for Y, m_n and m_{n_1}
        rmi = rmi + (imp - rmi) / (j - rNA)
        rsi = rsi + (imp - rmq) * (imp - rmi)
        rsb = rsb + (imp - rmi) * (tru - rmp)      
      endif
      
      ! column correlation
      if (commonrows - cNA(j) .eq. 1) then
        cmt(j) = tru
        cmi(j) = imp
      elseif (commonrows - cNA(j) .gt. 1) then
        cmp(j) = cmt(j)  ! for X, m_n and m_{n_1}
        cmt(j) = cmt(j) + (tru - cmt(j)) / (commonrows - cNA(j))
        cst(j) = cst(j) + (tru - cmp(j)) * (tru - cmt(j))
        cmq(j) = cmi(j)  ! for Y, m_n and m_{n_1}
        cmi(j) = cmi(j) + (imp - cmi(j)) / (commonrows - cNA(j))
        csi(j) = csi(j) + (imp - cmq(j)) * (imp - cmi(j))
        csb(j) = csb(j) + (imp - cmi(j)) * (tru - cmp(j))      
      endif        
      
      ! matrix correlation
      if (mn .eq. 1) then
        mmt = tru
        mmi = imp
      elseif (mn .gt. 1) then
        mmp = mmt  ! for X, m_n and m_{n_1}
        mmt = mmt + (tru - mmt) / (mn)
        mst = mst + (tru - mmp) * (tru - mmt)
        mmq = mmi  ! for Y, m_n and m_{n_1}
        mmi = mmi + (imp - mmi) / (mn)
        msi = msi + (imp - mmq) * (imp - mmi)
        msb = msb + (imp - mmi) * (tru - mmp)         
      endif
      
    end do ! end of j=1,nSNPs
    if (abs(rst) .lt. 1E-24) rst=0
    if (abs(rsi) .lt. 1E-24) rsi=0
    rowcors(i) = rsb / (sqrt(rst) * sqrt(rsi))
    i = i + 1
  enddo
  close(10)
  close(20)

  ! Evil 32-bit hacks to ensure "no" variance really is perceived as no variance.
  if (abs(mst) .lt. 1E-24) mst=0
  if (abs(msi) .lt. 1E-24) msi=0
  matcor = msb / (sqrt(mst) * sqrt(msi))
  
  where(abs(cst) .lt. 1E-24) cst=0
  where(abs(csi) .lt. 1E-24) csi=0
  colcors = csb / (sqrt(cst) * sqrt(csi))
  
  if (standardized == 1)  where (sds == 0.) colcors = 1/nan
  where (foundID .eqv. .false.)
    rowcors = 1/nan
    rowID = 1/nan
  end where

end subroutine
