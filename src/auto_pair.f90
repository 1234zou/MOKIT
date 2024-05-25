! written by jxzou at 20180416
! Automatically pair the active orbitals according to distances between any one active occupied orbital
!  and one active unoccupied orbital
! Note the order input orbitals must be Gaussian type: core, pair(occ), open, pair(vir), virtual

! modified by jxzou at 20180524: add subroutine pair_by_tdm (pair by transition dipole moments)
! modified by jxzou at 20180611: cycle in pairing
! modified by jxzou at 20180705: implement cases where occpuied LMO < unoccpuied LMO

! paired by minimizing the sum of orbitals distances (<i|r|i>-<j|r|j>)**2
subroutine pair_by_dis(ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif, coeff, mo_dipole, new_coeff)
 implicit none
 integer ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif
!f2py intent(in) :: ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif
 real(kind=8) coeff(nbf,nif), mo_dipole(3,nif,nif), new_coeff(nbf,nif)
!f2py intent(in,copy) :: coeff
!f2py depend(nbf,nif) :: coeff
!f2py intent(in) :: mo_dipole
!f2py depend(nif) :: mo_dipole
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: new_coeff

 if(nalpha-nopen-ncore>=npair .and. nvir_lmo==npair) then
  call pair_by_dis1(ncore, npair, nopen, nalpha, nbf, nif, coeff, mo_dipole)
  new_coeff = coeff
 else if(nalpha-nopen-ncore==npair .and. nvir_lmo>npair) then
  call pair_by_dis2(ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif, coeff, mo_dipole)
  new_coeff = coeff
 else
  write(6,'(A,I4)') 'ERROR in subroutine pair_by_dis: nalpha-nopen-ncore=', nalpha-nopen-ncore
  write(6,'(A,I4,A,I4)') 'And nvir_lmo=', nvir_lmo, ', But npair=', npair
  stop
 end if
end subroutine pair_by_dis

! for case where occpuied LMO >= unoccpuied LMO
subroutine pair_by_dis1(ncore, npair, nopen, nalpha, nbf, nif, coeff, mo_dipole)
 implicit none
 integer i, j, k, nocc_lmo, nmid
 integer tmp_idx(1), tmp_idx2(2)
 integer, intent(in) :: ncore, npair, nopen, nalpha, nbf, nif
 integer, allocatable :: pair_idx(:), opt_pair_idx(:)
 real(kind=8), intent(inout) :: coeff(nbf,nif)
 real(kind=8), intent(in) :: mo_dipole(3,nif,nif)
 real(kind=8), allocatable :: dis(:,:), dis1(:,:)
!                               occ vir
 real(kind=8) fBoys, tempv, sum_dis, min_sum
 real(kind=8), parameter :: rdiff = 1d-4
 real(kind=8) temp_dipole1(3), temp_dipole2(3)
 real(kind=8), allocatable :: coeff1(:,:)
 logical, allocatable :: used(:)

 nmid = nalpha - nopen - npair - ncore
 nocc_lmo = npair + nmid
 write(6,'(/,A)') 'Information of subroutine pair_by_dis1:'
 write(6,'(A9,I4,A)') 'There are', ncore,    ' core orbitals.'
 write(6,'(A9,I4,A)') 'There are', nocc_lmo, ' occupied LMOs.'
 write(6,'(A9,I4,A)') 'There are', nopen,    ' singly occupied MOs.'
 write(6,'(A9,I4,A)') 'There are', npair,    ' unoccupied LMOs.'
 write(6,'(A9,I4,A)') 'There are', nalpha,   ' alpha MOs.'

 ! calculate the modified Boys values in occ and vir LMO subspaces, respectively
 fBoys = 0d0
 j = ncore + nocc_lmo
 do i = ncore+1, j, 1
  temp_dipole1 =  mo_dipole(1:3,i,i)
  fBoys = fBoys + DOT_PRODUCT(temp_dipole1,temp_dipole1)
 end do
 fBoys = DSQRT(fBoys/DBLE(nocc_lmo))
 write(6,'(A,F12.6)') 'In occ LMO subspace, Modified f(Boys)=', fBoys
 fBoys = 0d0
 j = nalpha + npair
 do i = nalpha+1, j, 1
  temp_dipole1 =  mo_dipole(1:3,i,i)
  fBoys = fBoys + DOT_PRODUCT(temp_dipole1,temp_dipole1)
 end do
 fBoys = DSQRT(fBoys/DBLE(npair))
 write(6,'(A,F12.6)') 'In vir LMO subspace, Modified f(Boys)=', fBoys

 ! calculate the distances between any one active occ orbital and any one active vir orbital
 allocate(dis(nocc_lmo,npair), dis1(nocc_lmo,npair), used(nocc_lmo))
 dis = 0d0
 dis1 = 0d0
 used = .false.
 do i = nalpha+1, nalpha+npair, 1
  k = i - nalpha
  temp_dipole1 = mo_dipole(1:3,i,i)
  do j = ncore+1, ncore+nocc_lmo, 1
   temp_dipole2 = temp_dipole1 - mo_dipole(1:3,j,j)
   dis(j-ncore,k) = DOT_PRODUCT(temp_dipole2,temp_dipole2)
  end do
 end do
 tmp_idx2 = MINLOC(dis)
 write(6,'(A,F11.5)') 'Min pair distance**2=', dis(tmp_idx2(1),tmp_idx2(2))
 write(6,'(2I5)') tmp_idx2(1)+ncore, tmp_idx2(2)+nalpha

 ! find the closest vir orbital of an occ orbital
 allocate(pair_idx(npair), opt_pair_idx(npair))
 pair_idx = 0
 opt_pair_idx = 0
 sum_dis = 0d0
 ! the unoccupied LMOs are unchanged; pick proper number of occupied LMOs from occ space
 tempv = MAXVAL(dis) + 1d0
 min_sum = SUM(dis)

 do j = 1, npair, 1
  used = .false.
  tmp_idx = 0
  do i = j, npair, 1
   dis1 = dis
   do while(.true.)
    tmp_idx = MINLOC(dis1(:,i))
    if(.not. used(tmp_idx(1))) exit
    dis1(tmp_idx(1),i) = tempv
   end do
   pair_idx(i) = tmp_idx(1)
   used(tmp_idx(1)) = .true.
  end do
  do i = 1, j-1, 1
   dis1 = dis
   do while(.true.)
    tmp_idx = MINLOC(dis1(:,i))
    if(.not. used(tmp_idx(1))) exit
    dis1(tmp_idx(1),i) = tempv
   end do
   pair_idx(i) = tmp_idx(1)
   used(tmp_idx(1)) = .true.
  end do
  sum_dis = 0d0
  do i = 1, npair, 1
   sum_dis = sum_dis + dis(pair_idx(i),i)
  end do
  write(6,'(A8,F16.5)') 'sum_dis=', sum_dis
  if(sum_dis-min_sum < rdiff) then
   min_sum = sum_dis
   opt_pair_idx = pair_idx
  end if
 end do 

 write(6,'(A)') 'Final pairs:'
 do i = 1, npair, 1
  j = opt_pair_idx(i)
  write(6,'(3I4,F11.5)') i, i+nalpha, j+ncore, dis(j,i)
 end do
 write(6,'(A8,F16.5)') 'min_sum=', min_sum
 write(6,'(A)') 'Information of subroutine pair_by_dis2 printing done.'

 ! put new MO into the array coeff
 allocate(coeff1(nbf,nif))
 coeff1 = 0d0
 coeff1(:,1:ncore) = coeff(:,1:ncore)
 coeff1(:,nalpha-nopen+1:nif) = coeff(:,nalpha-nopen+1:nif)
 forall(i = 1:npair)
  coeff1(:,ncore+nmid+i) = coeff(:,opt_pair_idx(npair-i+1)+ncore)
 end forall
 j = 0
 do i = 1, nocc_lmo, 1
  if( ALL(opt_pair_idx(:)/=i) ) then
   j = j + 1
   coeff1(:,ncore+j) = coeff(:,ncore+i)
  end if
 end do
 coeff = coeff1
 deallocate(pair_idx, opt_pair_idx, coeff1, dis, dis1, used)
end subroutine pair_by_dis1

! for case where occpuied LMO < unoccpuied LMO
subroutine pair_by_dis2(ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif, coeff, mo_dipole)
 implicit none
 integer i, j, k
 integer tmp_idx(1), tmp_idx2(2)
 integer, intent(in) :: ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif
 integer, allocatable :: pair_idx(:), opt_pair_idx(:)
 real(kind=8), intent(inout) :: coeff(nbf,nif)
 real(kind=8), intent(in) :: mo_dipole(3,nif,nif)
 real(kind=8), allocatable :: dis(:,:), dis1(:,:)
!                               vir occ
 real(kind=8) fBoys, tempv, sum_dis, min_sum
 real(kind=8), parameter :: rdiff = 1d-4
 real(kind=8) temp_dipole1(3), temp_dipole2(3)
 real(kind=8), allocatable :: coeff1(:,:)
 logical, allocatable :: used(:)

 write(6,'(/,A)') 'Information of subroutine pair_by_dis2:'
 write(6,'(A9,I4,A)') 'There are', ncore,    ' core orbitals.'
 write(6,'(A9,I4,A)') 'There are', npair,    ' occupied LMOs.'
 write(6,'(A9,I4,A)') 'There are', nopen,    ' singly occupied MOs.'
 write(6,'(A9,I4,A)') 'There are', nvir_lmo, ' unoccupied LMOs.'
 write(6,'(A9,I4,A)') 'There are', nalpha,   ' alpha MOs.'

 ! calculate the modified Boys values in occ and vir LMO subspaces, respectively
 fBoys = 0d0
 j = ncore + npair
 do i = ncore+1, j, 1
  temp_dipole1 =  mo_dipole(1:3,i,i)
  fBoys = fBoys + DOT_PRODUCT(temp_dipole1,temp_dipole1)
 end do
 fBoys = DSQRT(fBoys/DBLE(npair))
 write(6,'(A,F12.6)') 'In occ LMO subspace, Modified f(Boys)=', fBoys
 fBoys = 0d0
 j = nalpha + nvir_lmo
 do i = nalpha+1, j, 1
  temp_dipole1 =  mo_dipole(1:3,i,i)
  fBoys = fBoys + DOT_PRODUCT(temp_dipole1,temp_dipole1)
 end do
 fBoys = DSQRT(fBoys/DBLE(nvir_lmo))
 write(6,'(A,F12.6)') 'In vir LMO subspace, Modified f(Boys)=', fBoys

 ! calculate the distances between any one active occ orbital and any one active vir orbital
 allocate(dis(nvir_lmo,npair), dis1(nvir_lmo,npair), used(nvir_lmo))
 dis = 0d0
 dis1 = 0d0
 used = .false.
 do i = ncore+1, ncore+npair, 1
  k = i - ncore
  temp_dipole1 = mo_dipole(1:3,i,i)
  do j = nalpha+1, nalpha+nvir_lmo, 1
   temp_dipole2 = temp_dipole1 - mo_dipole(1:3,j,j)
   dis(j-nalpha,k) = DOT_PRODUCT(temp_dipole2,temp_dipole2)
  end do
 end do
 tmp_idx2 = MINLOC(dis)
 write(6,'(A,F11.5)') 'Min pair distance**2=', dis(tmp_idx2(1),tmp_idx2(2))
 write(6,'(2I5)') tmp_idx2(1)+nalpha, tmp_idx2(2)+ncore

 ! find the closest vir orbital of an occ orbital
 allocate(pair_idx(npair), opt_pair_idx(npair))
 pair_idx = 0
 opt_pair_idx = 0
 sum_dis = 0d0
 ! the occupied LMOs are unchanged; pick proper number of unoccupied LMOs from occ space
 tempv = MAXVAL(dis) + 1d0
 min_sum = SUM(dis)

 do j = 1, npair, 1
  used = .false.
  tmp_idx = 0
  do i = j, npair, 1
   dis1 = dis
   do while(.true.)
    tmp_idx = MINLOC(dis1(:,i))
    if(.not. used(tmp_idx(1))) exit
    dis1(tmp_idx(1),i) = tempv
   end do
   pair_idx(i) = tmp_idx(1)
   used(tmp_idx(1)) = .true.
  end do
  do i = 1, j-1, 1
   dis1 = dis
   do while(.true.)
    tmp_idx = MINLOC(dis1(:,i))
    if(.not. used(tmp_idx(1))) exit
    dis1(tmp_idx(1),i) = tempv
   end do
   pair_idx(i) = tmp_idx(1)
   used(tmp_idx(1)) = .true.
  end do
  sum_dis = 0d0
  do i = 1, npair, 1
   sum_dis = sum_dis + dis(pair_idx(i),i)
  end do
  write(6,'(A8,F16.5)') 'sum_dis=', sum_dis
  if(sum_dis-min_sum < rdiff) then
   min_sum = sum_dis
   opt_pair_idx = pair_idx
  end if
 end do 

 write(6,'(A)') 'Final pairs:'
 do i = 1, npair, 1
  j = opt_pair_idx(i)
  write(6,'(3I4,F11.5)') i, i+ncore, j+nalpha, dis(j,i)
 end do
 write(6,'(A8,F16.5)') 'min_sum=', min_sum
 write(6,'(A)') 'Information of subroutine pair_by_dis2 printing done.'

 ! put new MO into the array coeff
 allocate(coeff1(nbf,nif))
 coeff1 = 0d0
 coeff1(:,1:nalpha) = coeff(:,1:nalpha)
 coeff1(:,nalpha+nvir_lmo+1:nif) = coeff(:,nalpha+nvir_lmo+1:nif)
 forall(i = 1:npair)
  coeff1(:,i+nalpha) = coeff(:,opt_pair_idx(npair-i+1)+nalpha)
 end forall
 j = 0
 do i = 1, nvir_lmo, 1
  if( ALL(opt_pair_idx(:)/=i) ) then
   j = j + 1
   coeff1(:,nalpha+npair+j) = coeff(:,nalpha+i)
  end if
 end do
 coeff = coeff1
 deallocate(pair_idx, opt_pair_idx, coeff1, dis, dis1, used)
end subroutine pair_by_dis2

! paired by maximizing the sum of orbital transition dipoles <i|r|j>**2
subroutine pair_by_tdm(ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif, coeff, mo_dipole, new_coeff)
 implicit none
 integer :: i
 integer :: ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif
!f2py intent(in) :: ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif
 real(kind=8) :: coeff(nbf,nif), mo_dipole(3,nif,nif), new_coeff(nbf,nif)
!f2py intent(in,copy) :: coeff
!f2py intent(in) :: mo_dipole
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: coeff, new_coeff
!f2py depend(nif) :: mo_dipole

 i = nalpha - nopen - ncore

 if(i>=npair .and. nvir_lmo==npair) then
  call pair_by_tdm1(ncore, npair, nopen, nalpha, nbf, nif, coeff, mo_dipole)
  new_coeff = coeff
 else if(i==npair .and. nvir_lmo>npair) then
  call pair_by_tdm2(ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif, coeff, mo_dipole)
  new_coeff = coeff
 else
  write(6,'(/,A,I0)') 'ERROR in subroutine pair_by_tdm: nalpha-nopen-ncore = ', i
  write(6,'(2(A,I0))') 'And nvir_lmo = ', nvir_lmo, ', But npair = ', npair
  stop
 end if
end subroutine pair_by_tdm

! for case where occpuied LMO >= unoccpuied LMO
subroutine pair_by_tdm1(ncore, npair, nopen, nalpha, nbf, nif, coeff, mo_dipole)
 implicit none
 integer i, j, k, nocc_lmo, nmid
 integer tmp_idx(1), tmp_idx2(2)
 integer, intent(in) :: ncore, npair, nopen, nalpha, nbf, nif
 integer, allocatable :: pair_idx(:), opt_pair_idx(:)
 real(kind=8), intent(inout) :: coeff(nbf,nif)
 real(kind=8), intent(in) :: mo_dipole(3,nif,nif)
 real(kind=8) fBoys, sum_tdm, max_sum, min_tdm
 real(kind=8) tempv, tempv1, temp_dipole(3)
 real(kind=8), parameter :: rdiff1 = 1d-4, rdiff2 = 1d-2
 real(kind=8), allocatable :: coeff1(:,:)
 real(kind=8), allocatable :: tdm(:,:), tdm1(:,:)
!                               occ vir
 logical, allocatable :: used(:)

 nmid = nalpha - nopen - npair - ncore
 nocc_lmo = npair + nmid
 write(6,'(/,A)') 'Information of subroutine pair_by_tdm1:'
 write(6,'(A9,I4,A)') 'There are', ncore,    ' core orbitals.'
 write(6,'(A9,I4,A)') 'There are', nocc_lmo, ' occupied LMOs.'
 write(6,'(A9,I4,A)') 'There are', nopen,    ' singly occupied MOs.'
 write(6,'(A9,I4,A)') 'There are', nalpha,   ' alpha MOs.'
 write(6,'(A9,I4,A)') 'There are', npair,    ' unoccupied LMOs.'

 ! calculate the modified Boys values in occ and vir LMO subspaces, respectively
 fBoys = 0d0
 j = ncore + nocc_lmo
 do i = ncore+1, j, 1
  temp_dipole =  mo_dipole(1:3,i,i)
  fBoys = fBoys + DOT_PRODUCT(temp_dipole,temp_dipole)
 end do
 fBoys = DSQRT(fBoys/DBLE(nocc_lmo))
 write(6,'(A,F12.6)') 'In occ LMO subspace, Modified f(Boys)=', fBoys
 fBoys = 0d0
 j = nalpha + npair
 do i = nalpha+1, j, 1
  temp_dipole =  mo_dipole(1:3,i,i)
  fBoys = fBoys + DOT_PRODUCT(temp_dipole,temp_dipole)
 end do
 fBoys = DSQRT(fBoys/DBLE(npair))
 write(6,'(A,F12.6)') 'In vir LMO subspace, Modified f(Boys)=', fBoys

 ! calculate the transition dipoles between any one active occ orbital and any one active vir orbital
 allocate(tdm(nocc_lmo,npair), tdm1(nocc_lmo,npair), used(nocc_lmo))
 tdm = 0d0
 tdm1 = 0d0
 used = .false.
 do i = nalpha+1, nalpha+npair, 1
  k = i - nalpha
  do j = ncore+1, ncore+nocc_lmo, 1
   temp_dipole =  mo_dipole(1:3,j,i)
   tdm(j-ncore,k) = DOT_PRODUCT(temp_dipole,temp_dipole)
  end do
 end do
 tmp_idx2 = MAXLOC(tdm)
 write(6,'(A,F11.5)') 'Max transition dipole**2=', tdm(tmp_idx2(1),tmp_idx2(2))
 write(6,'(2I5)') tmp_idx2(1)+ncore, tmp_idx2(2)+nalpha

 ! find the closest vir orbital of an occ orbital
 allocate(pair_idx(npair), opt_pair_idx(npair))
 pair_idx = 0
 opt_pair_idx = 0
 sum_tdm = 0d0
 ! the unoccupied LMOs are unchanged; pick proper number of occupied LMOs from occ space
 tempv = MINVAL(tdm) - 1d0
 max_sum = 0d0

 do j = 1, npair, 1
  used = .false.
  tmp_idx = 0

  do i = j, npair, 1
   tdm1 = tdm
   do while(.true.)
    tmp_idx = MAXLOC(tdm1(:,i))
    if(.not. used(tmp_idx(1))) exit
    tdm1(tmp_idx(1),i) = tempv
   end do ! for while
   pair_idx(i) = tmp_idx(1)
   used(tmp_idx(1)) = .true.
  end do ! for i

  do i = 1, j-1, 1
   tdm1 = tdm
   do while(.true.)
    tmp_idx = MAXLOC(tdm1(:,i))
    if(.not. used(tmp_idx(1))) exit
    tdm1(tmp_idx(1),i) = tempv
   end do ! for while
   pair_idx(i) = tmp_idx(1)
   used(tmp_idx(1)) = .true.
  end do ! for i

  sum_tdm = 0d0
  min_tdm = tdm(pair_idx(1),1) + 1d0
  do i = 1, npair, 1
   tempv1 = tdm(pair_idx(i),i)
   sum_tdm = sum_tdm + tempv1
   if(tempv1 < min_tdm) min_tdm = tempv1
  end do ! for i

  write(6,'(A8,F16.5)') 'sum_tdm=', sum_tdm
  if(j == 1) then
   max_sum = sum_tdm
   opt_pair_idx = pair_idx
  else
   if(sum_tdm-max_sum>rdiff1 .and. min_tdm>rdiff2) then
    max_sum = sum_tdm
    opt_pair_idx = pair_idx
   end if
  end if
 end do ! for j

 write(6,'(A)') 'Final pairs:'
 do i = 1, npair, 1
  j = opt_pair_idx(i)
  write(6,'(3I4,F11.5)') i, i+nalpha, j+ncore, tdm(j,i)
 end do
 write(6,'(A8,F16.5)') 'max_sum=', max_sum
 write(6,'(A)') 'Information of subroutine pair_by_tdm1 printing done.'

 ! put new MO into the array coeff
 allocate(coeff1(nbf,nif))
 coeff1 = 0d0
 coeff1(:,1:ncore) = coeff(:,1:ncore)
 coeff1(:,nalpha-nopen+1:nif) = coeff(:,nalpha-nopen+1:nif)
 forall(i = 1:npair)
  coeff1(:,ncore+nmid+i) = coeff(:,opt_pair_idx(npair-i+1)+ncore)
 end forall
 j = 0
 do i = 1, nocc_lmo, 1
  if( ALL(opt_pair_idx(:)/=i) ) then
   j = j + 1
   coeff1(:,ncore+j) = coeff(:,ncore+i)
  end if
 end do
 coeff = coeff1
 deallocate(pair_idx, opt_pair_idx, coeff1, tdm, tdm1, used)
end subroutine pair_by_tdm1

! for case where occpuied LMO < unoccpuied LMO
subroutine pair_by_tdm2(ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif, coeff, mo_dipole)
 implicit none
 integer i, j, k
 integer tmp_idx(1), tmp_idx2(2)
 integer, intent(in) :: ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif
 integer, allocatable :: pair_idx(:), opt_pair_idx(:)
 real(kind=8), intent(inout) :: coeff(nbf,nif)
 real(kind=8), intent(in) :: mo_dipole(3,nif,nif)
 real(kind=8) fBoys, sum_tdm, max_sum, min_tdm
 real(kind=8) tempv, tempv1, temp_dipole(3)
 real(kind=8), parameter :: rdiff1 = 1d-4, rdiff2 = 1d-2
 real(kind=8), allocatable :: coeff1(:,:)
 real(kind=8), allocatable :: tdm(:,:), tdm1(:,:)
!                               vir occ
 logical, allocatable :: used(:)

 write(6,'(/,A)') 'Information of subroutine pair_by_tdm2:'
 write(6,'(A9,I4,A)') 'There are', ncore,    ' core orbitals.'
 write(6,'(A9,I4,A)') 'There are', npair,    ' occupied LMOs.'
 write(6,'(A9,I4,A)') 'There are', nopen,    ' singly occupied MOs.'
 write(6,'(A9,I4,A)') 'There are', nalpha,   ' alpha MOs.'
 write(6,'(A9,I4,A)') 'There are', nvir_lmo, ' unoccupied LMOs.'

 ! calculate the modified Boys values in occ and vir LMO subspaces, respectively
 fBoys = 0d0
 j = ncore + npair
 do i = ncore+1, j, 1
  temp_dipole =  mo_dipole(1:3,i,i)
  fBoys = fBoys + DOT_PRODUCT(temp_dipole,temp_dipole)
 end do
 fBoys = DSQRT(fBoys/DBLE(npair))
 write(6,'(A,F12.6)') 'In occ LMO subspace, Modified f(Boys)=', fBoys
 fBoys = 0d0
 j = nalpha + nvir_lmo
 do i = nalpha+1, j, 1
  temp_dipole =  mo_dipole(1:3,i,i)
  fBoys = fBoys + DOT_PRODUCT(temp_dipole,temp_dipole)
 end do
 fBoys = DSQRT(fBoys/DBLE(nvir_lmo))
 write(6,'(A,F12.6)') 'In vir LMO subspace, Modified f(Boys)=', fBoys

 ! calculate the transition dipoles between any one active occ orbital and any one active vir orbital
 allocate(tdm(nvir_lmo,npair), tdm1(nvir_lmo,npair), used(nvir_lmo))
 tdm = 0d0
 tdm1 = 0d0
 used = .false.
 do i = ncore+1, ncore+npair, 1
  k = i - ncore
  do j = nalpha+1, nalpha+nvir_lmo, 1
   temp_dipole =  mo_dipole(1:3,j,i)
   tdm(j-nalpha,k) = DOT_PRODUCT(temp_dipole,temp_dipole)
  end do
 end do
 tmp_idx2 = MAXLOC(tdm)
 write(6,'(A,F11.5)') 'Max transition dipole**2=', tdm(tmp_idx2(1),tmp_idx2(2))
 write(6,'(2I5)') tmp_idx2(1)+nalpha, tmp_idx2(2)+ncore

 ! find the closest vir orbital of an occ orbital
 allocate(pair_idx(npair), opt_pair_idx(npair))
 pair_idx = 0
 opt_pair_idx = 0
 sum_tdm = 0d0
 ! the occupied LMOs are unchanged; pick proper number of unoccupied LMOs from occ space
 tempv = MINVAL(tdm) - 1d0
 max_sum = 0d0

 do j = 1, npair, 1
  used = .false.
  tmp_idx = 0
  do i = j, npair, 1
   tdm1 = tdm
   do while(.true.)
    tmp_idx = MAXLOC(tdm1(:,i))
    if(.not. used(tmp_idx(1))) exit
    tdm1(tmp_idx(1),i) = tempv
   end do
   pair_idx(i) = tmp_idx(1)
   used(tmp_idx(1)) = .true.
  end do
  do i = 1, j-1, 1
   tdm1 = tdm
   do while(.true.)
    tmp_idx = MAXLOC(tdm1(:,i))
    if(.not. used(tmp_idx(1))) exit
    tdm1(tmp_idx(1),i) = tempv
   end do
   pair_idx(i) = tmp_idx(1)
   used(tmp_idx(1)) = .true.
  end do
  sum_tdm = 0d0
  min_tdm = tdm(pair_idx(1),1) + 1d0
  do i = 1, npair, 1
   tempv1 = tdm(pair_idx(i),i)
   sum_tdm = sum_tdm + tempv1
   if(tempv1 < min_tdm) min_tdm = tempv1
  end do
  write(6,'(A8,F16.5)') 'sum_tdm=', sum_tdm
  if(j == 1) then
   max_sum = sum_tdm
   opt_pair_idx = pair_idx
  else
   if(sum_tdm-max_sum>rdiff1 .and. min_tdm>rdiff2) then
    max_sum = sum_tdm
    opt_pair_idx = pair_idx
   end if
  end if
 end do 

 write(6,'(A)') 'Final pairs:'
 do i = 1, npair, 1
  j = opt_pair_idx(i)
  write(6,'(3I4,F11.5)') i, i+ncore, j+nalpha, tdm(j,i)
 end do
 write(6,'(A8,F16.5)') 'max_sum=', max_sum
 write(6,'(A)') 'Information of subroutine pair_by_tdm2 printing done.'

 ! put new MO into the array coeff
 allocate(coeff1(nbf,nif))
 coeff1 = 0d0
 coeff1(:,1:nalpha) = coeff(:,1:nalpha)
 coeff1(:,nalpha+nvir_lmo+1:nif) = coeff(:,nalpha+nvir_lmo+1:nif)
 forall(i = 1:npair)
  coeff1(:,i+nalpha) = coeff(:,opt_pair_idx(npair-i+1)+nalpha)
 end forall
 j = 0
 do i = 1, nvir_lmo, 1
  if( ALL(opt_pair_idx(:)/=i) ) then
   j = j + 1
   coeff1(:,nalpha+npair+j) = coeff(:,nalpha+i)
  end if
 end do
 coeff = coeff1
 deallocate(pair_idx, opt_pair_idx, coeff1, tdm, tdm1, used)
end subroutine pair_by_tdm2

