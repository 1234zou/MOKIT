! written by jxzou at 20201209: move some subroutines from bas_gms2molcas.f90 to this file
! updated by jxzou at 20210222: merge subroutine read_prim_gau1 and read_prim_gau2

module pg
 implicit none
 integer :: natom     ! the number of atoms
 integer :: highest   ! highest angular momentum
 integer, allocatable :: nuc(:), ntimes(:)
 real(kind=8), allocatable :: coor(:,:)    ! Cartesian coordinates
 character(len=2), allocatable :: elem(:)  ! elements

 ! 'L' will be divided into two parts: 'S' and 'P'
 type :: primitive_gaussian
  character(len=1) :: stype = ' ' ! 'S','P','D','F','G','H','I'
  integer :: nline = 0
  integer :: ncol  = 0
  real(kind=8), allocatable :: coeff(:,:)
 end type primitive_gaussian

 ! 7 for 'S', 'P', 'D', 'F', 'G', 'H', 'I'
 !        1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7
 type(primitive_gaussian) :: prim_gau(7)

 ! --- below are arrays for ECP/PP ---
 type :: ecp_potential
  integer :: n = 0   ! size for col1, 2, and 3
  real(kind=8), allocatable :: col1(:)
  integer, allocatable :: col2(:)
  real(kind=8), allocatable :: col3(:)
 end type ecp_potential

 type :: ecp4atom
  logical :: ecp = .false. ! whether this atom has ECP/PP
  integer :: core_e  = 0   ! number of core electrons
  integer :: highest = 0   ! highest angular momentum
  type(ecp_potential), allocatable :: potential(:) ! size highest+1
 end type ecp4atom

 type(ecp4atom), allocatable :: all_ecp(:) ! size natom

 logical :: ecp_exist = .false.
 logical, allocatable :: ghost(:) ! size natom
end module pg

! Shared variables and arrays for serial and parallel Jacobi 2-by-2 rotations.
! The derivative type in this module cannot be recognized by f2py. So I have to
! put them here.
module lo_info
 implicit none
 integer :: np, npair, eff_npair, nsweep
 integer :: max_niter = 1999
 integer :: ndiis = 7
 ! TODO: invoke DIIS after a specified cycle where the function value change is
 ! small.
 integer :: nfile = 0 ! nfile = 2*ndiis - 1
 integer, allocatable :: ijmap(:,:), eff_ijmap(:,:), rrmap(:,:,:)

 real(kind=8), parameter :: QPI = DATAN(1d0)     ! PI/4
 real(kind=8), parameter :: HPI = 2d0*DATAN(1d0) ! PI/2
 real(kind=8), parameter :: upd_thres = 1d-9
 ! determine whether to rotate (and update MOs, dipole integrals)
 ! If upd_thres is set to 1d-7, the total change of the target function might
 ! have non-negligible difference with the total change obtained at 1d-8.

 real(kind=8), parameter :: pop_thres = 0.7d0
 real(kind=8) :: dis_thres, conv_thres
 ! dis_thres: MO distances beyond the threshold will not undergo 2*2 rotations
 ! conv_thres: determine whether rotation/localization is converged
 ! Values of these two variables will be copied from loc()/pbc_loc() in
 !  $MOKIT_ROOT/mokit/lib/gaussian.py

 real(kind=8), allocatable :: u(:,:), cayley_k(:,:), k_old(:,:), k_diis(:,:)

 character(len=240), allocatable :: binfile(:) ! size nfile

 type :: rotation_index
  integer :: npair
  integer, allocatable :: pair_idx(:,:) ! size (2,npair)
 end type rotation_index
 type(rotation_index), allocatable :: rot_idx(:) ! size nsweep

 type :: mo_centers
  integer :: nc   ! the number of centers
  integer, allocatable :: idx(:)
 end type mo_centers
 type(mo_centers), allocatable :: mo_cen(:) ! size nmo

 logical :: diis = .true. ! .true./.false. for turning on/off DIIS

contains

! Find effective i-j map using the distance matrix of MOs, i.e. create/update
! the eff_ijmap array.
subroutine find_eff_ijmap(nmo, mo_dis)
 implicit none
 integer :: i, n
 integer, intent(in) :: nmo
 real(kind=8), intent(in) :: mo_dis(nmo*(nmo-1)/2)
 logical, allocatable :: mask(:)

 if(nmo < 1) then
  write(6,'(A)') 'ERROR in subroutine find_eff_ijmap: nmo<1.'
  write(6,'(A,I0)') 'nmo=', nmo
  stop
 end if

 allocate(mask(npair), source=(mo_dis<dis_thres))
 eff_npair = COUNT(mask)
 allocate(eff_ijmap(2,eff_npair))
 i = 0

 do n = 1, npair, 1
  if(mask(n)) then
   i = i + 1
   eff_ijmap(:,i) = ijmap(:,n)
   if(i == eff_npair) exit
  end if
 end do ! for n

 deallocate(mask)
! The following code also works, but no acceleration is observed.
! eff_ijmap = RESHAPE(PACK(ijmap,RESHAPE([mask,mask],[2,npair])), [2,eff_npair])
end subroutine find_eff_ijmap

! generate effective ordering within the orbital distance threshold
subroutine screen_jacobi_idx_by_dis(nmo, mo_dis)
 implicit none
 integer :: i, j, k, k1, k2, k3, m
 integer, intent(in) :: nmo
 real(kind=8), intent(in) :: mo_dis(nmo*(nmo-1)/2)

 if(.not. allocated(rot_idx)) allocate(rot_idx(nsweep))
 m = 2*nmo

!$omp parallel do schedule(dynamic) default(shared) private(i,j,k,k1,k2,k3)
 do i = 1, nsweep, 1
  k = 0
  do j = 1, np, 1
   k1 = MIN(rrmap(1,j,i), rrmap(2,j,i))
   k2 = MAX(rrmap(1,j,i), rrmap(2,j,i))
   if(k1 > 0) then
    k3 = (m-k1)*(k1-1)/2 + k2 - k1
    if(mo_dis(k3) < dis_thres) k = k + 1
   end if
  end do ! for j
  rot_idx(i)%npair = k
  if(allocated(rot_idx(i)%pair_idx)) deallocate(rot_idx(i)%pair_idx)
  if(k > 0) allocate(rot_idx(i)%pair_idx(2,k))
 end do ! for i
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(shared) private(i,j,k,k1,k2,k3)
 do i = 1, nsweep, 1
  k = 0
  do j = 1, np, 1
   k1 = MIN(rrmap(1,j,i), rrmap(2,j,i))
   k2 = MAX(rrmap(1,j,i), rrmap(2,j,i))
   if(k1 > 0) then
    k3 = (m-k1)*(k1-1)/2 + k2 - k1
    if(mo_dis(k3) < dis_thres) then
     k = k + 1
     rot_idx(i)%pair_idx(:,k) = [k1,k2]
    end if
   end if
  end do ! for j
 end do ! for i
!$omp end parallel do
end subroutine screen_jacobi_idx_by_dis

subroutine para22berry_kernel(nbf, nmo, coeff, mo_zdip, change)
 implicit none
 integer :: i, j, m, n, npair0
 integer, intent(in) :: nbf, nmo
 integer, allocatable :: idx(:,:)
 real(kind=8) :: rtmp, Aij, Bij, alpha, cc, ss, sc, cos_2a, sin_2a, sin_4a, &
  rchange
 real(kind=8), intent(inout) :: coeff(nbf,nmo)
 real(kind=8), intent(out) :: change
 real(kind=8), allocatable :: cos_a(:), sin_a(:)
 real(kind=8), external :: cmplx_v3_square, cmplx_v3_dot_real
 complex(kind=8) :: vdiff(3), vtmp(3,4)
 complex(kind=8), intent(inout) :: mo_zdip(3,nmo,nmo)
 logical, allocatable :: skip(:)

 change = 0d0

 do m = 1, nsweep, 1
  npair0 = rot_idx(m)%npair
  if(npair0 == 0) cycle
  allocate(idx(2,npair0), source=rot_idx(m)%pair_idx)
  allocate(cos_a(npair0), sin_a(npair0), skip(npair0))
  skip = .false.; rchange = 0d0

  !$omp parallel do schedule(dynamic) default(private) reduction(+:rchange) &
  !$omp shared(npair0,nbf,cos_a,sin_a,idx,mo_zdip,coeff,skip)
  do n = 1, npair0, 1
   i = idx(1,n); j = idx(2,n) ! i<j ensured when generating rot_idx
   vtmp(:,1) = mo_zdip(:,i,i)
   vtmp(:,2) = mo_zdip(:,i,j)
   vtmp(:,3) = mo_zdip(:,j,j)
   vdiff = vtmp(:,1) - vtmp(:,3)
   Aij = cmplx_v3_square(vtmp(:,2)) - 0.25d0*cmplx_v3_square(vdiff)
   Bij = cmplx_v3_dot_real(vdiff, vtmp(:,2))
   rtmp = HYPOT(Aij, Bij)
   sin_4a = Bij/rtmp
   rtmp = rtmp + Aij
   if(rtmp < upd_thres) then
    skip(n) = .true.
    cycle
   end if

   rchange = rchange + rtmp
   alpha = 0.25d0*DASIN(MAX(-1d0, MIN(sin_4a, 1d0)))
   if(Aij > 0d0) then
    alpha = QPI - alpha
   else if(Aij<0d0 .and. Bij<0d0) then
    alpha = HPI + alpha
   end if
   if(alpha > QPI) alpha = alpha - HPI
   cos_a(n) = DCOS(alpha); sin_a(n) = DSIN(alpha)

   ! update two orbitals
   call rotate_mo_ij(cos_a(n), sin_a(n), nbf, coeff(:,i), coeff(:,j))

   ! update corresponding dipole integrals
   cc = cos_a(n)*cos_a(n); ss = sin_a(n)*sin_a(n); sc = sin_a(n)*cos_a(n)
   cos_2a = 2d0*cc-1d0; sin_2a = 2d0*sc
   vtmp(:,4) = sin_2a*vtmp(:,2)
   mo_zdip(:,i,i) = cc*vtmp(:,1) + ss*vtmp(:,3) + vtmp(:,4)
   mo_zdip(:,j,j) = ss*vtmp(:,1) + cc*vtmp(:,3) - vtmp(:,4)
   mo_zdip(:,i,j) = cos_2a*vtmp(:,2) - sc*vdiff
  end do ! for n
  !$omp end parallel do

  ! update remaining off-diagonal elements of mo_zdip(:,x,y)
  do n = 1, npair0, 1
   if(skip(n)) cycle
   call update_up_tri_ij_zdip(nmo, idx(1,n), idx(2,n), cos_a(n), sin_a(n), &
                              mo_zdip)
  end do ! for n

  deallocate(idx, cos_a, sin_a, skip)
  change = change + rchange
 end do ! for m
end subroutine para22berry_kernel

subroutine para22boys_kernel(nbf, nmo, coeff, mo_dip, change)
 implicit none
 integer :: i, j, m, n, npair0
 integer, intent(in) :: nbf, nmo
 integer, allocatable :: idx(:,:)
 real(kind=8) :: ddot, rtmp, Aij, Bij, alpha, cc, ss, sc, cos_2a, sin_2a, &
  sin_4a, rchange, vdiff(3), vtmp(3,4)
 real(kind=8), intent(inout) :: coeff(nbf,nmo), mo_dip(3,nmo,nmo)
 real(kind=8), intent(out) :: change
 real(kind=8), allocatable :: cos_a(:), sin_a(:), tmp_mo(:)
 logical, allocatable :: skip(:)

 change = 0d0
 allocate(tmp_mo(nbf))

 do m = 1, nsweep, 1
  npair0 = rot_idx(m)%npair
  if(npair0 == 0) cycle
  allocate(idx(2,npair0), source=rot_idx(m)%pair_idx)
  allocate(cos_a(npair0), sin_a(npair0), skip(npair0))
  skip = .false.; rchange = 0d0

  !$omp parallel do schedule(dynamic) default(private) reduction(+:rchange) &
  !$omp shared(npair0,nbf,cos_a,sin_a,idx,mo_dip,coeff,skip)
  do n = 1, npair0, 1
   i = idx(1,n); j = idx(2,n) ! i<j ensured when generating rot_idx
   vtmp(:,1) = mo_dip(:,i,i)
   vtmp(:,2) = mo_dip(:,i,j)
   vtmp(:,3) = mo_dip(:,j,j)
   vdiff = vtmp(:,1) - vtmp(:,3)
   Aij = ddot(3,vtmp(:,2),1,vtmp(:,2),1) - 0.25d0*ddot(3,vdiff,1,vdiff,1)
   Bij = ddot(3, vdiff, 1, vtmp(:,2), 1)
   rtmp = HYPOT(Aij, Bij)
   sin_4a = Bij/rtmp
   rtmp = rtmp + Aij
   if(rtmp < upd_thres) then
    skip(n) = .true.
    cycle
   end if

   rchange = rchange + rtmp
   alpha = 0.25d0*DASIN(MAX(-1d0, MIN(sin_4a, 1d0)))
   if(Aij > 0d0) then
    alpha = QPI - alpha
   else if(Aij<0d0 .and. Bij<0d0) then
    alpha = HPI + alpha
   end if
   if(alpha > QPI) alpha = alpha - HPI
   cos_a(n) = DCOS(alpha); sin_a(n) = DSIN(alpha)

   ! update two orbitals
   call rotate_mo_ij(cos_a(n), sin_a(n), nbf, coeff(:,i), coeff(:,j))

   ! update corresponding dipole integrals
   cc = cos_a(n)*cos_a(n); ss = sin_a(n)*sin_a(n); sc = sin_a(n)*cos_a(n)
   cos_2a = 2d0*cc-1d0; sin_2a = 2d0*sc
   vtmp(:,4) = sin_2a*vtmp(:,2)
   mo_dip(:,i,i) = cc*vtmp(:,1) + ss*vtmp(:,3) + vtmp(:,4)
   mo_dip(:,j,j) = ss*vtmp(:,1) + cc*vtmp(:,3) - vtmp(:,4)
   mo_dip(:,i,j) = cos_2a*vtmp(:,2) - sc*vdiff
  end do ! for n
  !$omp end parallel do

  ! update remaining off-diagonal elements of mo_dip(:,x,y)
  do n = 1, npair0, 1
   if(skip(n)) cycle
   call update_up_tri_ij_dip(nmo, idx(1,n), idx(2,n), cos_a(n), sin_a(n), &
                             mo_dip)
  end do ! for n

  deallocate(idx, cos_a, sin_a, skip)
  change = change + rchange
 end do ! for m

 deallocate(tmp_mo)
end subroutine para22boys_kernel

subroutine para22pm_kernel(nbf, nmo, natom, coeff, gross, change)
 implicit none
 integer :: i, j, m, n, npair0
 integer, intent(in) :: nbf, nmo, natom
 integer, allocatable :: idx(:,:)
 real(kind=8) :: ddot, rtmp, Aij, Bij, alpha, sin_4a, rchange
 real(kind=8), intent(inout) :: coeff(nbf,nmo), gross(natom,nmo,nmo)
 real(kind=8), intent(out) :: change
 real(kind=8), allocatable :: vdiff(:), vtmp(:,:), cos_a(:), sin_a(:)
 logical, allocatable :: skip(:)

 change = 0d0
 allocate(vdiff(natom), vtmp(natom,3))

 do m = 1, nsweep, 1
  npair0 = rot_idx(m)%npair
  if(npair0 == 0) cycle
  allocate(idx(2,npair0), source=rot_idx(m)%pair_idx)
  allocate(cos_a(npair0), sin_a(npair0), skip(npair0))
  skip = .false.; rchange = 0d0

  !$omp parallel do schedule(dynamic) default(private) reduction(+:rchange) &
  !$omp shared(natom,npair0,nbf,cos_a,sin_a,idx,gross,coeff,skip)
  do n = 1, npair0, 1
   i = idx(1,n); j = idx(2,n) ! i<j ensured when generating rot_idx
   call dcopy(natom, gross(:,i,i), 1, vtmp(:,1), 1)
   call dcopy(natom, gross(:,i,j), 1, vtmp(:,2), 1)
   call dcopy(natom, gross(:,j,j), 1, vtmp(:,3), 1)
   call dcopy(natom, vtmp(:,1)-vtmp(:,3), 1, vdiff, 1)
   Aij = ddot(natom,vtmp(:,2),1,vtmp(:,2),1) - 0.25d0*ddot(natom,vdiff,1,vdiff,1)
   Bij = ddot(natom, vdiff, 1, vtmp(:,2), 1)
   rtmp = HYPOT(Aij, Bij)
   sin_4a = Bij/rtmp
   rtmp = rtmp + Aij
   if(rtmp < upd_thres) then
    skip(n) = .true.
    cycle
   end if

   rchange = rchange + rtmp
   alpha = 0.25d0*DASIN(MAX(-1d0, MIN(sin_4a, 1d0)))
   if(Aij > 0d0) then
    alpha = QPI - alpha
   else if(Aij<0d0 .and. Bij<0d0) then
    alpha = HPI + alpha
   end if
   if(alpha > QPI) alpha = alpha - HPI
   cos_a(n) = DCOS(alpha); sin_a(n) = DSIN(alpha)

   call rotate_mo_ij(cos_a(n), sin_a(n), nbf, coeff(:,i), coeff(:,j))
   call update_gross_ii_jj_ji(cos_a(n), sin_a(n), natom, vtmp, vdiff, &
                              gross(:,i,i), gross(:,j,j), gross(:,i,j))
  end do ! for n
  !$omp end parallel do

  ! update remaining off-diagonal elements of gross(:,x,y)
  do n = 1, npair0, 1
   if(skip(n)) cycle
   call update_up_tri_ij_gross(natom, nmo, idx(1,n), idx(2,n), cos_a(n), &
                               sin_a(n), gross)
  end do ! for n

  deallocate(idx, cos_a, sin_a, skip)
  change = change + rchange
 end do ! for m

 deallocate(vdiff, vtmp)
end subroutine para22pm_kernel

end module lo_info

subroutine init_lo_diis(proname, nmo)
 use lo_info, only: ndiis, nfile, binfile, u, cayley_k, k_old, k_diis
 implicit none
 integer :: i
 integer, intent(in) :: nmo
 character(len=240), intent(in) :: proname

 nfile = 2*ndiis - 1
 allocate(binfile(nfile))

 write(binfile(1),'(A)') TRIM(proname)//'.N'

 do i = 2, ndiis
  write(binfile(i),'(A,I0)') TRIM(proname)//'.N-',i-1
 end do ! for i

 do i = ndiis+1, nfile, 1
  write(binfile(i),'(A,I0)') TRIM(proname)//'.aN-',i-ndiis+1
 end do ! for i

 call delete_files(nfile, binfile)
 allocate(u(nmo,nmo), cayley_k(nmo,nmo), k_old(nmo,nmo))
 call init_identity_mat(nmo, u)
 cayley_k = 0d0

 allocate(k_diis(ndiis+1,ndiis+1), source=0d0)
 k_diis(1:ndiis,ndiis+1) = 1d0
 k_diis(ndiis+1,1:ndiis) = 1d0
end subroutine init_lo_diis

! Get/Find the center of each MO using SCPA(C-squared Population Analysis)
subroutine get_mo_center_by_scpa(natom, nbf, nmo, bfirst, mo, mo_cen)
 use lo_info, only: pop_thres, mo_centers
 implicit none
 integer :: i, j, k, m, i1, i2, ak(1)
 integer, intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
 real(kind=8) :: r, min_v, sum_r
 real(kind=8), intent(in) :: mo(nbf,nmo)
 real(kind=8), allocatable :: tmp_mo(:), pop(:), gross(:,:)
 type(mo_centers), intent(out) :: mo_cen(nmo)

 allocate(tmp_mo(nbf), pop(natom), gross(natom,nmo))

!$omp parallel do schedule(static) default(private) &
!$omp shared(natom,nmo,mo,bfirst,gross)
 do i = 1, nmo, 1
  tmp_mo = mo(:,i)**2
  tmp_mo = tmp_mo/SUM(tmp_mo)
  do j = 1, natom, 1
   i1 = bfirst(j); i2 = bfirst(j+1) - 1
   pop(j) = SUM(tmp_mo(i1:i2))
  end do ! for j
  gross(:,i) = pop
 end do ! for i
!$omp end parallel do
 deallocate(tmp_mo)

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(natom,nmo,bfirst,gross,mo_cen)
 do i = 1, nmo, 1
  ! find the largest component on an atom of an MO
  pop = gross(:,i)
  mo_cen(i)%nc = 1; m = 1
  ak = MAXLOC(pop); k = ak(1); r = pop(k)

  ! check the 2nd largest component
  if(r < pop_thres) then
   ak = MINLOC(pop); m = ak(1); min_v = pop(m)
   pop(k) = min_v-0.5d0; m = 1; sum_r = r
   do j = 1, natom-1, 1
    m = m + 1; ak = MAXLOC(pop); k = ak(1)
    r = pop(k); sum_r = sum_r + r
    if(sum_r > pop_thres) exit
    pop(k) = min_v - 0.5d0
   end do ! for j
  end if

  mo_cen(i)%nc = m
  if(allocated(mo_cen(i)%idx)) deallocate(mo_cen(i)%idx)
  allocate(mo_cen(i)%idx(m))
 end do ! for i
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(natom,nmo,bfirst,gross,mo_cen)
 do i = 1, nmo, 1
  ! find the largest component on an atom of an MO
  pop = gross(:,i)
  ak = MAXLOC(pop); k = ak(1)
  r = pop(k); mo_cen(i)%idx(1) = k

  ! check the 2nd largest component
  if(r < pop_thres) then
   ak = MINLOC(pop); m = ak(1); min_v = pop(m)
   pop(k) = min_v-0.5d0; m = 1; sum_r = r
   do j = 1, natom-1, 1
    ak = MAXLOC(pop); k = ak(1); r = pop(k)
    m = m + 1; mo_cen(i)%idx(m) = k
    sum_r = sum_r + r
    if(sum_r > pop_thres) exit
    pop(k) = min_v - 0.5d0
   end do ! for j
  end if
 end do ! for i
!$omp end parallel do

 deallocate(pop, gross)
end subroutine get_mo_center_by_scpa

! This subroutine is the same as subroutine get_mo_center_from_diag_gross above,
!  except that here gross is a 3D array, and only gross(:,i,i) will be used.
subroutine get_mo_center_from_gross(natom, nmo, gross, mo_cen)
 use lo_info, only: pop_thres, mo_centers
 implicit none
 integer :: i, j, k, m, ak(1)
 integer, intent(in) :: natom, nmo
 real(kind=8) :: r, min_v, sum_r
 real(kind=8), intent(in) :: gross(natom,nmo,nmo)
 real(kind=8), allocatable :: pop(:)
 type(mo_centers), intent(out) :: mo_cen(nmo)

 allocate(pop(natom))

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(natom,nmo,gross,mo_cen)
 do i = 1, nmo, 1
  ! find the largest component on an atom of an MO
  pop = gross(:,i,i)
  mo_cen(i)%nc = 1; m = 1
  ak = MAXLOC(pop); k = ak(1); r = pop(k)

  ! check the 2nd largest component
  if(r < pop_thres) then
   ak = MINLOC(pop); m = ak(1); min_v = pop(m)
   pop(k) = min_v-0.5d0; m = 1; sum_r = r
   do j = 1, natom-1, 1
    m = m + 1; ak = MAXLOC(pop); k = ak(1)
    r = pop(k); sum_r = sum_r + r
    if(sum_r > pop_thres) exit
    pop(k) = min_v - 0.5d0
   end do ! for j
  end if

  mo_cen(i)%nc = m
  if(allocated(mo_cen(i)%idx)) deallocate(mo_cen(i)%idx)
  allocate(mo_cen(i)%idx(m))
 end do ! for i
!$omp end parallel do

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(natom,nmo,gross,mo_cen)
 do i = 1, nmo, 1
  ! find the largest component on an atom of an MO
  pop = gross(:,i,i)
  ak = MAXLOC(pop); k = ak(1)
  r = pop(k); mo_cen(i)%idx(1) = k

  ! check the 2nd largest component
  if(r < pop_thres) then
   ak = MINLOC(pop); m = ak(1); min_v = pop(m)
   pop(k) = min_v-0.5d0; m = 1; sum_r = r
   do j = 1, natom-1, 1
    ak = MAXLOC(pop); k = ak(1); r = pop(k)
    m = m + 1; mo_cen(i)%idx(m) = k
    sum_r = sum_r + r
    if(sum_r > pop_thres) exit
    pop(k) = min_v - 0.5d0
   end do ! for j
  end if
 end do ! for i
!$omp end parallel do

 deallocate(pop)
end subroutine get_mo_center_from_gross

! find/get the distance matrix of MOs from the distance matrix of atoms
subroutine atm_dis2mo_dis(natom, nmo, npair, dis, mo_cen, mo_dis)
 use lo_info, only: ijmap, mo_centers
 implicit none
 integer :: i, j, k, m, n, k1, nc1, nc2
 integer, intent(in) :: natom, nmo, npair
 real(kind=8) :: r, min_dis
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(out) :: mo_dis(npair)
 type(mo_centers), intent(in) :: mo_cen(nmo)

 if(ANY(mo_cen(:)%nc == 0)) then
  write(6,'(/,A)') 'ERROR in subroutine atm_dis2mo_dis: some MO does not have i&
                   &ts center.'
  write(6,'(A,2I7)') 'natom, nmo=', natom, nmo
  write(6,'(A)') 'mo_cen(:)%nc='
  write(6,'(12I7)') mo_cen(:)%nc
  stop
 end if

!$omp parallel do schedule(dynamic) default(private) &
!$omp shared(npair,ijmap,mo_cen,dis,mo_dis)
 do n = 1, npair, 1
  i = ijmap(1,n); j = ijmap(2,n)
  min_dis = dis(mo_cen(i)%idx(1), mo_cen(j)%idx(1))
  nc1 = mo_cen(i)%nc; nc2 = mo_cen(j)%nc

  do k = 1, nc1, 1
   k1 = mo_cen(i)%idx(k)
   do m = 1, nc2, 1
    r = dis(mo_cen(j)%idx(m), k1)
    if(r < min_dis) min_dis = r
   end do ! for m
  end do ! for k

  mo_dis(n) = min_dis
 end do ! for n
!$omp end parallel do
end subroutine atm_dis2mo_dis

! open a GAMESS .inp/.dat file and jump to the $DATA section
subroutine goto_data_section_in_gms_inp(inpname, fid)
 implicit none
 integer :: i
 integer, intent(out) :: fid
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical :: alive

 fid = 0
 inquire(file=TRIM(inpname),opened=alive)

 if(alive) then
  inquire(file=TRIM(inpname),number=fid)
  rewind(fid)
 else
  open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:6))
   if(buf(3:6) == 'DATA') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine goto_data_section_in_gms_inp: wrong for&
                   &mat in file '//TRIM(inpname)
  close(fid)
  stop
 end if
 ! Note: this file is still opened.
end subroutine goto_data_section_in_gms_inp

! copy GVB CI coefficients (or called pair coefficients) from a .dat file into
! another one 
subroutine copy_and_add_pair_coeff(addH_dat, datname, nopen)
 implicit none
 integer :: i, j, npair, fid1, fid2, fid3, RENAME
 integer, intent(in) :: nopen
!f2py intent(in) :: nopen
 character(len=240) :: buf, new_dat
 character(len=240), intent(in) :: addH_dat, datname
!f2py intent(in) :: addH_dat, datname

 i = INDEX(addH_dat, '.dat', back=.true.)
 new_dat = addH_dat(1:i-1)//'.t'

 call goto_data_section_in_gms_inp(addH_dat, fid1)

 open(newunit=fid2,file=TRIM(new_dat),status='replace')
 write(fid2,'(A)') ' $DATA'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:5) == '$VEC') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 open(newunit=fid3,file=TRIM(datname),status='old',position='rewind')
 do while(.true.)
  read(fid3,'(A)') buf
  if(buf(2:5) == '$SCF') exit
 end do ! for while

 i = INDEX(buf, 'CICOEF')
 if(i > 0) then
  write(fid2,'(A)') TRIM(buf)
  npair = 1
 else
  npair = 0
 end if

 do while(.true.)
  read(fid3,'(A)') buf
  i = INDEX(buf, 'CICOEF')
  if(i > 0) then
   j = INDEX(buf, '$END')
   if(j > 0) then
    buf(j:j+3) = '    '
    write(fid2,'(A)') TRIM(buf)
    exit
   else
    write(fid2,'(A)') TRIM(buf)
   end if
   npair = npair + 1
  else
   exit
  end if
 end do ! for while

 close(fid3)
 do i = 1, nopen, 1
  j = 2*(npair+i) - 1
  write(fid2,'(3X,A,I3,A)') 'CICOEF(',j,')= 0.7071067811865476,-0.7071067811865476'
 end do ! for i
 write(fid2,'(A,/,A)') ' $END', ' $VEC'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(new_dat), TRIM(addH_dat))
end subroutine copy_and_add_pair_coeff

! check whether UHF in a given GAMESS .inp file
subroutine check_uhf_in_gms_inp(inpname, uhf)
 implicit none
 integer :: fid
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical, intent(out) :: uhf

 uhf = .false.
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)') buf
  call upper(buf)
  if(INDEX(buf,'SCFTYP=UHF') > 0) then
   uhf = .true.
   exit
  end if
  if(INDEX(buf,'$END') > 0) exit
 end do ! for while

 close(fid)
end subroutine check_uhf_in_gms_inp

! read charge, spin multiplicity, uhf(bool) bohrs(bool) from a given GAMESS
! .inp/.dat file
subroutine read_charge_mult_isph_from_gms_inp(inpname, charge, mult, isph, uhf,&
                                              ghf, ecp)
 implicit none
 integer :: i, fid
 integer, intent(out) :: charge, mult, isph
 character(len=5) :: str5
 character(len=240) :: buf
 character(len=480) :: buf1
 character(len=240), intent(in) :: inpname
 logical, intent(out) :: uhf, ghf, ecp

 charge = 0; mult = 1; isph = -1
 uhf = .false.; ghf = .false.; ecp = .false.

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:8))
   if(buf(3:8) == 'CONTRL') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_charge_mult_isph_from_gms_inp: no &
                   &$CONTRL section'
  write(6,'(A)') 'found in file '//TRIM(inpname)
  close(fid)
  stop
 end if

 i = LEN_TRIM(buf)
 if(buf(i-3:i) == '$END') then
  buf1 = TRIM(buf)
 else
  read(fid,'(A)') buf1
  buf1 = TRIM(buf)//' '//TRIM(buf1)
 end if
 call upper(buf1)
 if(INDEX(buf1,'UHF') > 0) uhf = .true.
 if(INDEX(buf1,'GHF') > 0) ghf = .true.

 i = INDEX(buf1, 'ICHARG=')
 if(i > 0) read(buf1(i+7:),*) charge

 i = INDEX(buf1, 'MULT=')
 if(i > 0) read(buf1(i+5:),*) mult

 i = INDEX(buf1, 'ISPHER=')
 if(i > 0) read(buf1(i+7:),*) isph

 do while(.true.)
  read(fid,'(A)',iostat=i) str5
  if(i /= 0) exit

  if(str5(2:2) == '$') then
   call upper(str5(3:5))
   if(str5(3:5) == 'ECP') then
    ecp = .true.
    exit
   end if
  end if
 end do ! for while

 close(fid)
end subroutine read_charge_mult_isph_from_gms_inp

! read elements, nuclear charges and Cartesian coordinates from a GAMESS .inp file
subroutine read_elem_nuc_coor_from_gms_inp(inpname, natom, elem, nuc, coor, ghost)
 use phys_cons, only: Bohr_const
 implicit none
 integer :: i, k, fid, nline
 integer, intent(in) :: natom
 integer, intent(out) :: nuc(natom)
 real(kind=8), allocatable :: nuc1(:)
 real(kind=8), intent(out) :: coor(3,natom)
 character(len=1) :: str
 character(len=2), intent(out) :: elem(natom)
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical :: bohrs
 logical, intent(out) :: ghost(natom) ! size natom

 allocate(nuc1(natom), source=0d0)
 nuc = 0; coor = 0d0; elem = ' '; ghost = .false.

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 ! find in the first 6 lines whether the coordinates are in Angstrom or Bohr
 bohrs = .false.
 do i = 1, 6
  read(fid,'(A)') buf
  call upper(buf)
  if(index(buf,'UNITS=BOHR') /= 0) then
   bohrs = .true.
   exit
  end if
  if(index(buf,'$END') /= 0) exit
 end do ! for i
 ! Angstrom/Bohr determined

 call goto_data_section_in_gms_inp(inpname, fid)
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do k = 1, natom, 1
  read(fid,*) elem(k), nuc1(k), coor(1:3,k)
  if(nuc1(k) < 0d0) then ! ghost atom (0 charge, has basis function)
   nuc1(k) = -nuc1(k)
   ghost(k) = .true.
  end if

  do while(.true.)
   read(fid,'(A)') buf
   read(buf,*,iostat=i) str, nline
   if(i /= 0) exit

   do i = 1, nline, 1
    read(fid,'(A)') buf
   end do ! for do
  end do ! for while

 end do ! for k

 close(fid)

 call standardize_elem(natom, elem)
 forall(i = 1:natom) nuc(i) = NINT(nuc1(i))
 deallocate(nuc1)
 if(bohrs) coor = coor*Bohr_const
end subroutine read_elem_nuc_coor_from_gms_inp

! read nbf and nif from a given GAMESS .inp/.dat file
subroutine read_nbf_and_nif_from_gms_inp(inpname, nbf, nif)
 implicit none
 integer:: i, j, k, fid
 integer, intent(out) :: nbf, nif
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 nbf = 0; nif = 0
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf)
  j = INDEX(buf,'$GUESS'); k = INDEX(buf,'NORB=')
  if(j/=0 .and. k/=0) then
   read(buf(k+5:),*) nif
   exit
  end if
 end do ! for while

 if(i /= 0) then
  close(fid)
  return
 end if

 call goto_data_section_in_gms_inp(inpname, fid)
 ! read the Title Card line, find nbf in it
 read(fid,'(A)') buf
 i = INDEX(buf, 'nbf=')
 if(i == 0) then ! if not found, set to nif
  nbf = nif
 else
  read(buf(i+4:),*) nbf
 end if

 close(fid)
end subroutine read_nbf_and_nif_from_gms_inp

! read na, nb, nif and nbf from a given GAMESS .inp file
! Note: when spherical harmonic functions are used, the nbf here will <=
!  the number of basis functions in $VEC (where MOs are always expanded
!  on Cartesian functions)
subroutine read_na_nb_nif_nbf_from_gms_inp(inpname, na, nb, nif, nbf)
 implicit none
 integer :: i, fid
 integer, intent(out) :: na, nb, nif, nbf
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 na = 0; nb = 0; nif = 0; nbf = 0

 call goto_data_section_in_gms_inp(inpname, fid)
 read(fid,'(A)') buf ! read the Title Card line
 close(fid)

 i = INDEX(buf,'nbf=')
 read(buf(i+4:),*) nbf
 buf(i-1:) = ' '

 i = INDEX(buf,'nif=')
 read(buf(i+4:),*) nif
 buf(i-1:) = ' '

 i = INDEX(buf,'nb=')
 read(buf(i+3:),*) nb
 buf(i-1:) = ' '

 i = INDEX(buf,'na=')
 read(buf(i+3:),*) na
end subroutine read_na_nb_nif_nbf_from_gms_inp

! read type all_ecp from a given GAMESS .inp/.dat file
subroutine read_all_ecp_from_gms_inp(inpname)
 use pg, only: natom, all_ecp, ecp_exist
 implicit none
 integer :: i, j, k, m, n, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 if(natom == 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_all_ecp_from_gms_inp: natom = 0.'
  write(6,'(A)') 'The variable natom should be initialized before calling this &
                 &subroutine.'
  stop
 end if

 allocate(all_ecp(natom))
 all_ecp(:)%ecp = .false.
 if(.not. ecp_exist) return

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:2) == '$') then
   call upper(buf(3:5))
   if(buf(3:5) == 'ECP') exit
  end if
 end do ! for while

 do i = 1, natom, 1
  read(fid,'(A)') buf
  call upper(buf)
  if(index(buf,'NONE') /= 0) cycle

  all_ecp(i)%ecp = .true.
  k = INDEX(buf,'GEN')
  read(buf(k+3:),*) all_ecp(i)%core_e, m
  all_ecp(i)%highest = m
  allocate(all_ecp(i)%potential(m+1))

  do j = 1, m+1, 1
   read(fid,'(A)') buf
   if(j == 1) then
    k = INDEX(buf,'ul')
    if(k == 0) then
     write(6,'(/,A)') 'ERROR in subroutine read_all_ecp_from_gms_inp: ECP/PP do&
                      &s not start'
     write(6,'(A)') "with '-ul potential'. You should check the format of ECP/P&
                    &P data in file "//TRIM(inpname)
     close(fid)
     stop
    end if
   end if

   read(buf,*) n
   all_ecp(i)%potential(j)%n = n
   allocate(all_ecp(i)%potential(j)%col1(n), source=0d0)
   allocate(all_ecp(i)%potential(j)%col2(n), source=0)
   allocate(all_ecp(i)%potential(j)%col3(n), source=0d0)
   do k = 1, n, 1
    read(fid,*) all_ecp(i)%potential(j)%col1(k),all_ecp(i)%potential(j)%col2(k),&
    all_ecp(i)%potential(j)%col3(k)
   end do ! for k
  end do ! for j
 end do ! for i

 close(fid)
end subroutine read_all_ecp_from_gms_inp

! deallocate the allocatable arrays in array prim_gau
subroutine clear_prim_gau()
 use pg, only: prim_gau
 implicit none
 integer :: i

 do i = 1, 7, 1
  prim_gau(i)%stype = ' '
  prim_gau(i)%nline = 0
  prim_gau(i)%ncol = 0
  if(allocated(prim_gau(i)%coeff)) deallocate(prim_gau(i)%coeff)
 end do ! for i
end subroutine clear_prim_gau

! read this type of primitive gaussians, i.e., 'S', 'L', etc.
subroutine read_prim_gau(stype, nline, fid)
 use pg, only: prim_gau
 implicit none
 integer :: i, j, k, itmp, ncol, ncol1, nline0, nline1
 integer, intent(in) :: nline, fid
 real(kind=8), parameter :: zero = 1d-5
 real(kind=8) :: exp0, rtmp, rtmp1
 real(kind=8), allocatable :: coeff(:,:)
 character(len=1), intent(in) :: stype
 logical :: share

 call stype2itype(stype, k)
 ! two cases: 'L', or not 'L'

 if(k /= 0) then ! 'S','P','D','F','G','H','I', not 'L'

  ! two subcases: whether or not this angular momentum has occurred
  if(allocated(prim_gau(k)%coeff)) then
   read(fid,*) itmp, exp0, rtmp
   BACKSPACE(fid)

   nline0 = prim_gau(k)%nline
   ncol = prim_gau(k)%ncol
   allocate(coeff(nline0,ncol), source=prim_gau(k)%coeff)
   deallocate(prim_gau(k)%coeff)

   ! further two sub-subcases: whether sharing the same exponents
   share = .false. ! initialization
   if(DABS(coeff(1,1)-exp0) < zero) then
    nline1 = MAX(nline0, nline)
    share = .true.
   end if
   if((.not.share) .and. nline==1 .and. DABS(coeff(nline0,1)-exp0)<zero) then
    ! this rare case occurs in nobasistransform cc-pVDZ for Zn
    nline1 = nline0
    share = .true.
   end if
   if(.not. share) nline1 = nline0 + nline

   allocate(prim_gau(k)%coeff(nline1,ncol+1),source=0d0)
   prim_gau(k)%coeff(1:nline0,1:ncol) = coeff
   deallocate(coeff)
   ncol = ncol + 1
   prim_gau(k)%ncol = ncol
   prim_gau(k)%nline = nline1

   do i = 1, nline, 1
    read(fid,*) itmp, exp0, rtmp
    if(nline == 1) rtmp = DSIGN(1d0, rtmp)
    if(share) then
     if(i > nline0) then
      prim_gau(k)%coeff(i,1) = exp0
      prim_gau(k)%coeff(i,ncol) = rtmp
     else ! i <= nline0
      do j = 1, nline0, 1
       if(DABS(prim_gau(k)%coeff(j,1)-exp0) < zero) then
        prim_gau(k)%coeff(j,ncol) = rtmp
        exit
       end if
      end do ! for j
      if(j == nline0+1) then
       write(6,'(/,A)') 'ERROR in subroutine read_prim_gau: unexpected basis data.'
       write(6,'(A)') 'Did you forget to add `int=nobasistransform` in .gjf file?'
       write(6,'(3(A,I0))') 'stype='//stype//', i=', i, ', nline=', nline, &
                            ', nline0=',nline0
       write(6,'(2(A,E16.8))') 'exp0=', exp0, ', rtmp=', rtmp
       stop
      end if
     end if
    else ! not share
     prim_gau(k)%coeff(nline0+i,1) = exp0
     prim_gau(k)%coeff(nline0+i,ncol) = rtmp
    end if
   end do ! for i

  else ! never occurs before, read first time
   prim_gau(k)%stype = stype
   prim_gau(k)%nline = nline
   prim_gau(k)%ncol = 2
   allocate(prim_gau(k)%coeff(nline,2), source=0d0)
   do i = 1, nline, 1
    read(fid,*) itmp, prim_gau(k)%coeff(i,1), rtmp
    if(nline == 1) rtmp = DSIGN(1d0, rtmp)
    prim_gau(k)%coeff(i,2) = rtmp
   end do ! for i
  end if

 else ! k == 0, this is 'L' for Pople-type basis sets
  ! sharing exponents for 'L' is unlikely in Pople basis sets, don't consider

  if(allocated(prim_gau(1)%coeff)) then ! 'S'
   nline0 = prim_gau(1)%nline
   ncol = prim_gau(1)%ncol
   prim_gau(1)%nline = nline0 + nline
   prim_gau(1)%ncol = ncol + 1
   allocate(coeff(nline0,ncol), source=prim_gau(1)%coeff)
   deallocate(prim_gau(1)%coeff)
   allocate(prim_gau(1)%coeff(nline0+nline,ncol+1),source=0d0)
   prim_gau(1)%coeff(1:nline0,1:ncol) = coeff
   deallocate(coeff)
   ncol = ncol + 1
  else
   prim_gau(1)%stype = 'S'
   prim_gau(1)%nline = nline
   prim_gau(1)%ncol = 2
   allocate(prim_gau(1)%coeff(nline,2),source=0d0)
   nline0 = 0; ncol = 2
  end if

  if(allocated(prim_gau(2)%coeff)) then ! 'P'
   nline1 = prim_gau(2)%nline
   ncol1 = prim_gau(2)%ncol
   prim_gau(2)%nline = nline1 + nline
   prim_gau(2)%ncol = ncol1 + 1
   allocate(coeff(nline1,ncol1), source=prim_gau(2)%coeff)
   deallocate(prim_gau(2)%coeff)
   allocate(prim_gau(2)%coeff(nline1+nline,ncol1+1),source=0d0)
   prim_gau(2)%coeff(1:nline1,1:ncol1) = coeff
   deallocate(coeff)
   ncol1 = ncol1 + 1
  else
   prim_gau(2)%stype = 'P'
   prim_gau(2)%nline = nline
   prim_gau(2)%ncol = 2
   allocate(prim_gau(2)%coeff(nline,2),source=0d0)
   nline1 = 0; ncol1 = 2
  end if

  do i = 1, nline, 1
   read(fid,*) itmp, exp0, rtmp, rtmp1
   if(nline == 1) then
    rtmp = DSIGN(1d0, rtmp)
    rtmp1 = DSIGN(1d0, rtmp1)
   end if
   prim_gau(1)%coeff(nline0+i,1) = exp0
   prim_gau(1)%coeff(nline0+i,ncol) = rtmp
   prim_gau(2)%coeff(nline1+i,1) = exp0
   prim_gau(2)%coeff(nline1+i,ncol1) = rtmp1
  end do ! for i
 end if
end subroutine read_prim_gau

! determine the highest angular momentum quantum number
subroutine get_highest_am()
 use pg, only: highest, prim_gau
 implicit none

 if(allocated(prim_gau(7)%coeff)) then
  highest = 6
 else if(allocated(prim_gau(6)%coeff)) then
  highest = 5
 else if(allocated(prim_gau(5)%coeff)) then
  highest = 4
 else if(allocated(prim_gau(4)%coeff)) then
  highest = 3
 else if(allocated(prim_gau(3)%coeff)) then
  highest = 2
 else if(allocated(prim_gau(2)%coeff)) then
  highest = 1
 else
  highest = 0
 end if
end subroutine get_highest_am

! print primitive gaussians into an OpenMolcas .input file
subroutine prt_prim_gau2molcas_inp(iatom, fid)
 use pg, only: prim_gau, nuc, highest, all_ecp, elem, ecp_exist, ghost
 implicit none
 integer :: i, j, k, m, n, nline, ncol
 integer, intent(in) :: iatom, fid
 integer, allocatable :: list(:)
 character(len=1), parameter :: am(0:6) = ['S','P','D','F','G','H','I']

 call get_highest_am()

 if(ghost(iatom)) then
  write(fid,'(5X,A2,3X,I1)') '0.', highest
 else
  if(ecp_exist) then
   write(fid,'(5X,I0,A1,3X,I1)') nuc(iatom)-all_ecp(iatom)%core_e,'.',highest
  else
   write(fid,'(5X,I0,A1,3X,I1)') nuc(iatom), '.', highest
  end if
 end if

 do i = 1, 7, 1
  if(.not. allocated(prim_gau(i)%coeff)) cycle
  write(fid,'(A)') '* '//prim_gau(i)%stype//'-type functions'
  nline = prim_gau(i)%nline
  ncol = prim_gau(i)%ncol
  write(fid,'(2(1X,I4))') nline, ncol-1
  do j = 1, nline, 1
   write(fid,'(3X,ES16.9)') prim_gau(i)%coeff(j,1)
  end do ! for j
  do j = 1, nline, 1
   write(fid,'(10(ES16.9,2X))') (prim_gau(i)%coeff(j,k), k=2,ncol)
  end do ! for j
 end do ! for i

 if(.not. ecp_exist) return

 if(all_ecp(iatom)%ecp) then
  m = all_ecp(iatom)%highest
  write(fid,'(A2,1X,A2,1X,I0,1X,I0)') 'PP', elem(iatom), all_ecp(iatom)%core_e, m
  allocate(list(m+1))
  list(1) = m
  forall(i=2:m+1) list(i) = i-2

  do i = 1, m+1, 1
   n = all_ecp(iatom)%potential(i)%n
   if(i == 1) then
    write(fid,'(I0,A)') n,';!'//am(list(1))//' POTENTIAL'
   else ! i > 1
    write(fid,'(I0,A)') n,';!'//am(list(i))//'-'//am(list(1))//' POTENTIAL'
   end if

   do j = 1, n, 1
    write(fid,'(I0,2(1X,F16.8))') all_ecp(iatom)%potential(i)%col2(j),&
     all_ecp(iatom)%potential(i)%col3(j), all_ecp(iatom)%potential(i)%col1(j)
   end do ! for j
  end do ! for i
  write(fid,'(A1,/,A,/,A)') '*', 'Spectral', 'End of Spectral'
 end if
end subroutine prt_prim_gau2molcas_inp

! update the number of times each atom occurred
subroutine update_ntimes(iatom)
 use pg, only: elem, ntimes
 implicit none
 integer :: i
 integer, intent(in) :: iatom
 character(len=2) :: ctmp

 ctmp = elem(iatom)

 do i = iatom-1, 1, -1
  if(ctmp == elem(i)) then
   ntimes(iatom) = ntimes(i) + 1
   exit
  end if
 end do ! for i
end subroutine update_ntimes

! update the number of times each atom occurred
subroutine calc_ntimes(natom, elem, ntimes)
 implicit none
 integer :: i, j
 integer, intent(in) :: natom
 integer, intent(out) :: ntimes(natom)
 character(len=2) :: tmp
 character(len=2), intent(in) :: elem(natom)

 ntimes = 1

 do i = 2, natom, 1
  tmp = elem(i)

  do j = i-1, 1, -1
   if(tmp == elem(j)) then
    ntimes(i) = ntimes(j) + 1
    exit
   end if
  end do ! for j

 end do ! for i
end subroutine calc_ntimes

! generate contracted string, e.g. 5s3p1d -> 3s2p1d
subroutine gen_contracted_string(nline, ncol, str1, str2)
 implicit none
 integer :: i
 integer, intent(in) :: nline(7), ncol(7)
 character(len=1), parameter :: am(7) = ['s','p','d','f','g','h','i']
 character(len=3) :: str
 character(len=21), intent(out) :: str1, str2
 !10s10p10d10f10g10h10i

 str1 = ' '; str2 = ' '

 do i = 1, 7, 1
  if(nline(i) > 0) then
   write(str,'(I0,A1)') nline(i), am(i)
   str1 = TRIM(str1)//TRIM(str)
   write(str,'(I0,A1)') ncol(i)-1,am(i)
   str2 = TRIM(str2)//TRIM(str)
  end if
 end do ! for i
end subroutine gen_contracted_string

! read Alpha or (both Alpha and Beta) MOs from a GAMESS .dat or .inp file
! Note: if you want to read both Alpha and Beta MOs, just double the variable
! nif
subroutine read_mo_from_dat(datname, nbf, nif, coeff)
 implicit none
 integer i, j, k, nline, nleft, fid
 integer, intent(in) :: nbf, nif
 character(len=5) :: str1
 character(len=30) :: str2
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
 real(kind=8), intent(out) :: coeff(nbf,nif)

 coeff = 0d0
 open(newunit=fid,file=TRIM(datname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  call upper(buf(3:5))
  if(buf(2:5) == '$VEC') exit
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mo_from_dat: No '$VEC' section in&
                 & file "//TRIM(datname)
  close(fid)
  stop
 end if

 nline = nbf/5
 nleft = nbf - nline*5

 do i = 1, nif, 1
  k = 1
  do j = 1, nline, 1
   read(fid,'(A)') buf
   buf = buf(6:)
   read(buf,'(5ES15.8)') coeff(k:k+4,i)
   k = k + 5
  end do ! for j

  if(nleft > 0) then
   read(fid,'(A)') buf
   buf = buf(6:)
   str1 = ' '
   write(str1,'(I5)') nleft
   str1 = ADJUSTL(str1)
   str2 = '('//TRIM(str1)//'ES15.8)'
   read(buf,TRIM(str2)) coeff(k:nbf,i)
  end if
 end do ! for i

 close(fid)
end subroutine read_mo_from_dat

! read the number of GVB pairs from a GAMESS .dat file
subroutine read_npair_from_dat(datname, npair)
 implicit none
 integer :: i, fid
 integer, intent(out) :: npair
!f2py intent(out) :: npair
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
!f2py intent(in) :: datname

 npair = 0
 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:5) == '$SCF') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_npair_from_dat: no '$SCF' found in&
                   & file "//TRIM(datname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 do while(.true.)
  read(fid,'(A)') buf
  i = INDEX(buf, 'CICOEF')
  if(i > 0) npair = npair + 1
  if(index(buf,'END')>0 .or. index(buf,'VEC')>0) exit
 end do ! for while

 close(fid)
end subroutine read_npair_from_dat

! read CI coefficients from a GAMESS .dat or .inp file
subroutine read_ci_coeff_from_dat(fname, npair, coeff)
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: npair
!f2py intent(in) :: npair
 character(len=240) :: buf
 character(len=240), intent(in) :: fname
!f2py intent(in) :: fname
 real(kind=8), intent(out) :: coeff(2,npair)
!f2py intent(out) :: coeff
!f2py depend(npair) :: coeff

 buf = ' '; coeff = 0d0
 open(newunit=fid,file=TRIM(fname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  j = INDEX(buf,'CICOEF(')
  if(j == 0) j = INDEX(buf,'cicoef(')
  if(j /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ci_coeff_from_dat: no GVB CI coeff&
                   &icients found in'
  write(6,'(A)') 'file '//TRIM(fname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 do i = 1, npair, 1
  read(fid,'(A)') buf
  j = INDEX(buf,'=')
  k = INDEX(buf,',')
  read(buf(j+1:k-1),*) coeff(1,i)
  read(buf(k+1:),*) coeff(2,i)
 end do ! for i

 close(fid)
end subroutine read_ci_coeff_from_dat

! print MOs into .dat file
! if replace is .true., the MOs in the original file will be replaced
! if replace is .false., a new *_new.dat file will be generated
subroutine write_mo_into_dat(datname, nbf, nif, coeff, replace)
 implicit none
 integer :: i, j, k, nline, nleft, fid1, fid2, RENAME
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif)
 character(len=240), intent(in) :: datname
 character(len=240) :: newdat, buf
 logical, intent(in) :: replace

 buf = ' '; newdat = ' '
 i = INDEX(datname,'.dat',.true.)
 if(i == 0) i = INDEX(datname,'.inp',.true.)
 newdat = datname(1:i-1)//'_new.dat'

 open(newunit=fid1,file=TRIM(datname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(newdat),status='replace')
 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  call upper(buf(3:5))
  if(buf(2:5)=='$VEC') exit
 end do ! for while

 ! print MOs
 nline = nbf/5
 nleft = nbf - 5*nline

 do i = 1, nif, 1
  k = MOD(i,100)
  do j = 1, nline, 1
   write(fid2,'(I2,I3,5ES15.8)') k, MOD(j,1000), coeff(5*j-4:5*j,i)
  end do ! for j
  if(nleft > 0) then
   write(fid2,'(I2,I3,5ES15.8)') k, MOD(j,1000), coeff(5*j-4:nbf,i)
  end if
 end do ! for i
 write(fid2,'(A)') ' $END'
 ! print MOs done

 ! skip the MOs in datname
 do while(.true.)
  read(fid1,'(A)') buf
  call upper(buf(3:5))
  if(buf(2:5) == '$END') exit
 end do

 ! copy remaining contens
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid2)
 ! copy done

 if(replace) then
  close(fid1,status='delete')
  i = RENAME(TRIM(newdat), TRIM(datname))
 else
  close(fid1)
 end if
end subroutine write_mo_into_dat

! delete the first $VEC section in a specified .dat file
subroutine del_vec_in_dat(datname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, datname1
 character(len=240), intent(in) :: datname

 datname1 = TRIM(datname)//'.t'
 open(newunit=fid,file=TRIM(datname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(datname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:5) == '$VEC') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:5) == '$END') exit
 end do ! for while

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(datname1), TRIM(datname))
end subroutine del_vec_in_dat

! read the number of alpha/beta electrons from a GAMESS .inp file
subroutine read_na_and_nb_from_gms_inp(inpname, na, nb)
 implicit none
 integer :: i, fid
 integer, intent(out) :: na, nb
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 na = 0; nb = 0
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'GAMESS inp') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_nb_from_gms_inp: failed to read nb&
                   & from file '//TRIM(inpname)
  stop
 end if

 i = INDEX(buf, 'na=')
 read(buf(i+3:),*) na

 i = INDEX(buf, 'nb=')
 read(buf(i+3:),*) nb
end subroutine read_na_and_nb_from_gms_inp

! read ncontr from a GAMESS .dat/.inp file
subroutine read_ncontr_from_gms_inp(inpname, ncontr)
 implicit none
 integer :: i, k, fid, nline
 integer, intent(out) :: ncontr
 character(len=1) :: str1
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 ncontr = 0
 call goto_data_section_in_gms_inp(inpname, fid)
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do while(.true.)
  read(fid,'(A)',iostat=k) buf ! elem(i), nuc(i), coor(:,i)
  if(k /= 0) exit
  if(buf(2:2) == '$') exit

  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
   read(buf,*) str1, nline
   if(str1 == 'L') then
    ncontr = ncontr + 2
   else
    ncontr = ncontr + 1
   end if
   do i = 1, nline, 1
    read(fid,'(A)') buf
   end do ! for i
  end do ! for while
 end do ! for while

 close(fid)
 if(k /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ncontr_from_gms_inp: problematic f&
                   &ile '//TRIM(inpname)
  stop
 end if
end subroutine read_ncontr_from_gms_inp

! read arrays shell_type and shell2atom_map from a GAMESS .dat/.inp file
subroutine read_shltyp_and_shl2atm_from_gms_inp(inpname, ncontr,shltyp,shl2atm)
 implicit none
 integer :: i, j, k, nline, natom, fid
 integer, intent(in) :: ncontr
 integer, intent(out) :: shltyp(ncontr), shl2atm(ncontr)
 character(len=1) :: stype
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 j = 0; natom = 0; shltyp = 0; shl2atm = 0
 call goto_data_section_in_gms_inp(inpname, fid)
 read(fid,'(A)') buf
 read(fid,'(A)') buf

 do while(.true.)
  read(fid,'(A)',iostat=k) buf ! elem(i), nuc(i), coor(:,i)
  if(k /= 0) exit
  if(buf(2:2) == '$') exit
  natom = natom + 1

  do while(.true.)
   read(fid,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
   read(buf,*) stype, nline
   j = j + 1
   call stype2itype(stype, shltyp(j))
   shltyp(j) = shltyp(j) - 1 ! remember to minus 1
   if(shltyp(j) == -1) then ! 'L' -> 'S' & 'P'
    shltyp(j:j+1) = [0,1]
    shl2atm(j:j+1) = [natom,natom]
    j = j + 1
   else
    shl2atm(j) = natom
   end if
   do i = 1, nline, 1
    read(fid,'(A)') buf
   end do ! for i
  end do ! for while
 end do ! for while

 close(fid)
 if(k /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_shltyp_and_shl2atm_from_gms_inp: p&
                   &roblematic'
  write(6,'(A)') 'file '//TRIM(inpname)
  stop
 end if
end subroutine read_shltyp_and_shl2atm_from_gms_inp

! find LenNCZ from the type array all_ecp(natom)
subroutine find_LenNCZ_in_all_ecp(natom, all_ecp, LenNCZ)
 use pg, only: ecp4atom
 implicit none
 integer :: i
 integer, intent(in) :: natom
 integer, intent(out) :: LenNCZ
 type(ecp4atom), intent(in) :: all_ecp(natom)

 LenNCZ = 0

 do i = 1, natom, 1
  if(.not. all_ecp(i)%ecp) cycle
  LenNCZ = LenNCZ + SUM(all_ecp(i)%potential(:)%n)
 end do ! for i
end subroutine find_LenNCZ_in_all_ecp

subroutine all_ecp2ecp_arrays(all_ecp, natom, LenNCZ, KFirst, KLast, Lmax, &
                              LPSkip, NLP, RNFroz, CLP, ZLP)
 use pg, only: ecp4atom
 implicit none
 integer :: i, j, k, m, n
 integer, intent(in) :: natom, LenNCZ
 integer, intent(out) :: KFirst(natom,10), KLast(natom,10), Lmax(natom), &
  LPSkip(natom), NLP(LenNCZ)
 real(kind=8), intent(out) :: RNFroz(natom), CLP(LenNCZ), ZLP(LenNCZ)
 type(ecp4atom), intent(in) :: all_ecp(natom)

 KFirst = 0; KLast = 0; Lmax = 0; LPSkip = 1; NLP = 0
 RNFroz = 0d0; CLP = 0d0; ZLP = 0d0

 forall(i=1:natom, all_ecp(i)%ecp)
  LPSkip(i) = 0
  RNFroz(i) = DBLE(all_ecp(i)%core_e)
  Lmax(i) = all_ecp(i)%highest
 end forall

 m = 0
 do i = 1, natom, 1
  if(.not. all_ecp(i)%ecp) cycle
  k = all_ecp(i)%highest + 1
  do j = 1, k, 1
   n = all_ecp(i)%potential(j)%n
   KFirst(i,j) = m + 1
   KLast(i,j) = m + n
   CLP(m+1:m+n) = all_ecp(i)%potential(j)%col1
   NLP(m+1:m+n) = all_ecp(i)%potential(j)%col2
   ZLP(m+1:m+n) = all_ecp(i)%potential(j)%col3
   m = m + n
  end do ! for j
 end do ! for i
end subroutine all_ecp2ecp_arrays

! read Mulliken charges from a given GAMESS .dat file
subroutine read_mul_char_from_dat(datname, natom, charge, alive)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: natom
 real(kind=8) :: r
 real(kind=8), intent(out) :: charge(natom)
 character(len=2) :: elem
 character(len=240) :: buf
 character(len=240), intent(in) :: datname
 logical, intent(out) :: alive
 ! whether the Mulliken charges exists in .dat file

 charge = 0d0; alive = .false.
 open(newunit=fid,file=TRIM(datname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(2:2) == '$') exit
 end do ! for while

 if(buf(2:2) == '$') then
  close(fid)
  return
 end if

 alive = .true.
 do i = 1, natom, 1
  read(fid,*,iostat=k) elem, r, charge(i)
  if(k /= 0) exit
 end do ! for i

 close(fid)
 if(k /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mul_char_from_dat: no enough charg&
                   &e data to read!'
  write(6,'(A)') 'File='//TRIM(datname)
  stop
 end if
end subroutine read_mul_char_from_dat

! read the relativistic setting from a GAMESS .inp file
subroutine read_irel_from_gms_inp(inpname, irel)
 implicit none
 integer :: i, fid
 integer, intent(out) :: irel
 ! -3/-2/-1/0/2/4 for sfX2C/RESC/None/DKH0/DKH2/DKH4
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 character(len=480) :: buf0

 irel = -1
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:8))
   if(buf(3:8) == 'CONTRL') exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'Warning from subroutine read_irel_from_gms_inp: no $CONTRL &
                   &secton found'
  write(6,'(A)') 'in file '//TRIM(inpname)//'. Treating as non-relativistic.'
  close(fid)
  return
 end if

 i = INDEX(buf, '$END')
 if(i == 0) then
  read(fid,'(A)') buf0
  buf0 = TRIM(buf)//' '//TRIM(buf0)
 else
  buf0 = TRIM(buf)
 end if

 call upper(buf0)
 i = INDEX(buf0, 'RELWFN=DK')
 if(i > 0) irel = 2
 i = INDEX(buf0, 'RELWFN=RESC')
 if(i > 0) irel = -2
 i = INDEX(buf0, 'X2C')
 if(i > 0) irel = -3

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:2) == '$') then
   call upper(buf(3:8))
   if(buf(3:8) == 'RELWFN') exit
  end if
 end do ! for while

 close(fid)

 if(i == 0) then
  call upper(buf)
  i = INDEX(buf, 'NORDER=')
  if(i > 0) then
   if(irel < -1) then
    write(6,'(/,A)') 'ERROR in subroutine read_irel_from_gms_inp: relativistic &
                     &settings are'
    write(6,'(A)') 'inconsistent in file '//TRIM(inpname)
    stop
   end if
   read(buf(i+7:),*) irel
   if(irel == 1) irel = 0 ! DKH0
  end if
 end if
end subroutine read_irel_from_gms_inp

