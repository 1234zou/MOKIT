! written by jxzou at 20220524: move some subroutines from rwwfn.f90 and
!  pop.f90 to this file

module population
 implicit none
 integer :: ncontr, nbf, nif, natom, nmo, i1, i2
! nmo <= nif, the number of MOs to be analyzed. nmo = i2-i1+1
! i1: the beginning index of MOs to be analyzed
! i2: the final index of MOs to be analyzed
 integer, allocatable :: shl2atm(:), shltyp(:), bfirst(:), mo_center(:,:)
! shl2atm: Shell to atom map. Atom begins from 1, not 0
! shltyp: Shell types
!   Spherical     |     Cartesian
! -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
!  H  G  F  D  L  S  P  D  F  G  H
! bfirst: the beginning index of basis func. of each atom, size natom+1
! mo_center: the center(s) of each MOs (multiple centers allowed), (0:natom,nmo)
! mo_center(0,j) is number of atomic centers of the j-th MO
! mo_center(1:,j) is the atomic centers of the j-th MO

 real(kind=8), allocatable :: mo(:,:) ! molecular orbital coefficients
 real(kind=8), allocatable :: ao_ovlp(:,:) ! AO overlap integral matrix
 real(kind=8), allocatable :: mo_dis(:,:) ! distances between MOs, defined as
 ! the shortest distances of two atomic centers
 logical :: cart

 type :: mo_cluster    ! an MO cluster
  integer :: nocc = 0  ! number of occupied MOs
  integer :: nvir = 0  ! number of virtual MOs
  integer, allocatable :: occ_idx(:) ! indices of occupied MOs, size nocc
  integer, allocatable :: vir_idx(:) ! indices of virtual MOs, size nvir
 end type mo_cluster

contains

! initialize arrays shell_type, shell2atom and bfirst
subroutine init_shltyp_shl2atm_bfirst(fchname)
 implicit none
 character(len=240), intent(in) :: fchname

 call read_ncontr_from_fch(fchname, ncontr)
 allocate(shltyp(ncontr), shl2atm(ncontr))
 call read_shltyp_and_shl2atm_from_fch(fchname, ncontr, shltyp, shl2atm)

 natom = shl2atm(ncontr)
 if(allocated(bfirst)) deallocate(bfirst)
 allocate(bfirst(natom+1))
 call get_bfirst(ncontr, shltyp, shl2atm, natom, bfirst)
 deallocate(shltyp, shl2atm)
end subroutine init_shltyp_shl2atm_bfirst

! find centers of each MOs from a specified .fch(k) file
subroutine get_mo_center_from_fch(fchname, ibegin, iend, popm)
 implicit none
 integer :: k
 integer, intent(in) :: ibegin, iend
 real(kind=8), allocatable :: pop(:,:)
 character(len=8) :: pop_meth
 character(len=8), optional :: popm
 character(len=240), intent(in) :: fchname

 if(PRESENT(popm)) then
  pop_meth = popm
 else
  pop_meth = 'lowdin'
 end if

 ! get integer array bfirst (natom would be initialized in this subroutine)
 call init_shltyp_shl2atm_bfirst(fchname)

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 if(allocated(mo)) deallocate(mo)
 if(allocated(ao_ovlp)) deallocate(ao_ovlp)
 allocate(mo(nbf,nif), ao_ovlp(nbf,nbf))

 ! read MOs and AO-basis overlap integrals, respectively
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)
 call get_ao_ovlp_using_fch(fchname, nbf, ao_ovlp)

 if(allocated(mo_center)) deallocate(mo_center)
 allocate(mo_center(0:4,ibegin:iend))
 ! mo_center(0,i) is the number of centers of the i-th MO

 k = iend - ibegin + 1
 allocate(pop(natom,k))
 call calc_diag_gross_pop(natom, nbf, k, bfirst, ao_ovlp, mo(:,ibegin:iend), &
                          pop_meth, pop)
 call get_mo_center_from_pop(natom, k, pop, mo_center)

 deallocate(pop)
 ! arrays mo and ao_ovlp are not deallocated here, since they may be re-used by
 ! other subroutines in this file
end subroutine get_mo_center_from_fch

! calculate distances among MOs using MO centers
subroutine gen_mo_dis_from_mo_center(fchname)
 implicit none
 integer :: i, j, k, m, p, q, i3, i4
 real(kind=8) :: dis0, coor0(3), coor1(3)
 real(kind=8), allocatable :: coor(:,:), dis(:,:)
 character(len=240), intent(in) :: fchname

 allocate(coor(3,natom))
 call read_coor_from_fch(fchname, natom, coor)
 allocate(dis(natom,natom), source=0d0)

 do i = 1, natom-1, 1
  coor0 = coor(:,i)
  do j = i+1, natom, 1
   coor1 = coor(:,j) - coor0
   dis(j,i) = DSQRT(DOT_PRODUCT(coor1, coor1))
   dis(i,j) = dis(j,i)
  end do ! for j
 end do ! for i

 deallocate(coor)
 allocate(mo_dis(i1:i2,i1:i2), source=0d0)

 do i = i1, i2-1, 1
  i3 = mo_center(0,i)
  do j = i+1, i2, 1
   i4 = mo_center(0,j)
   dis0 = dis(mo_center(1,i), mo_center(1,j))

   do k = 1, i3, 1
    p = mo_center(k,i)
    do m = 1, i4, 1
     q = mo_center(m,j)
     if(dis(q,p) < dis0) dis0 = dis(q,p)
    end do ! for m
   end do ! for k

   mo_dis(j,i) = dis0
   mo_dis(i,j) = mo_dis(j,i)
  end do ! for j
 end do ! for i

 deallocate(dis, mo_center)
end subroutine gen_mo_dis_from_mo_center

! perform Mulliken population for a set of MOs, get their MO centers and
! calculate MO distances using their centers. An MO is allowed to have multiple
! centers since we may deal with diradical orbitals.
subroutine get_mo_dis_from_fch(fchname, ibegin, iend)
 implicit none
 integer, intent(in) :: ibegin, iend
 character(len=240), intent(in) :: fchname

 i1 = ibegin; i2 = iend
 if(i2<i1 .or. i1<1 .or. i2<1) then
  write(6,'(/,A)') 'ERROR in subroutine get_mo_dis_from_fch: invalid i1, i2.'
  write(6,'(2(A,I0))') 'i1=', i1, ', i2=', i2 
  stop
 end if

 call get_mo_center_from_fch(fchname, i1, i2)
 call gen_mo_dis_from_mo_center(fchname)
end subroutine get_mo_dis_from_fch

end module population

! delete irrelevant truncated PAOs
subroutine del_irrel_tr_pao(tr_nbf, nbf, s_v, s_d, cpao, tr_cpao, nleft)
 implicit none
 integer :: i, k
 integer, intent(in) :: tr_nbf, nbf
 integer, intent(out) :: nleft
 real(kind=8) :: r
 real(kind=8), parameter :: thres = 0.08d0
 real(kind=8), intent(in) :: s_v(tr_nbf,nbf), s_d(tr_nbf,tr_nbf), cpao(nbf,tr_nbf)
 real(kind=8), intent(inout) :: tr_cpao(tr_nbf,tr_nbf)
 real(kind=8), external :: calc_CiTSCi, calc_CiTSCip
 real(kind=8), allocatable :: new_mo(:,:), mo_i(:)
 logical, allocatable :: del(:)

 allocate(mo_i(tr_nbf), del(tr_nbf))
 del = .false.

 do i = 1, tr_nbf, 1
  mo_i = tr_cpao(:,i)
  r = 1d0 + calc_CiTSCi(tr_nbf,mo_i,s_d) - 2d0*calc_CiTSCip(tr_nbf,nbf,mo_i,s_v,cpao(:,i))
  if(r > thres) del(i) = .true.
 end do ! for j

 deallocate(mo_i)
 k = COUNT(del .eqv. .true.)
 nleft = tr_nbf - k

 if(k == 0) then
  deallocate(del)
  return
 else if(k == tr_nbf) then
  write(6,'(/,A)') 'ERROR in subroutine del_irrel_tr_pao: all truncated PAOs ar&
                   &e irrelevant.'
  write(6,'(A)') 'Maybe because the basis set is too small such that virtual sp&
                 &ace is of poor'
  write(6,'(A)') 'quality. You can use a larger basis set like 6-31G(d,p) or cc&
                 &-pVDZ.'
  stop
 end if

 allocate(new_mo(tr_nbf,nleft))
 k = 0

 do i = 1, tr_nbf, 1
  if(del(i)) cycle
  k = k + 1
  new_mo(:,k) = tr_cpao(:,i)
 end do ! for i
 deallocate(del)

 tr_cpao(:,1:nleft) = new_mo
 deallocate(new_mo)
 if(nleft < tr_nbf) tr_cpao(:,nleft+1:tr_nbf) = 0d0
end subroutine del_irrel_tr_pao

! find an MO which resembles a given MO from a set of MOs
subroutine svd_get_one_mo(nbf, nif, ao_ovlp, mo_k, mo)
 implicit none
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf), mo_k(nbf)
 real(kind=8), intent(inout) :: mo(nbf,nif)
 real(kind=8), allocatable :: mo_k1(:,:), mo_ovlp(:,:), u(:,:), vt(:,:), sv(:),&
  new_mo(:,:)

 allocate(mo_k1(nbf,1))
 mo_k1(:,1) = mo_k

 allocate(mo_ovlp(nif,1))
 call calc_CTSCp2(nbf, nif, 1, mo, ao_ovlp, mo_k1, mo_ovlp)

 allocate(u(nif,nif), vt(1,1), sv(nif))
 call do_svd(nif, 1, mo_ovlp, u, vt, sv)
 write(6,'(A,F20.10)') 'sv=', sv(1)
 deallocate(mo_ovlp, vt, sv)

 allocate(new_mo(nbf,nif), source=0d0)
 call dgemm('N','N', nbf,nif,nif, 1d0,mo,nbf, u,nif, 0d0,new_mo,nbf)
 mo = new_mo
 deallocate(new_mo)
end subroutine svd_get_one_mo

subroutine find_antibonding_orb(fchname, i1, i2, i3)
 use population, only: nbf, nif, bfirst, mo, ao_ovlp, mo_center, &
  get_mo_center_from_fch
 implicit none
 integer :: i, j, k, k1, k2, m, p, tr_nbf, npair, ncenter, nzero
 integer, intent(in) :: i1, i2, i3
!f2py intent(in) :: i1, i2, i3
 integer, allocatable :: chosen(:)
 real(kind=8), allocatable :: ao_dip(:,:,:), dip(:,:,:), mo_dip(:,:,:), &
  pao(:,:), cpao(:,:), tr_cpao(:,:), tr_ocpao(:,:), s_d(:,:), s_v(:,:), &
  svc(:,:), svc1(:)
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call find_specified_suffix(fchname, '.fch', i)
 new_fch = fchname(1:i-1)//'_a.fch'

 ! Find centers of each MO. Arrays bfirst, mo, ao_ovlp are allocated in memory.
 call get_mo_center_from_fch(fchname, i1, i2)

 npair = i2 - i1 + 1
 if(i2<i1 .or. i3<=i2 .or. (nif < i3+npair-1)) then
  write(6,'(/,A)') 'ERROR in subroutine find_antibonding_orb: wrong orbital ind&
                   &ices.'
  write(6,'(A)') 'fchname='//TRIM(fchname)
  write(6,'(3(A,I0))') 'i1=', i1, ', i2=', i2, ', i3=', i3
  stop
 end if

 ! get AO dipole integral matrix
 allocate(ao_dip(nbf,nbf,3))
 call get_gau_ao_dip_from_pyscf(fchname, nbf, ao_dip)
 call check_orthonormal(nbf, nif, mo, ao_ovlp)

 do i = i2, i1, -1
  ! construct PAOs using remaining virtual MOs
  k = i3 + i2 - i - 1
  allocate(pao(nbf,nbf))
  call get_normalized_pao(nbf, k, ao_ovlp, mo(:,1:k), pao)

  ! get the number of PAOs located at the centers of the i-th MO
  tr_nbf = 0
  ncenter = mo_center(0,i)
  do j = 1, ncenter, 1
   k = mo_center(j,i)
   tr_nbf = tr_nbf + bfirst(k+1) - bfirst(k)
  end do ! for j

  allocate(cpao(nbf,tr_nbf), s_v(tr_nbf,nbf), chosen(tr_nbf))
  tr_nbf = 0

  ! select PAOs which are located at the centers of the i-th MO
  do j = 1, ncenter, 1
   k = mo_center(j,i); k1 = bfirst(k); k2 = bfirst(k+1)
   m = k2 - k1
   cpao(:,tr_nbf+1:tr_nbf+m) = pao(:,k1:k2-1)
   s_v(tr_nbf+1:tr_nbf+m,:) = ao_ovlp(k1:k2-1,:)
   forall(p = 1:m) chosen(tr_nbf+p) = k1 + p - 1
   tr_nbf = tr_nbf + m
  end do ! for j
  deallocate(pao)

  ! AO overlap integral matrix within the truncated basis range
  allocate(s_d(tr_nbf,tr_nbf))
  forall(j=1:tr_nbf, k=1:tr_nbf) s_d(j,k) = ao_ovlp(chosen(j),chosen(k))

  ! AO dipole integral matrices within the truncated basis range
  allocate(dip(tr_nbf,tr_nbf,3))
  forall(j=1:tr_nbf,k=1:tr_nbf,m=1:3) dip(j,k,m) = ao_dip(chosen(j),chosen(k),m)

  ! calculate (S_v)C
  allocate(svc(tr_nbf,tr_nbf), source=0d0)
  call dgemm('N','N',tr_nbf,tr_nbf,nbf,1d0,s_v,tr_nbf,cpao,nbf,0d0,svc,tr_nbf)

  ! perform Boughton-Pulay projection to truncate the basis range of PAOs
  allocate(tr_cpao(tr_nbf,tr_nbf))
  call solve_multi_lin_eqs(tr_nbf, tr_nbf, s_d, tr_nbf, svc, tr_cpao)
  deallocate(svc)

  ! Delete irrelevant truncated PAOs. There are several PAOs whose truncation
  ! residues > 0.08, which is possible since not all PAOs can be truncated into
  ! this small basis range. Deleting them is necessary for cases e.g. two
  ! monomers which are far away from each other.
  call del_irrel_tr_pao(tr_nbf, nbf, s_v, s_d, cpao, tr_cpao, k)
  deallocate(cpao)

  ! Perform canonical orthogonalization on truncated PAOs.
  ! Q: Why not simply perform canonical orthogonalization on selected original
  !    PAOs, but on truncated PAOs?
  ! A: Canonical orthogonalization makes MOs extremely delocalized among its
  !    basis range. The basis range of the selected original PAOs are the whole
  !    molecule. So we need to firstly truncate selected PAOs to shrink their
  !    basis range. Then the canonical-orthogonalized PAOs are only delocalized
  !    on centers of the i-th MO (around 2 atoms usually).
  allocate(tr_ocpao(tr_nbf,0:k))
  call orthonormalize_orb(.false., .true., tr_nbf, k, s_d, tr_cpao(:,1:k), &
                          tr_ocpao(:,1:k))
  deallocate(tr_cpao)
  call detect_zero_mo(tr_nbf, k, tr_ocpao(:,1:k), nzero)
  k = k - nzero + 1

  ! perform Boughton-Pulay projection to truncate the i-th MO
  allocate(svc1(tr_nbf), source=0d0)
  call dgemv('N', tr_nbf, nbf, 1d0, s_v, tr_nbf, mo(:,i), 1, 0d0, svc1, 1)
  deallocate(s_v)
  call solve_lin_eqs(tr_nbf, tr_nbf, s_d, svc1, tr_ocpao(:,0))
  deallocate(svc1)
  call normalize_mo(tr_nbf, s_d, tr_ocpao(:,0))
  deallocate(s_d)
  ! tr_ocpao(:,0) is non-orthogonal to tr_ocpao(:,1:k), although the orthogonality
  ! is not required here

  ! calculate MO dipole integral matrices within this basis range
  allocate(mo_dip(3,k,k))
  do j = 1, 3
   call calc_CTSC(tr_nbf, k, tr_ocpao(:,0:k-1), dip(:,:,j), mo_dip(j,:,:))
  end do ! for j
  deallocate(dip)

  ! find the antibonding orbital via 2*2 orbital rotations
  allocate(tr_cpao(tr_nbf,k))
  call assoc_loc2(tr_nbf, k, 0,1, 1,k, tr_ocpao(:,0:k-1), mo_dip, tr_cpao)
  deallocate(tr_ocpao, mo_dip)

  ! map this antibonding orbital back to the basis range of the whole molecule
  allocate(svc1(nbf), source=0d0)
  forall(j = 1:tr_nbf) svc1(chosen(j)) = tr_cpao(j,2)
  deallocate(chosen, tr_cpao)
  call normalize_mo(nbf, ao_ovlp, svc1)

  ! This antibonding orbital is not orthogonal to mo(:,1:i3+i2-i-1). We can get
  !  an eligible virtual MO which resembles this antibonding orbital via SVD
  !  between remaining virtual MOs and this (only one) antibonding orbital.
  ! Note: if this step is not done within the loop, but firstly collect all
  !  antibonding orbitals and try to get orthogonal ones out of the loop. It
  !  would be not easy to obtain orthogonal antibonding orbitals via simple
  !  SVD, and more steps are needed.
  k = i3 + i2 - i
  m = nif - k + 1
  call svd_get_one_mo(nbf, m, ao_ovlp, svc1, mo(:,k:nif))
  deallocate(svc1)
 end do ! for i

 deallocate(bfirst, mo_center, ao_dip)
 call check_orthonormal(nbf, nif, mo, ao_ovlp)
 deallocate(ao_ovlp)

 call sys_copy_file(fchname, new_fch, .false.)
 call write_mo_into_fch(new_fch, nbf, nif, 'a', mo)
end subroutine find_antibonding_orb

! construct antibonding orbitals for a set of bonding orbitals
! Note:
! 1) mo(:,i3:nif) will be updated, where mo(:,i3:i3+i2-i1) are generated anti-
!  bonding orbitals, and mo(:,i3+i2-i1+1:nif) are remaining virtual orbitals.
! 2) all MOs are still orthonormalized.
subroutine find_antibonding_orbitals(i1, i2, i3, natom, nbf, nif, bfirst, &
                                     mo_center, ao_ovlp, ao_dip, mo)
 implicit none
 integer :: i, j, k, k1, k2, m, p, tr_nbf, ncenter, nzero
 integer, intent(in) :: i1, i2, i3, natom, nbf, nif
!f2py intent(in) :: i1, i2, i3, natom, nbf, nif
 integer, intent(in) :: bfirst(natom+1)
!f2py intent(in) :: bfirst
!f2py depend(natom) :: bfirst
 real(kind=8), intent(in) :: mo_center(0:4,i1:i2), ao_ovlp(nbf,nbf), &
  ao_dip(3,nbf,nbf)
!f2py intent(in) :: mo_center, ao_ovlp, ao_dip
!f2py depend(i1,i2) :: mo_center
!f2py depend(nbf) :: ao_ovlp, ao_dip
 integer, allocatable :: chosen(:)
 real(kind=8), allocatable :: dip(:,:,:), mo_dip(:,:,:), pao(:,:), cpao(:,:), &
  tr_cpao(:,:), tr_ocpao(:,:), s_d(:,:), s_v(:,:), svc(:,:), svc1(:)
 real(kind=8), intent(inout) :: mo(nbf,nif)
!f2py intent(in,out) :: mo
!f2py depend(nbf,nif) :: mo

 do i = i2, i1, -1
  ! construct PAOs using remaining virtual MOs
  k = i3 + i2 - i - 1
  allocate(pao(nbf,nbf))
  call get_normalized_pao(nbf, k, ao_ovlp, mo(:,1:k), pao)

  ! get the number of PAOs located at the centers of the i-th MO
  tr_nbf = 0
  ncenter = mo_center(0,i)
  do j = 1, ncenter, 1
   k = mo_center(j,i)
   tr_nbf = tr_nbf + bfirst(k+1) - bfirst(k)
  end do ! for j

  allocate(cpao(nbf,tr_nbf), s_v(tr_nbf,nbf), chosen(tr_nbf))
  tr_nbf = 0

  ! select PAOs which are located at the centers of the i-th MO
  do j = 1, ncenter, 1
   k = mo_center(j,i); k1 = bfirst(k); k2 = bfirst(k+1)
   m = k2 - k1
   cpao(:,tr_nbf+1:tr_nbf+m) = pao(:,k1:k2-1)
   s_v(tr_nbf+1:tr_nbf+m,:) = ao_ovlp(k1:k2-1,:)
   forall(p = 1:m) chosen(tr_nbf+p) = k1 + p - 1
   tr_nbf = tr_nbf + m
  end do ! for j
  deallocate(pao)

  ! AO overlap integral matrix within the truncated basis range
  ! AO dipole integral matrices within the truncated basis range
  allocate(s_d(tr_nbf,tr_nbf), dip(3,tr_nbf,tr_nbf))
!$omp parallel do schedule(dynamic) default(shared) private(j,k)
  do j = 1, tr_nbf, 1
   do k = 1, tr_nbf, 1
    s_d(k,j) = ao_ovlp(chosen(k),chosen(j))
    dip(:,k,j) = ao_dip(:,chosen(k),chosen(j))
   end do ! for k
  end do ! for j
!$omp end parallel do

  ! calculate (S_v)C
  allocate(svc(tr_nbf,tr_nbf), source=0d0)
  call dgemm('N','N',tr_nbf,tr_nbf,nbf,1d0,s_v,tr_nbf,cpao,nbf,0d0,svc,tr_nbf)

  ! perform Boughton-Pulay projection to truncate the basis range of PAOs
  allocate(tr_cpao(tr_nbf,tr_nbf))
  call solve_multi_lin_eqs(tr_nbf, tr_nbf, s_d, tr_nbf, svc, tr_cpao)
  deallocate(svc)

  ! Delete irrelevant truncated PAOs. There are several PAOs whose truncation
  ! residues > 0.08, which is possible since not all PAOs can be truncated into
  ! this small basis range. Deleting them is necessary for cases e.g. two
  ! monomers which are far away from each other.
  call del_irrel_tr_pao(tr_nbf, nbf, s_v, s_d, cpao, tr_cpao, k)
  deallocate(cpao)

  ! Perform canonical orthogonalization on truncated PAOs.
  ! Q: Why not simply perform canonical orthogonalization on selected original
  !    PAOs, but on truncated PAOs?
  ! A: Canonical orthogonalization makes MOs extremely delocalized among its
  !    basis range. The basis range of the selected original PAOs are the whole
  !    molecule. So we need to firstly truncate selected PAOs to shrink their
  !    basis range. Then the canonical-orthogonalized PAOs are only delocalized
  !    on centers of the i-th MO (around 2 atoms usually).
  allocate(tr_ocpao(tr_nbf,0:k))
  call orthonormalize_orb(.false., .true., tr_nbf, k, s_d, tr_cpao(:,1:k), &
                          tr_ocpao(:,1:k))
  deallocate(tr_cpao)
  call detect_zero_mo(tr_nbf, k, tr_ocpao(:,1:k), nzero)
  k = k - nzero + 1

  ! perform Boughton-Pulay projection to truncate the i-th MO
  allocate(svc1(tr_nbf), source=0d0)
  call dgemv('N', tr_nbf, nbf, 1d0, s_v, tr_nbf, mo(:,i), 1, 0d0, svc1, 1)
  deallocate(s_v)
  call solve_lin_eqs(tr_nbf, tr_nbf, s_d, svc1, tr_ocpao(:,0))
  deallocate(svc1)
  call normalize_mo(tr_nbf, s_d, tr_ocpao(:,0))
  deallocate(s_d)
  ! tr_ocpao(:,0) is non-orthogonal to tr_ocpao(:,1:k), although the orthogonality
  ! is not required here

  ! calculate MO dipole integral matrices within this basis range
  allocate(mo_dip(3,k,k))
  do j = 1, 3
   call calc_CTSC(tr_nbf, k, tr_ocpao(:,0:k-1), dip(j,:,:), mo_dip(j,:,:))
  end do ! for j
  deallocate(dip)

  ! find the antibonding orbital via 2*2 orbital rotations
  allocate(tr_cpao(tr_nbf,k))
  call assoc_loc2(tr_nbf, k, 0,1, 1,k, tr_ocpao(:,0:k-1), mo_dip, tr_cpao)
  deallocate(tr_ocpao, mo_dip)

  ! map this antibonding orbital back to the basis range of the whole molecule
  allocate(svc1(nbf), source=0d0)
  forall(j = 1:tr_nbf) svc1(chosen(j)) = tr_cpao(j,2)
  deallocate(chosen, tr_cpao)
  call normalize_mo(nbf, ao_ovlp, svc1)

  ! This antibonding orbital is not orthogonal to mo(:,1:i3+i2-i-1). We can get
  !  an eligible virtual MO which resembles this antibonding orbital via SVD
  !  between remaining virtual MOs and this (only one) antibonding orbital.
  ! Note: if this step is not done within the loop, but firstly collect all
  !  antibonding orbitals and try to get orthogonal ones out of the loop. It
  !  would be not easy to obtain orthogonal antibonding orbitals via simple
  !  SVD, and more steps are needed.
  k = i3 + i2 - i
  m = nif - k + 1
  call svd_get_one_mo(nbf, m, ao_ovlp, svc1, mo(:,k:nif))
  deallocate(svc1)
 end do ! for i
end subroutine find_antibonding_orbitals

! calculate GVB bond orders using information in _s.fch and _s.dat files
subroutine get_gvb_bond_order_from_fch(fchname)
 use population, only: mo_center, get_mo_center_from_fch
 implicit none
 integer :: i, j, k, m, i1, j1, npair, na, nb, nopen, ibegin
 integer, allocatable :: idx(:)
 real(kind=8), allocatable :: ci_coeff(:,:), bo(:)
 character(len=20) :: str
 character(len=240) :: datname
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 i = INDEX(fchname, '.fch', back=.true.)
 datname = fchname(1:i-1)//'.dat'
 call read_npair_from_dat(datname, npair)
 if(npair > 0) then
  allocate(ci_coeff(2,npair), bo(npair))
  call read_ci_coeff_from_dat(datname, npair, ci_coeff)
  forall(i = 1:npair)
   bo(i) = DABS(ci_coeff(1,i)*ci_coeff(1,i) - ci_coeff(2,i)*ci_coeff(2,i))
  end forall
  deallocate(ci_coeff)
 end if

 call read_na_and_nb_from_fch(fchname, na, nb)
 nopen = na - nb
 ibegin = nb - npair + 1
 ! Lowdin populations are performed for bonding orbitals only
 call get_mo_center_from_fch(fchname, ibegin, nb)

 do i = ibegin, nb, 1
  k = mo_center(0,i)
  if(k < 2) cycle
  allocate(idx(k))
  call sort_int_array(k, mo_center(1:k,i), .true., idx)
  deallocate(idx)
 end do ! for i

 write(6,'(/,A)') '--- GVB Bond Order Analysis ---'
 write(6,'(A)') 'GVB bond order: (n_{i,1} - n_{i,2})/2, see J. Phys. Chem. A 20&
                &20, 124, 8321.'

 write(6,'(/,A)') 'GVB pairs identified as atomic core orbitals or lone pairs:'
 write(6,'(A)') ' Pair   Center'
 do i = 1, npair, 1
  if(mo_center(0,ibegin+i-1) /= 1) cycle
  write(6,'(2I5)') i, mo_center(1,ibegin+i-1)
 end do ! for i

 write(6,'(/,A)') 'GVB pairs identified as chemical bonds:'
 write(6,'(A)') ' Pair   Bond order  Center'
 do i = 1, npair, 1
  j = ibegin + i - 1
  if(mo_center(0,j) < 2) cycle
  k = mo_center(0,j)
  write(str,'(A,I0,A)') '(I5,3X,F5.3,4X,', k, 'I5)'
  write(6,TRIM(str)) i, bo(i), mo_center(1:k,j)
 end do ! for i

 write(6,'(/,A)') 'Merging chemical bonds which share the same atomic centers:'
 write(6,'(A)') ' Bond order  Center'
 do i = 1, npair, 1
  i1 = ibegin + i - 1
  k = mo_center(0,i1)
  if(k < 2) cycle
  do j = i+1, npair, 1
   j1 = ibegin + j - 1
   m = mo_center(0,j1)
   if(m<2 .or. k/=m) cycle
   if(ALL(mo_center(1:k,i1) == mo_center(1:k,j1))) then
    bo(i) = bo(i) + bo(j); bo(j) = 0d0
    mo_center(:,j1) = 0
   end if
  end do ! for j
  write(str,'(A,I0,A)') '(F7.3,3X,', k, 'I5)'
  write(6,TRIM(str)) bo(i), mo_center(1:k,i1)
 end do ! for i
 if(npair > 0) deallocate(bo, mo_center)

 write(6,'(/,A)') 'Remark: only GVB pair orbitals are analyzed here. Be careful&
                  & with extra bond'
 write(6,'(A)') 'orders coming from the doubly occupied space.'
 write(6,'(/,A)') '--- End of GVB Bond Order Analysis ---'
end subroutine get_gvb_bond_order_from_fch

! get integer array bfirst (the beginning index of basis func. of each atom)
subroutine get_bfirst(ncontr, shltyp, shl2atm, natom, bfirst)
 implicit none
 integer :: i, j
 integer, allocatable :: ang0(:)
 integer, intent(in) :: ncontr, natom
 integer, intent(in) :: shltyp(ncontr), shl2atm(ncontr)
 integer, intent(out) :: bfirst(natom+1)
 logical :: cart

 bfirst = 0; cart = .false.
 if( ANY(shltyp>1) ) cart = .true. ! 6D 10F

 allocate(ang0(ncontr), source=0)
 if(cart) then
  forall(i = 1:ncontr) ang0(i) = (shltyp(i)+1)*(shltyp(i)+2)/2
 else
  where(shltyp == -1)
   ang0 = 4
  elsewhere
   ang0 = 2*IABS(shltyp) + 1
  end where
 end if

 bfirst = 0; bfirst(1) = 1

 do i = 1, ncontr, 1
  j = shl2atm(i) + 1
  bfirst(j) = bfirst(j) + ang0(i)
 end do ! for i

 deallocate(ang0)

 do i = 2, natom+1, 1
  bfirst(i) = bfirst(i) + bfirst(i-1)
 end do ! for i
end subroutine get_bfirst

! get the integer array bfirst from a specified .fch file
subroutine get_bfirst_from_fch(fchname, natom, bfirst)
 implicit none
 integer :: ncontr
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 integer, intent(out) :: bfirst(natom+1)
!f2py intent(out) :: bfirst
!f2py depend(natom) :: bfirst
 integer, allocatable :: shltyp(:) , shl2atm(:)
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call read_ncontr_from_fch(fchname, ncontr)
 allocate(shltyp(ncontr), shl2atm(ncontr))
 call read_shltyp_and_shl2atm_from_fch(fchname, ncontr, shltyp, shl2atm)

 call get_bfirst(ncontr, shltyp, shl2atm, natom, bfirst)
 deallocate(shltyp, shl2atm)
end subroutine get_bfirst_from_fch

! calculate the Mulliken gross populations based on provided density matrix
subroutine mulliken_pop_of_dm(natom, bfirst, nbf, dm, ovlp, gross)
 implicit none
 integer :: i, u
 integer, intent(in) :: natom, nbf
 integer, intent(in) :: bfirst(natom+1)
 ! bfirst: the beginning index of basis func. of each atom
 real(kind=8), intent(in) :: dm(nbf,nbf), ovlp(nbf,nbf)
 real(kind=8) :: rtmp, ddot
 real(kind=8), intent(out) :: gross(natom)

!$omp parallel do schedule(dynamic) default(shared) private(i,u,rtmp)
 do i = 1, natom, 1
  rtmp = 0d0
  do u = bfirst(i), bfirst(i+1)-1, 1
   rtmp = rtmp + ddot(nbf, dm(u,:), 1, ovlp(:,u), 1)
  end do ! for u
  gross(i) = rtmp
 end do ! for i
!$omp end parallel do
end subroutine mulliken_pop_of_dm

! calculate the Lowdin gross populations based on provided density matrix
subroutine lowdin_pop_of_dm(natom, bfirst, nbf, dm, ovlp, is_sqrt, gross)
 implicit none
 integer :: i, u
 integer, intent(in) :: natom, nbf
 integer, intent(in) :: bfirst(natom+1)
 ! bfirst: the beginning index of basis func. of each atom
 real(kind=8) :: rtmp, ddot
 real(kind=8), intent(in) :: dm(nbf,nbf), ovlp(nbf,nbf)
 real(kind=8), intent(out) :: gross(natom)
 real(kind=8), allocatable :: sqrt_s(:,:), n_sqrt_s(:,:)
 logical, intent(in) :: is_sqrt
 ! True : input ovlp is S^(1/2) actually
 ! False: input ovlp is S and we need to calculate S^(1/2)
 ! When this subroutine is frequently called, repeated calculation of S^(1/2)
 ! is time-consuming, so the input of S^(1/2) is allowed. We can calculate
 ! S^(1/2) outside of this subroutine and use it repeatedly.

 if(is_sqrt) then
  allocate(sqrt_s(nbf,nbf), source=ovlp)
  allocate(n_sqrt_s(nbf,nbf), source=0d0)
 else
  allocate(sqrt_s(nbf,nbf), n_sqrt_s(nbf,nbf))
  call mat_dsqrt(nbf, ovlp, .false., sqrt_s, n_sqrt_s)
  n_sqrt_s = 0d0
 end if

 ! use n_sqrt_s to store (S^(1/2))P
 call dsymm('L', 'L', nbf, nbf, 1d0,sqrt_s,nbf, dm,nbf, 0d0,n_sqrt_s,nbf)

!$omp parallel do schedule(dynamic) default(shared) private(i,u,rtmp)
 do i = 1, natom, 1
  rtmp = 0d0
  do u = bfirst(i), bfirst(i+1)-1, 1
   rtmp = rtmp + ddot(nbf, n_sqrt_s(u,:), 1, sqrt_s(:,u), 1)
  end do ! for u
  gross(i) = rtmp
 end do ! for i
!$omp end parallel do
end subroutine lowdin_pop_of_dm

! Rotate one pair of occ/vir orbitals each time, to make the Mulliken(+Lowdin)
! gross populations match the given gross populations.
subroutine rot_mo2match_gross_pop(nocc, nbf, nif, natom, bfirst, mo, ovlp, &
                                  m_gross0, l_gross0, new_mo)
 implicit none
 integer :: i, j, a, u, v, niter
 integer, parameter :: niter_max = 999
 integer, intent(in) :: nocc, nbf, nif, natom
 integer, intent(in) :: bfirst(natom+1)
 real(kind=8) :: decr, r1, r2, H, II, JJ, KK, cos_t, sin_t, y, ddot
 real(kind=8), intent(in) :: mo(nbf,nif), ovlp(nbf,nbf), m_gross0(natom), &
  l_gross0(natom)
 real(kind=8), intent(out) :: new_mo(nbf,nif)
 real(kind=8), parameter :: threshold1 = 1d-8, threshold2 = 1d-6
 real(kind=8), allocatable :: sqrt_s(:,:), n_sqrt_s(:,:), sqrt_s_b(:,:), &
  sqrt_s_d(:,:), dm(:,:), m_gross(:), l_gross(:), B(:,:), D(:,:), E(:), F(:), &
  G(:), ri(:), ra(:)

 new_mo = mo
 allocate(dm(nbf,nbf))
 call calc_cct(nbf, nocc, mo(:,1:nocc), dm)
 dm = 2d0*dm
 allocate(m_gross(natom))
 call mulliken_pop_of_dm(natom, bfirst, nbf, dm, ovlp, m_gross)
 write(6,'(/,A)') 'Initial Mulliken gross populations:'
 write(6,'(5(1X,ES15.8))') m_gross

 allocate(sqrt_s(nbf,nbf), n_sqrt_s(nbf,nbf))
 call mat_dsqrt(nbf, ovlp, .false., sqrt_s, n_sqrt_s)
 deallocate(n_sqrt_s)
 allocate(l_gross(natom))
 call lowdin_pop_of_dm(natom, bfirst, nbf, dm, sqrt_s, .true., l_gross)
 write(6,'(A)') 'Initial Lowdin gross populations:'
 write(6,'(5(1X,ES15.8))') l_gross

 allocate(E(natom))
 E = m_gross - m_gross0 + l_gross - l_gross0
 y = ddot(natom, E, 1, E, 1)
 if(y < threshold2) then
  write(6,'(/,A)') 'Initial gross populations already equal to target ones.'
!  deallocate(dm, gross, E)
!  return
 end if

 allocate(B(nbf,nbf), D(nbf,nbf), sqrt_s_b(nbf,nbf), sqrt_s_d(nbf,nbf), &
          F(natom), G(natom), ri(nbf), ra(nbf))

 ! perform 2*2 rotation
 niter = 0
 do while(niter <= niter_max)
  decr = 0d0

  do i = nocc, 1, -1
   do a = nocc+1, nif, 1
    !$omp parallel do schedule(dynamic) default(shared) private(u,v,r1,r2)
    do v = 1, nbf, 1
     r1 = new_mo(v,a); r2 = new_mo(v,i)
     do u = 1, nbf, 1
      B(u,v) = new_mo(u,a)*r1 - new_mo(u,i)*r2
      D(u,v) = 2d0*(new_mo(u,i)*r1 + new_mo(u,a)*r2)
     end do ! for u
    end do ! for v
    !$omp end parallel do

    sqrt_s_b = 0d0; sqrt_s_d = 0d0
    call dsymm('L', 'L', nbf, nbf, 1d0,sqrt_s,nbf, B,nbf, 0d0,sqrt_s_b,nbf)
    call dsymm('L', 'L', nbf, nbf, 1d0,sqrt_s,nbf, D,nbf, 0d0,sqrt_s_d,nbf)

    do j = 1, natom, 1
     r1 = 0d0; r2 = 0d0
     do u = bfirst(j), bfirst(j+1)-1, 1
      r1 = r1 + ddot(nbf, B(u,:), 1, ovlp(:,u), 1) + &
                ddot(nbf, sqrt_s_b(u,:), 1, sqrt_s(:,u), 1)
      r2 = r2 + ddot(nbf, D(u,:), 1, ovlp(:,u), 1) + &
                ddot(nbf, sqrt_s_d(u,:), 1, sqrt_s(:,u), 1)
     end do ! for u
     E(j) = r1; F(j) = r2
    end do ! for j

    G = E + m_gross - m_gross0 + l_gross - l_gross0
    H = 2d0*ddot(natom, E, 1, G, 1)
    II = -2d0*ddot(natom, F, 1, G, 1)
    JJ = 0.5d0*(ddot(natom, F, 1, F, 1) - ddot(natom, E, 1, E, 1))
    KK = ddot(natom, E, 1, F, 1)
    ! this subroutine finds the maximum of function
    ! f(x) = Acos2x - Bsin2x + Ccos4x - Dsin4x - A - C
    call find_cos_quartic_poly_maximum(H, II, JJ, KK, cos_t, sin_t, y)
    ! if the change of the function value is too tiny, no need to rotate
    if(y < threshold1) cycle

    write(6,'(2I4,3F20.10)') i, a, cos_t, sin_t, y
    decr = decr + y
    ri = new_mo(:,i); ra = new_mo(:,a)
    new_mo(:,i) = cos_t*ri - sin_t*ra
    new_mo(:,a) = sin_t*ri + cos_t*ra
    call calc_cct(nbf, nocc, new_mo(:,1:nocc), dm)
    dm = 2d0*dm
    call mulliken_pop_of_dm(natom, bfirst, nbf, dm, ovlp, m_gross)
    call lowdin_pop_of_dm(natom, bfirst, nbf, dm, sqrt_s, .true., l_gross)
   end do ! for a
  end do ! for i

  niter = niter + 1
  write(6,'(A,I3,A,F13.8)') 'niter=', niter, ', decrease=', decr
  if(decr < threshold2) exit
 end do ! for while

 deallocate(sqrt_s_b, sqrt_s_d, B, D, E, F, G, ri, ra)
 if(niter <= niter_max) then
  write(6,'(A)') 'Rotations converged successfully.'
 else
  write(6,'(A,I5)') 'niter_max=', niter_max
  write(6,'(A)') 'Rotations failed to converge.'
 end if

 call calc_cct(nbf, nocc, new_mo(:,1:nocc), dm)
 dm = 2d0*dm
 call mulliken_pop_of_dm(natom, bfirst, nbf, dm, ovlp, m_gross)
 write(6,'(/,A)') 'Final Mulliken gross populations:'
 write(6,'(5(1X,ES15.8))') m_gross
 call lowdin_pop_of_dm(natom, bfirst, nbf, dm, sqrt_s, .true., l_gross)
 write(6,'(A)') 'Final Lowdin gross populations:'
 write(6,'(5(1X,ES15.8))') l_gross

 deallocate(dm, sqrt_s, m_gross, l_gross)
end subroutine rot_mo2match_gross_pop

! read Mulliken gross populations from a specified .fch file
subroutine read_mull_gross_pop_from_fch(fchname, natom, gross)
 implicit none
 integer :: i, fid
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 real(kind=8), intent(out) :: gross(natom)
!f2py intent(out) :: gross
!f2py depend(natom) :: gross
 real(kind=8), allocatable :: eff_nuc(:)
 character(len=10) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 gross = 0d0
 allocate(eff_nuc(natom), source=0d0)
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'Nuclear ch') exit
 end do ! for i

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mull_gross_pop_from_fch: 'Nuclear &
                   &ch' not found in"
  write(6,'(A)') 'file '//TRIM(fchname)
  close(fid)
  stop
 end if

 read(fid,*) eff_nuc

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'Mulliken C') exit
 end do ! for i

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mull_gross_pop_from_fch: 'Mulliken&
                   & C' not found"
  write(6,'(A)') 'in file '//TRIM(fchname)
  close(fid)
  stop
 end if

 read(fid,*) gross ! Mulliken charges
 close(fid)

 gross = eff_nuc - gross ! effect nuclear charges considered
 ! Be careful about the sign of gross. The electrons possess negative charges.
 ! The gross is usually taken as positive. -gross contain the electric charges
 ! of each atom. eff_nuc-gross contain the atomic partial charges.
 deallocate(eff_nuc)
end subroutine read_mull_gross_pop_from_fch

subroutine solve_mo_from_gross_pop(fchname, natom, m_gross0, l_gross0)
 implicit none
 integer :: i, nbf, nif, na, nb
 integer, allocatable :: bfirst(:)
 integer, intent(in) :: natom
!f2py intent(in) :: natom
 real(kind=8), allocatable :: mo(:,:), new_mo(:,:), ovlp(:,:)
 real(kind=8), intent(in) :: m_gross0(natom), l_gross0(natom)
!f2py intent(in) :: m_gross0, l_gross0
!f2py depend(natom) :: m_gross0, l_gross0
 character(len=240) :: new_fch
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call find_specified_suffix(fchname, '.fch', i)
 new_fch = fchname(1:i-1)//'_new.fch'

 call read_na_and_nb_from_fch(fchname, na, nb)
 if(na /= nb) then
  write(6,'(/,A)') 'ERROR in subroutine solve_mo_from_gross_pop: currently only&
                   & the closed-shell'
  write(6,'(A)') 'system is supported.'
  stop
 end if

 allocate(bfirst(natom+1))
 call get_bfirst_from_fch(fchname, natom, bfirst)

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(mo(nbf,nif))
 call read_mo_from_fch(fchname, nbf, nif, 'a', mo)

 allocate(ovlp(nbf,nbf))
 call get_ao_ovlp_using_fch(fchname, nbf, ovlp)

 write(6,'(A)') 'Target Mulliken gross populations:'
 write(6,'(5(1X,ES15.8))') m_gross0
 write(6,'(A)') 'Target Lowdin gross populations:'
 write(6,'(5(1X,ES15.8))') l_gross0

 allocate(new_mo(nbf,nif))
 call rot_mo2match_gross_pop(na, nbf, nif, natom, bfirst, mo, ovlp, m_gross0, &
                             l_gross0, new_mo)
 deallocate(bfirst, mo, ovlp)

 call sys_copy_file(TRIM(fchname), TRIM(new_fch), .false.)
 call write_mo_into_fch(new_fch, nbf, nif, 'a', new_mo)
 deallocate(new_mo)
end subroutine solve_mo_from_gross_pop
 
! calculate the number of unpaired electrons and generate unpaired electron
! density .fch file
! Note: the input fchname must include natural orbitals and corresponding
! orbital occupation numbers.
subroutine calc_unpaired_from_fch(fchname, wfn_type, gen_dm, unpaired_e)
 implicit none
 integer :: i, j, k, ne, nbf, nif, mult, fid
 integer, intent(in) :: wfn_type ! 1/2/3 for UNO/GVB/CASSCF NOs
!f2py intent(in) :: wfn_type
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 real(kind=8) :: t0, t1, y0, y1, upe(3)
 real(kind=8), intent(out) :: unpaired_e
!f2py intent(out) :: unpaired_e
 real(kind=8), allocatable :: noon(:,:), coeff(:,:), dm(:,:)
 logical, intent(in) :: gen_dm
!f2py intent(in) :: gen_dm
 ! True/False: generate unpaired/odd electron density or not

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(noon(nif,5))
 call read_eigenvalues_from_fch(fchname, nif, 'a', noon(:,1))

 write(6,'(A)') REPEAT('-',23)//' Radical index '//REPEAT('-',23)
 call read_mult_from_fch(fchname, mult)

 if(mult == 1) then
  open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:14) == 'Number of elec') exit
  end do ! for while
  close(fid)
  read(buf(50:),*) ne
  i = ne/2   ! assuming NOs are ordered in decreasing occupation number
  if(wfn_type == 1) then ! UNO
   t0 = (noon(i,1) - noon(i+1,1))*0.5d0
   y0 = 1d0 - 2d0*t0/(1d0+t0*t0)
   t1 = (noon(i-1,1) - noon(i+2,1))*0.5d0
   y1 = 1d0 - 2d0*t1/(1d0+t1*t1)
   write(6,'(A,F7.3)') 'biradical character   (1-2t/(1+t^2)) y0=', y0
   write(6,'(A,F7.3)') 'tetraradical character(1-2t/(1+t^2)) y1=', y1
  else ! GVB/CASSCF NOs
   ! For GVB/CAS NOs, there is no unique way to define radical index.
   ! Here we adopt the occupation numbers of LUNO and LUNO+1.
   ! You can adopt the way of calculating y0/y1 in UHF, if you like.
   y0 = noon(i+1,1); y1 = noon(i+2,1)
   write(6,'(A,F7.3)') 'biradical character y0 = n_LUNO =', y0
   write(6,'(A,F7.3)') 'tetraradical character y1 = n_{LUNO+1} =', y1
  end if
 else
  write(6,'(A)') 'Not spin singlet. Biradical character will not be computed.'
 end if

 call prt_unpaired_e(nif, noon, upe)
 unpaired_e = upe(3)

 if(gen_dm) then
  i = INDEX(fchname, '.fch')
  fchname1 = fchname(1:i-1)//'_unpaired.fch'
  call copy_file(fchname, fchname1, .false.)

  allocate(coeff(nbf,nif))
  call read_mo_from_fch(fchname, nbf, nif, 'a', coeff)
  allocate(dm(nbf,nbf), source=0d0) ! initialization

  do i = 1, nbf, 1
   do j = 1, i, 1
    do k = 1, nif, 1
     if(noon(k,5)<1d-5 .or. (2d0-noon(k,5)<1d-5)) cycle
     dm(j,i) = dm(j,i) + noon(k,5)*coeff(j,k)*coeff(i,k)
    end do ! for k
   end do ! for j
  end do ! for i

  deallocate(coeff)
  call write_dm_into_fch(fchname1, .true., nbf, dm)
  deallocate(dm)
  call write_eigenvalues_to_fch(fchname1, nif, 'a', noon(:,5), .true.)
 end if

 deallocate(noon)
end subroutine calc_unpaired_from_fch

! calculate the number of unpaired electrons using a GAMESS GVB .dat file
subroutine calc_unpaired_from_dat(datname, mult, unpaired_e)
 implicit none
 integer :: nopen, npair, nif
 integer, intent(in) :: mult
!f2py intent(in) ::  mult
 character(len=240), intent(in) :: datname
!f2py intent(in) :: datname
 real(kind=8) :: upe(3)
 real(kind=8), intent(out) :: unpaired_e
!f2py intent(out) :: unpaired_e
 real(kind=8), allocatable :: pair_coeff(:,:), noon(:,:)

 if(mult < 1) then
  write(6,'(/,A)') 'ERROR in subroutine calc_unpaired_from_dat: mult<1.'
  write(6,'(A)') 'Wrong parameter for spin multiplicity!'
  stop
 end if

 nopen = mult - 1
 call read_npair_from_dat(datname, npair)
 allocate(pair_coeff(2,npair))
 call read_ci_coeff_from_dat(datname, npair, pair_coeff)
 
 nif = npair*2 + nopen
 allocate(noon(nif,5), source=0d0)
 call gen_noon_from_pair_coeff(npair, pair_coeff, nopen, nif, noon(:,1))
 deallocate(pair_coeff)

 call prt_unpaired_e(nif, noon, upe)
 deallocate(noon)
 unpaired_e = upe(3)
end subroutine calc_unpaired_from_dat

! calculate the number of unpaired electrons using a GAMESS GVB .gms file
subroutine calc_unpaired_from_gms_out(outname, unpaired_e)
 implicit none
 integer :: ncore, nopen, npair, nif
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname
 real(kind=8) :: upe(3)
 real(kind=8), intent(out) :: unpaired_e
!f2py intent(out) :: unpaired_e
 real(kind=8), allocatable :: pair_coeff(:,:), noon(:,:)

 call read_npair_from_gms(outname, ncore, nopen, npair)
 allocate(pair_coeff(2,npair))
 call read_ci_coeff_from_gms(outname, npair, pair_coeff)
 
 nif = npair*2 + nopen
 allocate(noon(nif,5), source=0d0)
 call gen_noon_from_pair_coeff(npair, pair_coeff, nopen, nif, noon(:,1))
 deallocate(pair_coeff)

 call prt_unpaired_e(nif, noon, upe)
 deallocate(noon)
 unpaired_e = upe(3)
end subroutine calc_unpaired_from_gms_out

! generate GVB NOON from pair coefficients
subroutine gen_noon_from_pair_coeff(npair, pair_coeff, nopen, nif, noon)
 implicit none
 integer :: i
 integer, intent(in) :: npair, nopen, nif
 real(kind=8), intent(in) :: pair_coeff(2,npair)
 real(kind=8), intent(out) :: noon(nif)

 noon = 0d0
 if(nopen > 0) noon(1:nopen) = 1d0
 if(npair == 0) return

 forall(i = 1:npair)
  noon(nopen+2*i-1) = 2d0*pair_coeff(1,i)*pair_coeff(1,i)
  noon(nopen+2*i) = 2d0*pair_coeff(2,i)*pair_coeff(2,i)
 end forall
end subroutine gen_noon_from_pair_coeff

! print unpaired electron information
subroutine prt_unpaired_e(nif, noon, upe)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 real(kind=8), intent(inout) :: noon(nif,5)
 real(kind=8), intent(out) :: upe(3)

 forall(i = 1:nif) noon(i,2) = 2d0 - noon(i,1)
 forall(i = 1:nif)
  noon(i,3) = noon(i,1)*noon(i,2)
  noon(i,4) = MIN(noon(i,1), noon(i,2))
 end forall
 forall(i = 1:nif) noon(i,5) = noon(i,3)*noon(i,3)

 upe(1) = SUM(noon(:,3))
 upe(2) = SUM(noon(:,4))
 upe(3) = SUM(noon(:,5))
 write(6,'(A,F7.3)') "Yamaguchi's unpaired electrons  (sum_n n(2-n)      ):",upe(1)
 write(6,'(A,F7.3)') "Head-Gordon's unpaired electrons(sum_n min(n,(2-n))):",upe(2)
 write(6,'(A,F7.3)') "Head-Gordon's unpaired electrons(sum_n (n(2-n))^2  ):",upe(3)
 write(6,'(A)') REPEAT('-',61)
end subroutine prt_unpaired_e

! read CI coefficients from a GAMESS .gms file
! Note: if there exist multiple sets of CI coefficients in the file,
!       only the last set will be read
subroutine read_ci_coeff_from_gms(fname, npair, coeff)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: npair
!f2py intent(in) :: npair
 character(len=240) :: buf
 character(len=240), intent(in) :: fname
!f2py intent(in) :: fname
 real(kind=8), intent(out) :: coeff(2,npair)
!f2py intent(out) :: coeff
!f2py depend(npair) :: coeff

 buf = ' '; coeff = 0d0
 open(newunit=fid,file=TRIM(fname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(i/=0 .or. buf(7:18)=='ORBITAL   CI') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_ci_coeff_from_gms: no GVB CI coeff&
                   &icients found'
  write(6,'(A)') 'in file '//TRIM(fname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf   ! skip one line

 do i = 1, npair, 1
  read(fid,'(A)') buf
  read(buf,*) j, j, j, coeff(1,i), coeff(2,i)
 end do ! for i

 close(fid)
end subroutine read_ci_coeff_from_gms

! perform Mulliken population analysis based on density matrix
!subroutine mulliken_pop_of_dm(nshl, shl2atm, ang, ibas, cart, nbf, P, S, natom, eff_nuc)
! implicit none
! integer :: i, j
! integer, intent(in) :: nshl, nbf, natom
! integer, intent(in) :: shl2atm(nshl), ang(nshl), ibas(nshl), eff_nuc(natom)
! integer, allocatable :: bfirst(:) ! size natom+1
! ! bfirst: the beginning index of basis func. of each atom
! real(kind=8), intent(in) :: P(nbf,nbf), S(nbf,nbf)
! real(kind=8) :: rtmp, ddot
! real(kind=8), allocatable :: gross(:) ! size natom
! logical, intent(in) :: cart
!
! allocate(bfirst(natom+1))
! call get_bfirst(nshl, shl2atm, ang, ibas, cart, natom, bfirst)
!
! allocate(gross(natom), source=0d0)
!
! do i = 1, natom, 1
!  rtmp = 0d0
!  do j = bfirst(i), bfirst(i+1)-1, 1
!   rtmp = rtmp + ddot(nbf, P(j,:), 1, S(:,j), 1)
!  end do ! for j
!  gross(i) = rtmp
! end do ! for i
!
! deallocate(bfirst)
! gross = eff_nuc - gross
! write(6,'(/,A)') 'Mulliken population:'
!
! do i = 1, natom, 1
!  write(6,'(I5,1X,F11.6)') i, gross(i)
! end do ! for i
!
! deallocate(gross)
!end subroutine mulliken_pop_of_dm

! perform lowdin/lowedin population analysis based on density matrix
!subroutine lowdin_pop_of_dm(nshl, shl2atm, ang, ibas, cart, nbf, P, S, natom, eff_nuc)
! implicit none
! integer :: i, j, m, lwork, liwork
! integer, intent(in) :: nshl, nbf, natom
! integer, intent(in) :: shl2atm(nshl), ang(nshl), ibas(nshl), eff_nuc(natom)
! integer, allocatable :: bfirst(:) ! size natom+1
! ! bfirst: the beginning index of basis func. of each atom
! integer, allocatable :: iwork(:), isuppz(:)
! real(kind=8), intent(in) :: P(nbf,nbf), S(nbf,nbf)
! real(kind=8) :: rtmp, ddot
! real(kind=8), allocatable :: gross(:) ! size natom
! real(kind=8), allocatable :: work(:), e(:), ev(:,:), sqrt_S_P(:,:), S0(:,:)
! ! e: eigenvalues, ev: eigenvectors, sqrt_S_P: S^(1/2)*P
! logical, intent(in) :: cart
!
! allocate(bfirst(natom+1))
! call get_bfirst(nshl, shl2atm, ang, ibas, cart, natom, bfirst)
!
! allocate(e(nbf), ev(nbf, nbf), isuppz(2*nbf))
! lwork = -1; liwork = -1
! allocate(work(1), iwork(1))
! allocate(S0(nbf,nbf), source=S)
! call dsyevr('V', 'A', 'L', nbf, S0, nbf, 0d0, 0d0, 0, 0, 1d-6, m, e, ev,&
!             nbf, isuppz, work, lwork, iwork, liwork, i)
! lwork = CEILING(work(1))
! liwork = iwork(1)
! deallocate(work, iwork)
! allocate(work(lwork), iwork(liwork))
! call dsyevr('V', 'A', 'L', nbf, S0, nbf, 0d0, 0d0, 0, 0, 1d-6, m, e, ev,&
!             nbf, isuppz, work, lwork, iwork, liwork, i)
! deallocate(isuppz, work, iwork)
!
! S0 = 0d0
! forall(i = 1:nbf) S0(i,i) = DSQRT(DABS(e(i)))
! deallocate(e)
! allocate(sqrt_S_P(nbf,nbf))
! call dsymm('R', 'L', nbf, nbf, 1d0, S0, nbf, ev, nbf, 0d0, sqrt_S_P, nbf)
! call dgemm('N', 'T', nbf, nbf, nbf, 1d0, sqrt_S_P, nbf, ev, nbf, 0d0, S0, nbf)
! deallocate(ev, sqrt_S_P)
!
! allocate(gross(natom), source=0d0)
! allocate(e(nbf))
!
! do i = 1, natom, 1
!  rtmp = 0d0
!
!  do j = bfirst(i), bfirst(i+1)-1, 1
!   e = S0(j,:)
!   do m = 1, nbf, 1
!    rtmp = rtmp + ddot(nbf,e,1,P(:,m),1)*S0(m,j)
!   end do ! for m
!  end do ! for j
!
!  gross(i) = rtmp
! end do ! for i
!
! deallocate(bfirst, e, S0)
! gross = eff_nuc - gross
!
! write(6,'(/,A)') 'Lowdin population:'
! do i = 1, natom, 1
!  write(6,'(I5,1X,F11.6)') i, gross(i)
! end do ! for i
!
! deallocate(gross)
!end subroutine lowdin_pop_of_dm

! perform SVD on MOs in two .fch(k) files
subroutine mo_svd_in2fch(fchname1, fchname2, idx1, idx2)
 implicit none
 integer :: nbf, nif, nmo
 integer, intent(in) :: idx1, idx2
!f2py intent(in) :: idx1, idx2
 character(len=240), intent(in) :: fchname1, fchname2
!f2py intent(in) :: fchname1, fchname2
 real(kind=8), allocatable :: S(:,:), mo_ovlp(:,:), mo1(:,:), mo2(:,:), u(:,:),&
  vt(:,:), ev(:)

 call read_nbf_and_nif_from_fch(fchname1, nbf, nif)
 nmo = idx2 - idx1 + 1
 allocate(mo1(nbf,nif), mo2(nbf,nif), S(nbf,nbf))
 call read_mo_from_fch(fchname1, nbf, nif, 'a', mo1)
 call read_mo_from_fch(fchname2, nbf, nif, 'a', mo2)
 call get_ao_ovlp_using_fch(fchname1, nbf, S)

 ! compute the MO-basis overlap matrix (C1^T)SC2
 allocate(mo_ovlp(nmo,nmo))
 call calc_CTSCp(nbf, nmo, mo1(:,idx1:idx2), S, mo2(:,idx1:idx2), mo_ovlp)
 deallocate(mo1, mo2, S)

 allocate(u(nmo,nmo), vt(nmo,nmo), ev(nmo))
 call do_svd(nmo, nmo, mo_ovlp, u, vt, ev)
 deallocate(mo_ovlp, u, vt)

 write(6,'(/,A)') 'SVD analysis of two sets of MOs:'
 write(6,'(A)') 'All singular values:'
 write(6,'(5(1X,ES15.8))') ev
 write(6,'(/,A,ES15.8)') 'The smallest singular value:', MINVAL(ev)
 write(6,'(A,I0)') 'No. singular values > 0.9: ', COUNT(ev>0.9d0)
 write(6,'(A,I0)') 'No. singular values < 0.1: ', COUNT(ev<0.1d0)
 write(6,'(A,I0)') 'No. singular values <0.01: ', COUNT(ev<0.01d0)
 deallocate(ev)
end subroutine mo_svd_in2fch

