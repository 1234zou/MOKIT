! written by jxzou at 20180331: generate UHF natural orbtials (UNO) from UHF canonical orbitals
! updated by jxzou at 20180420: to support na > nb
! updated by jxzou at 20180520: modify the intent of array alpha_coeff from (inout) to (in,copy)
! updated by jxzou at 20180825: fix the bug when only 1 pair
! updated by jxzou at 20191215: delete SVD generation of virtual/inactive MOs (PAO outside recommended)
! updated by jxzou at 20200426: add NOON output
! updated by jxzou at 20210518: add an intent(in) parameter ON_thres
! updated by jxzou at 20220711: change ON_thres to uno_thres

subroutine check_unity(n, a, maxv, abs_mean)
 implicit none
 integer :: i
 integer, intent(in) :: n
 real(kind=8), intent(in) :: a(n,n)
 real(kind=8), intent(out) :: maxv, abs_mean
 real(kind=8), allocatable :: b(:,:)

 allocate(b(n,n), source=a)
 forall(i = 1:n) b(i,i) = b(i,i) - 1d0
 b = DABS(b)
 maxv = MAXVAL(b)
 abs_mean = SUM(b)/DBLE(n*n)
 deallocate(b)
end subroutine check_unity

subroutine svd_and_rotate(nbf, na, nb, ab_ovlp, mo_a, mo_b, sv, reverse)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, na, nb
 real(kind=8), intent(in) :: ab_ovlp(na,nb)
 real(kind=8), intent(inout) :: mo_a(nbf,na), mo_b(nbf,nb)
 real(kind=8), intent(out) :: sv(na)
 real(kind=8), allocatable :: u(:,:), vt(:,:), u2(:,:), vt2(:,:),  new_mo(:,:),&
  sv2(:)
 logical, intent(in) :: reverse

 ! perform SVD on alpha_beta_ovlp of occupied or virtual orbitals
 allocate(u(na,na), vt(nb,nb))
 call do_svd(na, nb, ab_ovlp, u, vt, sv)

 ! rotate occupied or virtual orbitals
 if(reverse) then
  allocate(sv2(na))
  do i = 1, na, 1
   sv2(i) = sv(na-i+1)
  end do ! for i
  sv = sv2
  deallocate(sv2)
  allocate(u2(na,na))
  do i = 1, na, 1
   u2(:,i) = u(:,na-i+1)
  end do ! for i
  u = u2
  deallocate(u2)
  allocate(vt2(nb,nb))
  do i = 1, nb, 1
   vt2(i,:) = vt(nb-i+1,:)
  end do ! for i
  vt = vt2
  deallocate(vt2)
 end if

 allocate(new_mo(nbf,na), source=0d0)
 call dgemm('N', 'N', nbf, na, na, 1d0, mo_a, nbf, u, na, 0d0, new_mo, nbf)
 mo_a = new_mo
 deallocate(u, new_mo)

 allocate(new_mo(nbf,nb), source=0d0)
 call dgemm('N', 'T', nbf, nb, nb, 1d0, mo_b, nbf, vt, nb, 0d0, new_mo, nbf)
 mo_b = new_mo
 deallocate(vt, new_mo)
end subroutine svd_and_rotate

! This subroutine is designed to be imported as a module in Python.
subroutine uno(outname, nbf, nif, na, nb, mo_a, mo_b, ao_ovlp, uno_thres, idx, noon, &
               uno_coeff)
 implicit none
 integer :: i, ndb, nact, nact0, nopen, nocc, fid
 ! ndb: the number of doubly occupied MOs
 ! nact: the number of active occupied orbitals
 ! nact = nact0 + nopen
 ! nocc = na + nact0
 integer, intent(in) :: nbf, nif, na, nb
!f2py intent(in) :: nbf, nif, na, nb
 ! nbf: the number of basis functions
 ! nif: the number of independent functions, i.e., the number of MOs
 ! na: the number of alpha electrons
 ! nb: the number of beta electrons

 ! array idx is in Fortran convention, i.e., begins from 1, not 0
 integer, intent(out) :: idx(3) ! 1/2/3: start of pair/start of virtual/nopen
!f2py intent(out) :: idx
 integer, allocatable :: idx1(:), idx2(:)

 real(kind=8) :: maxv, abs_mean
 real(kind=8), intent(in) :: uno_thres ! UNO occupation number threshold
!f2py intent(in) :: uno_thres
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
!f2py intent(in) :: ao_ovlp
!f2py depend(nbf) :: ao_ovlp
 real(kind=8), intent(in) :: mo_a(nbf,nif), mo_b(nbf,nif)
!f2py intent(in) :: mo_a, mo_b
!f2py depend(nbf,nif) :: mo_a, mo_b
 real(kind=8), intent(out) :: noon(nif), uno_coeff(nbf,nif)
!f2py intent(out) :: noon, uno_coeff
!f2py depend(nif) :: noon
!f2py depend(nbf,nif) :: uno_coeff
 real(kind=8), parameter :: ON_criteria = 1d-5
 real(kind=8), allocatable :: mo_ovlp(:,:), occ_a(:,:), occ_b(:,:), sv_occ0(:),&
  sv_occ(:)
 character(len=*), intent(in) :: outname
 character(len=63), parameter :: on_warn1 = 'Warning in subroutine uno: uno_th&
                                            &res deviates from ON_criteria.'
 character(len=35), parameter :: on_warn2 = 'You better know what you are doing.'

 noon = 0d0
 nopen = na - nb
 if(nopen < 0) then
  write(6,'(/,A)') 'ERROR in subroutine uno: na < nb.'
  write(6,'(2(A,I0))') 'na=', na, ', nb=', nb
  stop
 end if

 open(newunit=fid,file=TRIM(outname),status='replace')
 write(fid,'(A4,I5)') 'nbf=', nbf
 write(fid,'(A4,I5)') 'nif=', nif
 write(fid,'(A,F11.7)') 'ON_criteria=', ON_criteria
 write(fid,'(A,F11.7)') 'uno_thres=', uno_thres
 if(DABS(uno_thres - ON_criteria) > 1d-5) then
  write(fid,'(/,A)') on_warn1
  write(fid,'(A)') on_warn2
  write(6,'(/,A)') on_warn1
  write(6,'(A)') on_warn2
 end if

 ! check the orthonormality of initial Alpha and Beta MO, respectively
 allocate(mo_ovlp(nif,nif))
 call calc_CTSC(nbf, nif, mo_a, ao_ovlp, mo_ovlp)
 call check_unity(nif, mo_ovlp, maxv, abs_mean)
 write(fid,'(/,A)') 'The orthonormality of initial Alpha MO:'
 write(fid,'(A,F16.10)') 'maxv=', maxv
 write(fid,'(A,F16.10)') 'abs_mean=', abs_mean

 call calc_CTSC(nbf, nif, mo_b, ao_ovlp, mo_ovlp)
 call check_unity(nif, mo_ovlp, maxv, abs_mean)
 deallocate(mo_ovlp)
 write(fid,'(/,A)') 'The orthonormality of initial Beta MO:'
 write(fid,'(A,F16.10)') 'maxv=', maxv
 write(fid,'(A,F16.10)') 'abs_mean=', abs_mean
 ! check orthonormality done

 uno_coeff = mo_a
 if(nb == 0) then ! no beta electrons, return
  forall(i = 1:nopen) noon(i) = 1d0
  idx = [1, nopen+1, nopen]
  write(fid,'(/,A6,I5)') 'ndb  =', 0
  write(fid,'(A6,I5)')   'nact =', na
  write(fid,'(A6,I5)')   'nact0=', 0
  write(fid,'(A6,3I5)')  'idx  =', idx
  close(fid)
  return
 end if

 allocate(occ_a(nbf,na), source=mo_a(:,1:na))
 allocate(occ_b(nbf,nb), source=mo_b(:,1:nb))

 ! calculate the overlap between alpha and beta occupied spatial orbitals
 allocate(mo_ovlp(na,nb))
 call calc_CTSCp2(nbf, na, nb, occ_a, ao_ovlp, occ_b, mo_ovlp)
 ! calculate done

 ! do SVD on the alpha_beta_ovlp of occupied spatial orbitals
 allocate(sv_occ(na))
 call svd_and_rotate(nbf, na, nb, mo_ovlp, occ_a, occ_b, sv_occ, .false.)
 deallocate(mo_ovlp)
 !write(6,'(/,A)') 'Singular values from SVD of Alpha/Beta MOs:'
 !write(6,'(5(1X,ES15.8))') (sv_occ(i), i=1,na)
 ! SVD done in occ space

 allocate(sv_occ0(na), source=sv_occ)
 nact = COUNT(sv_occ < 1d0-ON_criteria)
 ndb = na - nact
 nact0 = nact - nopen
 nocc = na + nact0
 idx = [ndb+1, nocc+1, nopen]
 write(fid,'(/,A6,I5)') 'ndb  =', ndb
 write(fid,'(A6,I5)')   'nact =', nact
 write(fid,'(A6,I5)')   'nact0=', nact0
 write(fid,'(A6,3I5)')  'idx  =', idx

 ! generate NOON (Natural Orbital Occupation Number)
 forall(i = 1:na) noon(i) = 1d0 + sv_occ(i)
 forall(i = 1:nact0) noon(na+i) = 2d0 - noon(nb-i+1)
 deallocate(sv_occ)

 ! copy the doubly occupied MOs
 uno_coeff(:,1:ndb) = occ_a(:,1:ndb)

 ! copy the singly occupied MO
 if(nopen > 0) uno_coeff(:,nb+1:na) = occ_a(:,nb+1:na)

 ! transform the corresponding orbitals to UNOs in occ space
 allocate(idx1(nact0), idx2(nact0))
 forall(i = 1:nact0)
  idx1(i) = ndb + i
  idx2(i) = nocc + 1 - i
 end forall
 forall(i = 1:nact0)
  uno_coeff(:,idx1(i)) = (occ_a(:,idx1(i)) + occ_b(:,idx1(i)))/DSQRT(2d0*noon(idx1(i)))
  uno_coeff(:,idx2(i)) = (occ_a(:,idx1(i)) - occ_b(:,idx1(i)))/DSQRT(2d0*noon(idx2(i)))
 end forall
 deallocate(idx1, idx2, occ_a, occ_b)
 ! done transform in occ space

 ! Set virtual MOs to be zero. They are supposed to be calculated in subsequent
 ! PAO constructions
 if(nocc < nif) uno_coeff(:,nocc+1:nif) = 0d0

 ! check the orthonormality of final Alpha MO
 allocate(mo_ovlp(nocc,nocc))
 call calc_CTSC(nbf, nocc, uno_coeff(:,1:nocc), ao_ovlp, mo_ovlp)
 call check_unity(nocc, mo_ovlp, maxv, abs_mean)
 deallocate(mo_ovlp)
 write(fid,'(/,A)') 'The orthonormality of final Alpha MO:'
 write(fid,'(A,F16.10)') 'maxv=', maxv
 write(fid,'(A,F16.10)') 'abs_mean=', abs_mean

 ! now let's update ndb, nact, nact0, idx according to input uno_thres
 if(DABS(uno_thres - ON_criteria) > 1d-5) then
  nact = COUNT(sv_occ0 < 1d0-uno_thres)
  ndb = na - nact
  nact0 = nact - nopen
  nocc = na + nact0
  idx = [ndb+1, nocc+1, nopen]
 end if

 write(fid,'(/,A6,I5)') 'ndb  =', ndb
 write(fid,'(A6,I5)')   'nact =', nact
 write(fid,'(A6,I5)')   'nact0=', nact0
 write(fid,'(A6,3I5)')  'idx  =', idx
 close(fid)
 deallocate(sv_occ0)
end subroutine uno

subroutine ovlp_sv2c1_c2(d, c1, c2)
 implicit none
 real(kind=8) :: d1, r1, r2, r3, r4
 real(kind=8), intent(in) :: d
 real(kind=8), intent(out) :: c1, c2

 d1 = MAX(-1d0, MIN(d, 1d0))
 r1 = 1d0 + d1; r2 = 1d0 - d1
 r3 = DSQRT(r1/r2); r4 = DSQRT(r2/r1)
 c1 = r3 + r4
 c2 = r3 - r4
end subroutine ovlp_sv2c1_c2

! Generate UNOs from UHF MOs, then swap the specified UNO pair (e.g. HONO/LUNO)
! and back-transform to get new UHF MOs. Singly occ. MOs will not be swapped.
! `swap(nb)=.True.` means swapping HONO and LUNO.
! `swap(nb-1)=.True.` means swapping HONO-1 and LUNO+1.
subroutine swap_pair_in_uno(na, nb, nbf, nif, swap, ao_ovlp, mo_a0, mo_b0, &
                            mo_a, mo_b)
 implicit none
 integer :: i, nopen
 integer, intent(in) :: na, nb, nbf, nif
!f2py intent(in) :: na, nb, nbf, nif
 real(kind=8) :: maxv, abs_mean, c1, c2
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf), mo_a0(nbf,nif), mo_b0(nbf,nif)
!f2py intent(in) :: ao_ovlp, mo_a0, mo_b0
!f2py depend(nbf) :: ao_ovlp
!f2py depend(nbf,nif) :: mo_a0, mo_b0
 real(kind=8), intent(out) :: mo_a(nbf,nif), mo_b(nbf,nif)
!f2py intent(out) :: mo_a, mo_b
!f2py depend(nbf,nif) :: mo_a, mo_b
 real(kind=8), allocatable :: mo_ovlp(:,:), occ_a(:,:), occ_b(:,:), sv_occ(:)
 logical, intent(in) :: swap(nb)
!f2py intent(in) :: swap
!f2py depend(nb) :: swap

 nopen = na - nb
 if(nopen < 0) then
  write(6,'(/,A)') 'ERROR in subroutine swap_pair_in_uno: na < nb.'
  write(6,'(2(A,I0))') 'na=', na, ', nb=', nb
  stop
 end if

 ! check the orthonormality of initial Alpha and Beta MO, respectively
 allocate(mo_ovlp(nif,nif))
 call calc_CTSC(nbf, nif, mo_a0, ao_ovlp, mo_ovlp)
 call check_unity(nif, mo_ovlp, maxv, abs_mean)
 write(6,'(/,A)') 'The orthonormality of initial Alpha MO:'
 write(6,'(A,F16.10)') 'maxv=', maxv
 write(6,'(A,F16.10)') 'abs_mean=', abs_mean

 call calc_CTSC(nbf, nif, mo_b0, ao_ovlp, mo_ovlp)
 call check_unity(nif, mo_ovlp, maxv, abs_mean)
 deallocate(mo_ovlp)
 write(6,'(/,A)') 'The orthonormality of initial Beta MO:'
 write(6,'(A,F16.10)') 'maxv=', maxv
 write(6,'(A,F16.10)') 'abs_mean=', abs_mean
 ! check orthonormality done

 allocate(occ_a(nbf,na), source=mo_a0(:,1:na))
 allocate(occ_b(nbf,nb), source=mo_b0(:,1:nb))
 allocate(mo_ovlp(na,nb))
 call calc_CTSCp2(nbf, na, nb, occ_a, ao_ovlp, occ_b, mo_ovlp)

 allocate(sv_occ(na))
 call svd_and_rotate(nbf, na, nb, mo_ovlp, occ_a, occ_b, sv_occ, .false.)
 deallocate(mo_ovlp)
 mo_a(:,1:na) = occ_a
 mo_b(:,1:nb) = occ_b

 do i = 1, nb, 1
  if(swap(i)) then
   call ovlp_sv2c1_c2(sv_occ(i), c1, c2)
   mo_a(:,i) = 0.5d0*(c1*occ_a(:,i) - c2*occ_b(:,i))
   mo_b(:,i) = 0.5d0*(c2*occ_a(:,i) - c1*occ_b(:,i))
  end if
 end do ! for i
 deallocate(sv_occ, occ_a, occ_b)

 ! Assuming mo_a(:,k) is updated in the above loop, mo_a(:,1:na) are orthonormal,
 ! but mo_a(:,k) is not orthogonal to alpha virtual MOs mo_a(:,na+1:nif), so we
 ! need to construct alpha virtual MOs. Similarly, mo_b(:,1:nb) are orthonormal,
 ! but mo_b(:,k) is not orthogonal to beta virtual MOs mo_b(:,nb+1:nif), so we
 ! need to construct beta virtual MOs.
 allocate(occ_a(nbf,nif), source=mo_a)
 call construct_vir(nbf, nif, na+1, occ_a, ao_ovlp, mo_a)
 occ_a = mo_b
 call construct_vir(nbf, nif, nb+1, occ_a, ao_ovlp, mo_b)
 deallocate(occ_a)
end subroutine swap_pair_in_uno

! Enlarge the active space by utilizing two sets of singular values:
!  1) SVD of {doubly occupied space of mo1, active space of mo2}
!  2) SVD of {virtual space of mo1, active space of mo2}
! mo1 and mo2 are assumed to be expanded on the same set of basis functions,
!  and thus share the same AO-based overlap integrals.
! A commonly used example:
!  mo1 = singlet ground state CASSCF MOs
!  mo2 = triplet UNO or CASSCF NOs
subroutine enlarge_as_by_svd(idx, nbf, nif, mo1, mo2, ao_ovlp, new_mo)
 implicit none
 integer :: ndb1, nvir1, nact2, ndb_add, nvir_add
 integer, intent(in) :: idx(4), nbf, nif
!f2py intent(in) :: idx, nbf, nif
! idx(1)/(2): the index of the 1st/last active orbital in mo1
! idx(3)/(4): the index of the 1st/last active orbital in mo2
! all four integers are in Fortran convention, i.e. start from 1
 real(kind=8), parameter :: thres = 0.75d0
 real(kind=8), intent(in) :: mo1(nbf,nif), mo2(nbf,nif), ao_ovlp(nbf,nbf)
!f2py intent(in) :: mo1, mo2, ao_ovlp
!f2py depend(nbf,nif) :: mo1, mo2
!f2py depend(nbf) :: ao_ovlp
 real(kind=8), intent(out) :: new_mo(nbf,nif)
!f2py depend(nbf,nif) :: new_mo
!f2py intent(out) :: new_mo
 real(kind=8), allocatable :: mo3(:,:), mo4(:,:), mo_ovlp(:,:), sv(:)

 new_mo = mo1
 ndb1 = idx(1) - 1
 nvir1 = nif - idx(2)
 nact2 = idx(4) - idx(3) + 1

 allocate(mo3(nbf,ndb1), source=mo1(:,1:ndb1))
 allocate(mo4(nbf,nact2), source=mo2(:,idx(3):idx(4)))
 allocate(mo_ovlp(ndb1,nact2))
 call calc_CTSCp2(nbf, ndb1, nact2, mo3, ao_ovlp, mo4, mo_ovlp)
 allocate(sv(ndb1))
 call svd_and_rotate(nbf, ndb1, nact2, mo_ovlp, mo3, mo4, sv, .true.)
 deallocate(mo_ovlp)
 write(6,'(A)') 'Singular values of docc1 v.s. act2:'
 write(6,'(5(1X,ES15.8))') sv
 ndb_add = COUNT(sv > thres)
 write(6,'(A,I0)') 'ndb_add=', ndb_add
 if(ndb_add > 0) new_mo(:,1:ndb1) = mo3
 deallocate(sv, mo3)

 mo4 = mo2(:,idx(3):idx(4))
 allocate(mo3(nbf,nvir1), source=mo1(:,idx(2)+1:nif))
 allocate(mo_ovlp(nvir1,nact2))
 call calc_CTSCp2(nbf, nvir1, nact2, mo3, ao_ovlp, mo4, mo_ovlp)
 allocate(sv(nvir1))
 call svd_and_rotate(nbf, nvir1, nact2, mo_ovlp, mo3, mo4, sv, .false.)
 write(6,'(/,A)') 'Singular values of vir1 v.s. act2:'
 write(6,'(5(1X,ES15.8))') sv
 deallocate(mo4, mo_ovlp)
 nvir_add = COUNT(sv > thres)
 write(6,'(A,I0)') 'nvir_add=', nvir_add
 if(nvir_add > 0) new_mo(:,idx(2)+1:nif) = mo3
 deallocate(sv, mo3)
end subroutine enlarge_as_by_svd

! calculated the SVD singular values of overlap of two sets of MOs
! Note: the input coeff1 and coeff2 must have the same dimension (nbf,nif)
subroutine svd_of_two_mo(nbf, nif, ao_ovlp, old_mo1, old_mo2, new_mo1, new_mo2)
 implicit none
 integer :: i
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf), old_mo1(nbf,nif), old_mo2(nbf,nif)
!f2py intent(in) :: ao_ovlp, old_mo1, old_mo2
!f2py depend(nbf,nif) :: old_mo1, old_mo2
!f2py depend(nbf) :: ao_ovlp
 real(kind=8), intent(out) :: new_mo1(nbf,nif), new_mo2(nbf,nif)
!f2py intent(out) :: new_mo1, new_mo2
!f2py depend(nbf,nif) :: new_mo1, new_mo2
 real(kind=8), allocatable :: mo_ovlp(:,:), ev(:)

 allocate(mo_ovlp(nif,nif))
 call calc_CTSCp(nbf, nif, old_mo1, ao_ovlp, old_mo2, mo_ovlp)
 allocate(ev(nif))
 new_mo1 = old_mo1
 new_mo2 = old_mo2
 call svd_and_rotate(nbf, nif, nif, mo_ovlp, new_mo1, new_mo2, ev, .False.)
 deallocate(mo_ovlp)

 write(6,'(/,A)') 'SVD analysis of two sets of MOs:'
 write(6,'(A,ES15.8)') 'The smallest singular value:', MINVAL(ev)
 i = COUNT(ev < 1d-1)
 write(6,'(A,I0)') 'Number of singular values< 0.1: ', i
 i = COUNT(ev < 1d-2)
 write(6,'(A,I0)') 'Number of singular values<0.01: ', i

 write(6,'(A)') 'All singular values:'
 write(6,'(5(1X,ES15.8))') (ev(i),i=1,nif)
 deallocate(ev)
end subroutine svd_of_two_mo

