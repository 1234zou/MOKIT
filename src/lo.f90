! written by jxzou at 20190723
! updated by jxzou at 20200411: add Pipek-Mezey orbital localization (DOI: 10.1063/1.456588)
! updated by jxzou at 20200413: add Cholesky decomposition LMOs (DOI: 10.1063/1.2360264)
! updated by jxzou at 20220520: generate NO from a NSO .fch file

! Note: before PySCF-1.6.4, its dumped .molden file is wrong when using Cartesian functions.

! Please use subroutine gen_no_from_density_and_ao_ovlp in rwwfn.f90, which has
! the same functionality as the subroutine no

! localize singly occupied orbitals in a .fch file
subroutine localize_singly_occ_orb(fchname, pm_loc)
 implicit none
 integer :: na, nb
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: pm_loc

 call read_na_and_nb_from_fch(fchname, na, nb)
 if(na < nb) then
  write(6,'(/,A)') 'ERROR in subroutine localize_singly_occ_orb: na<nb.'
  write(6,'(A)') 'Something must be wrong. Check file '//TRIM(fchname)
  stop
 else if(na == nb) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine localize_singly_occ_orb: no singly oc&
                 &cupied orbital'
  write(6,'(A)') 'to be localized.'
  write(6,'(A)') REPEAT('-',79)
 else if(na == nb+1) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine localize_singly_occ_orb: only one sin&
                 &gly occupied'
  write(6,'(A)') 'orbital, no need to perform orbital localization.'
  write(6,'(A)') REPEAT('-',79)
 else
  write(6,'(A)') 'Perform orbital localization on singly occupied orbitals...'
  call localize_orb(fchname, nb+1, na, pm_loc)
 end if
end subroutine localize_singly_occ_orb

! the loc() function in mokit.lib.gaussian is called to localize specified
! orbitals in a .fch(k) file using the PM method
subroutine localize_orb(fchname, i1, i2, pm_loc)
 implicit none
 integer :: i, fid
 integer, intent(in) :: i1, i2 ! Fortran convention
 character(len=240) :: lmofch, pyname, outname
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: pm_loc

 call find_specified_suffix(fchname, '.fch', i)
 lmofch = fchname(1:i-1)//'_LMO.fch'
 pyname = fchname(1:i-1)//'.py'
 outname = fchname(1:i-1)//'.out'

 if(i2 < i1) then
  write(6,'(/,A)') 'ERROR in subroutine localize_orb: i2<i1. Invalid values.'
  write(6,'(2(A,I0))') 'i1=', i1, ', i2=', i2
  stop
 end if

 open(newunit=fid,file=TRIM(pyname),status='replace')
 write(fid,'(A)') 'from shutil import copyfile'
 write(fid,'(A)') 'from mokit.lib.gaussian import load_mol_from_fch'
 write(fid,'(A)') 'from mokit.lib.rwwfn import read_nbf_and_nif_from_fch, \'
 write(fid,'(A)') ' read_eigenvalues_from_fch'
 write(fid,'(A)') 'from mokit.lib.fch2py import fch2py'
 write(fid,'(A)') 'from mokit.lib.py2fch import py2fch'
 if(pm_loc) then
  write(fid,'(A)') 'from mokit.lib.lo import pm'
 else
  write(fid,'(A)') 'from pyscf.lo.boys import dipole_integral'
  write(fid,'(A)') 'from mokit.lib.lo import boys'
 end if
 write(fid,'(A)') 'import os'

 write(fid,'(/,A)') "fchname = '"//TRIM(fchname)//"'"
 write(fid,'(A)') "lmofch = '"//TRIM(lmofch)//"'"
 write(fid,'(A)') 'mol = load_mol_from_fch(fchname)'
 write(fid,'(A)') 'nbf, nif = read_nbf_and_nif_from_fch(fchname)'
 write(fid,'(A)') "mo_coeff = fch2py(fchname, nbf, nif, 'a')"
 write(fid,'(2(A,I0),A)') 'idx = range(',i1-1,',',i2,')'
 write(fid,'(A)') 'nmo = len(idx)'
 if(pm_loc) then
  write(fid,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
  write(fid,'(A)') 'loc_orb = pm(mol.nbas, mol._bas[:,0],mol._bas[:,1],mol._bas&
                   &[:,3], mol.cart,'
  write(fid,'(A)') "             nbf, nmo, mo_coeff[:,idx], S, 'mulliken')"
 else
  write(fid,'(A)') 'mo_dipole = dipole_integral(mol, mo_coeff[:,idx])'
  write(fid,'(A)') 'loc_orb = boys(nbf, nmo, mo_coeff[:,idx], mo_dipole)'
 end if
 write(fid,'(A)') 'mo_coeff[:,idx] = loc_orb.copy()'
 write(fid,'(A)') "noon = read_eigenvalues_from_fch(fchname, nif, 'a')"
 write(fid,'(A)') 'copyfile(fchname, lmofch)'
 write(fid,'(A)') "py2fch(lmofch,nbf,nif,mo_coeff,'a',noon,False,False)"
 write(fid,'(A)') 'os.rename(lmofch, fchname)'
 close(fid)

 call submit_pyscf_job(pyname, .true.)
 call delete_files(2, [pyname, outname])
end subroutine localize_orb

! generate AO-basis density matrix based on a .fch file
! Note: be careful about the itype!
subroutine gen_ao_dm_from_fch(fchname, itype, nbf, dm)
 implicit none
 integer :: i, nbf0, nif, na, nb
 integer, intent(in) :: itype, nbf
 real(kind=8), allocatable :: noon(:), n(:,:), mo(:,:), Cn(:,:)
 real(kind=8), intent(out) :: dm(nbf,nbf)
 character(len=240), intent(in) :: fchname

 dm = 0d0 ! initialization

 call read_nbf_and_nif_from_fch(fchname, nbf0, nif)
 if(nbf0 /= nbf) then
  write(6,'(/,A)') 'ERROR in subroutine gen_ao_dm_from_fch: nbf0/=nbf.'
  write(6,'(A,2I6)') 'nbf0, nbf=', nbf0, nbf
  stop
 end if

 call read_na_and_nb_from_fch(fchname, na, nb)
 allocate(noon(nif), source=0d0)

 select case(itype)
 case(0) ! R(O)HF orbitals, occupation 0/1/2
  noon(1:nb) = 2d0
  if(na > nb) noon(nb+1:na) = 1d0
 case(1) ! UHF alpha MOs, occupation 0/1
  noon(1:na) = 1d0
 case(2) ! UHF beta MOs, occupation 0/1
  noon(1:nb) = 1d0
 case(3) ! NSO alpha, fractional occupation
  call read_eigenvalues_from_fch(fchname, nif, 'a', noon)
 case(4) ! NSO beta, fractional occupation
  call read_eigenvalues_from_fch(fchname, nif, 'b', noon)
 case default
  write(6,'(A)') 'ERROR in subroutine gen_ao_dm_from_fch: itype out of range.'
  write(6,'(A,I0)') 'Only 0~4 are allowed. But got itype=', itype
  stop
 end select

 allocate(n(nif,nif),source=0d0)
 forall(i = 1:nif) n(i,i) = noon(i)
 deallocate(noon)

 allocate(mo(nbf,nif), source=0d0)
 select case(itype)
 case(0,1,3)
  call read_mo_from_fch(fchname, nbf, nif, 'a', mo)
 case(2,4)
  call read_mo_from_fch(fchname, nbf, nif, 'b', mo)
 end select

 allocate(Cn(nbf,nif), source=0d0)
 ! P = Cn(C^T)
 call dsymm('R', 'L', nbf, nif, 1d0, mo, nbf, n, nif, 0d0, Cn, nbf)
 deallocate(n)
 call dgemm('N', 'T', nbf, nbf, nif, 1d0, Cn, nbf, mo, nbf, 0d0, dm, nbf)
 deallocate(Cn, mo)
end subroutine gen_ao_dm_from_fch

! generate spatial natural orbitals from natural spin orbitals
subroutine gen_no_from_nso(fchname)
 use util_wrapper, only: fch_u2r_wrap
 implicit none
 integer :: i, nbf, nif
 character(len=240) :: no_fch
 character(len=240), intent(in) :: fchname ! must have NSO in it
!f2py intent(in) :: fchname
 real(kind=8), allocatable :: noon(:), dm(:,:), mo(:,:), S(:,:)

 i = INDEX(fchname, '.fch', back=.true.)
 no_fch = fchname(1:i-1)//'_NO.fch'
 call fch_u2r_wrap(fchname, no_fch)

 call read_nbf_and_nif_from_fch(fchname, nbf, nif)
 allocate(dm(nbf,nbf))
 call read_dm_from_fch(fchname, 1, nbf, dm)

 allocate(S(nbf,nbf))
 call get_ao_ovlp_using_fch(fchname, nbf, S)

 allocate(noon(nif), mo(nbf,nif))
 call gen_no_from_density_and_ao_ovlp(nbf, nif, dm, S, noon, mo)
 deallocate(dm, S)

 call write_mo_into_fch(no_fch, nbf, nif, 'a', mo)
 deallocate(mo)

 call write_eigenvalues_to_fch(no_fch, nif, 'a', noon, .true.)
 deallocate(noon)
end subroutine gen_no_from_nso

! compute MO-based density matrix
! Note: the result of this subroutine is identical to that of subroutine
!       solve_ON_matrix, but using a different formula
subroutine get_mo_based_dm(nbf, nif, coeff, S, P, dm)
 implicit none
 integer :: i, j
 integer :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8) :: coeff(nbf,nif), S(nbf,nbf), P(nbf,nbf), dm(nif,nif)
!f2py intent(in) :: coeff, S, P
!f2py intent(out) :: dm
!f2py depend(nbf,nif) :: coeff
!f2py depend(nbf) :: S, P
!f2py depend(nif) :: dm
 real(kind=8), allocatable :: SC(:,:), PSC(:,:)

 allocate(SC(nbf,nif), source=0d0)
 call dsymm('L', 'L', nbf, nif, 1d0, S, nbf, coeff, nbf, 0d0, SC, nbf)

 allocate(PSC(nbf,nif), source=0d0)
 call dsymm('L', 'L', nbf, nif, 1d0, P, nbf, SC, nbf, 0d0, PSC, nbf)

 call dgemm('T', 'N', nif, nif, nbf, 1d0, SC, nbf, PSC, nbf, 0d0, dm, nif)

 write(6,'(A)') 'MO-based density matrix (Final one electron symbolic density matrix):'
 do i = 1, nif, 1
  do j = i, nif, 1
   write(6,'(2I4,F15.8)') j, i, dm(j,i)
  end do ! for j
 end do ! for i

 deallocate(SC, PSC)
end subroutine get_mo_based_dm

! generate Coulson-Fischer orbitals from GVB natural orbitals
! this is actually orthogonal -> non-orthogonal orbital transformation
subroutine gen_cf_orb(datname, ndb, nopen)
 implicit none
 integer :: i, j, k, nbf, nif, npair
 integer, intent(in) :: ndb, nopen
 ! ndb and nopen cannot be determined from .dat file, manual input required
 real(kind=8) :: a, a2, fac
 real(kind=8), allocatable :: coeff(:,:), ci_coeff(:,:), rtmp(:,:)
 ! coeff: MO coefficients
 ! ci_coeff: GVB CI coefficients, pair coefficients
 ! rtmp: temporary array to hold two orbitals
 character(len=240), intent(in) :: datname

 if(ndb<0 .or. nopen<0) then
  write(6,'(A)') 'ERROR in subroutine gen_cf_orb: ndb<0 or nopen<0 found.'
  write(6,'(A)') 'Correct values should be assgined to these two parameters.'
  stop
 end if

 ! read Cartesian-type nbf and nif from GAMESS .dat file
 ! No matter spherical harmonic/Cartesian-type basis you use, the MOs in
 ! GAMESS .inp/.dat file are in Cartesian-type. So this subroutine reads
 ! Cartesian-type nbf
 call read_cart_nbf_nif_from_dat(datname, nbf, nif)

 allocate(coeff(nbf,nif))
 call read_mo_from_dat(datname, nbf, nif, coeff)

 call read_npair_from_dat(datname, npair)
 allocate(ci_coeff(2,npair))
 call read_ci_coeff_from_dat(datname, npair, ci_coeff)

 if(ANY(ci_coeff(1,:)<0d0) .or. ANY(ci_coeff(2,:)>0d0)) then
  write(6,'(A)') 'ERROR in subroutine gen_cf_orb: pair coefficients in file '&
                  //TRIM(datname)//' violate'
  write(6,'(A)') 'the rule C1>0, C2<0 for each pair. If your MOs and pair coef&
                 &ficients are correct,'
  write(6,'(A)') 'you need to swap the bonding orbital with the anti-bonding o&
                 &ne, and swap two pair coefficients.'
  write(6,'(A)') 'Then try to call this subroutine again.'
  stop
 end if

 k = ndb + nopen
 if(k+2*npair > nif) then
  write(6,'(A)') 'ERROR in subroutine gen_cf_orb: probably wrong ndb or nopen.'
  write(6,'(A,2I4)') 'Your input ndb, nopen=', ndb, nopen
  stop
 end if

 allocate(rtmp(nbf,2), source=0d0)
 write(6,'(/,A)') 'Non-orthogonal overlap for each pair:'

 do i = 1, npair, 1
  a2 = -ci_coeff(1,i)/ci_coeff(2,i)
  ! always set b=1, so no need to calculated b
  fac = 1d0/DSQRT(a2 + 1d0)
  a = DSQRT(a2)

  j = 2*i - 1
  rtmp = coeff(:,k+j:k+j+1)
  rtmp(:,1) = a*rtmp(:,1)
  coeff(:,k+j) = fac*(rtmp(:,1) + rtmp(:,2))
  coeff(:,k+j+1) = fac*(rtmp(:,1) - rtmp(:,2))
  write(6,'(A,I3,A,F10.6)') 'i=', i, ', S_i=', (a2-1d0)/(a2+1d0)
 end do ! for i

 deallocate(rtmp, ci_coeff)

 call write_mo_into_dat(datname, nbf, nif, coeff, .false.)
 deallocate(coeff)
end subroutine gen_cf_orb

! use Cholesky factorization/decomposition of the density matrix to generate LMOs
subroutine cholesky(nbf, nif, coeff, new_coeff)
 implicit none
 integer :: i, j, t0, t1, time
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif)
!f2py intent(in) :: coeff
!f2py depend(nbf,nif) :: coeff
 real(kind=8), intent(out) :: new_coeff(nbf,nif)
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: new_coeff
 real(kind=8) :: ddot, rtmp1, rtmp2
 real(kind=8), allocatable :: P(:,:)

 t0 = time()
 new_coeff = coeff
 if(nif == 1) then
  write(6,'(/,A)') 'ERROR in subroutine cholesky: only 1 orbital. Nothing to do.'
  stop
 end if

 allocate(P(nbf,nbf), source=0d0)
 call dgemm('N', 'T', nbf, nbf, nif, 1d0, coeff, nbf, coeff, nbf, 0d0, P, nbf)

 do i = 1, nif, 1
  rtmp1 = P(i,i) - ddot(i-1, new_coeff(i,1:i-1), 1, new_coeff(i,1:i-1), 1)

  if(rtmp1 < 1d-12) then
   write(6,'(/,A)') 'ERROR in subroutine cholesky: density matrix not positive &
                    &semidefinite.'
   write(6,'(A)') 'rtmp1 < 1d-12.'
   stop
  end if
  rtmp1 = DSQRT(rtmp1)

  do j = 1, nbf, 1
   rtmp2 = ddot(i-1, new_coeff(j,1:i-1), 1, new_coeff(i,1:i-1), 1)
   new_coeff(j,i) = (P(j,i)-rtmp2)/rtmp1
  end do ! for j
 end do ! for i

 deallocate(P)
 t1 = time()
 write(6,'(A,I0)') 'Decomposition time(sec):', t1-t0
end subroutine cholesky

! perform Boys orbital localization (Jacobian 2*2 rotations) on a set of MOs
subroutine boys(nbf, nif, coeff, mo_dipole, new_coeff)
 implicit none
 integer :: t0, t1, time
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif)
!f2py intent(in) :: coeff
!f2py depend(nbf,nif) :: coeff
 real(kind=8), intent(out) :: new_coeff(nbf,nif)
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: new_coeff
 real(kind=8) :: mo_dipole(3,nif,nif)
!f2py intent(in,copy) :: mo_dipole
!f2py depend(nif) :: mo_dipole

 t0 = time()
 write(6,'(/,A)') 'Boys orbital localization begins:'
 new_coeff = coeff

 if(nif == 1) then
  write(6,'(A)') 'Warning in subroutine boys: only 1 orbital. No rotation.'
  return
 end if

 call serial2by2(nbf, nif, new_coeff, 3, mo_dipole)

 t1 = time()
 write(6,'(A,I0)') 'Localization time(sec):', t1-t0
end subroutine boys

! perform Pipek-Mezey orbital localization (Jacobian 2*2 rotations) on a set of MOs
subroutine pm(nshl, shl2atm, ang, ibas, cart, nbf, nif, coeff, S, pop, new_coeff)
 implicit none
 integer :: i, j, k, natom, t0, t1, time, i1, i2, i3, lwork, liwork
 integer, intent(in) :: nshl, nbf, nif
!f2py intent(in) :: nshl, nbf, nif
 integer, intent(in) :: shl2atm(nshl), ibas(nshl)
!f2py intent(in) :: shl2atm, ibas
!f2py depend(nshl) :: shl2atm, ibas
 integer :: ang(nshl)
!f2py intent(in,copy) :: ang
!f2py depend(nshl) :: ang
 integer, allocatable :: bfirst(:) ! size natom
 ! bfirst: the beginning index of basis func. of each atom
 integer, allocatable :: iwork(:), isuppz(:)
 character(len=*), intent(in) :: pop ! 'mulliken' or 'lowdin'
!f2py intent(in) :: pop
 real(kind=8), intent(in) :: coeff(nbf,nif), S(nbf,nbf)
!f2py intent(in) :: coeff, S
!f2py depend(nbf) :: S
!f2py depend(nbf,nif) :: coeff
 real(kind=8), intent(out) :: new_coeff(nbf,nif)
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nif) :: new_coeff
 real(kind=8) :: ddot
 real(kind=8), allocatable :: gross(:,:,:), SC(:,:), rootS(:,:), e(:), ev(:,:), work(:)
 ! gross: gross population of an orthonormalized MO, (natom,nif,nif)
 ! SC: MATMUL(S,C)
 ! rootS: S^-1/2
 ! e: eigenvalues of S after diagonalization
 ! ev: eigenvectors
 logical, intent(in) :: cart
!f2py intent(in) :: cart

 t0 = time()
 write(6,'(/,A)') 'PM orbital localization begins: using '//TRIM(pop)//' population'
 new_coeff = coeff

 if(nif == 1) then
  write(6,'(A)') 'Warning in subroutine pm: only 1 orbital. No rotation.'
  return
 end if

 if(ANY(ang<0)) then
  write(6,'(/,A)') 'ERROR in subroutine pm: there exists some ang(i)<0.'
  stop
 end if

 if(cart) then
  forall(i = 1:nshl) ang(i) = (ang(i)+1)*(ang(i)+2)/2
 else
  forall(i = 1:nshl) ang(i) = 2*ang(i) + 1
 end if

 j = DOT_PRODUCT(ang, ibas)
 if(j /= nbf) then
  write(6,'(/,A)') 'ERROR in subroutine pm: number of basis function is inconsi&
                   &stent between j and nbf.'
  write(6,'(2(A,I0))') 'j=', j, ', nbf=', nbf
  stop
 end if

 natom = shl2atm(nshl) + 1
 write(6,'(A,I0)') 'natom=', natom
 allocate(bfirst(natom+1), source=0)
 bfirst(1) = 1
 do i = 1, nshl, 1
  j = shl2atm(i) + 2
  bfirst(j) = bfirst(j) + ang(i)*ibas(i)
 end do ! for i

 do i = 2, natom+1, 1
  bfirst(i) = bfirst(i) + bfirst(i-1)
 end do ! for i

 select case(TRIM(pop))
 case('mulliken')
  allocate(SC(nbf,nif), source=0d0)
  call dsymm('L', 'L', nbf, nif, 1d0, S, nbf, coeff, nbf, 0d0, SC, nbf)
  allocate(gross(natom,nif,nif))
 
  do i = 1, nif, 1
   do j = i, nif, 1
    do k = 1, natom, 1
     i1 = bfirst(k); i2 = bfirst(k+1) - 1
     i3 = i2 - i1 + 1
     gross(k,j,i) = ddot(i3, coeff(i1:i2,i), 1, SC(i1:i2,j), 1) + &
                    ddot(i3, coeff(i1:i2,j), 1, SC(i1:i2,i), 1)
     gross(k,j,i) = 0.5d0*gross(k,j,i)
    end do ! for k
   end do ! for j
  end do ! for i

 case('lowdin')
  allocate(e(nbf), ev(nbf, nbf), isuppz(2*nbf))
  allocate(SC(nbf,nbf), source=S) ! use SC to temporarily store S
  lwork = -1
  liwork = -1
  allocate(work(1), iwork(1))
  call dsyevr('V', 'A', 'L', nbf, SC, nbf, 0d0, 0d0, 0, 0, 1d-8, k, e, &
              ev, nbf, isuppz, work, lwork, iwork, liwork, i)
  lwork = CEILING(work(1))
  liwork = iwork(1)
  deallocate(work, iwork)
  allocate(work(lwork), iwork(liwork))
  call dsyevr('V', 'A', 'L', nbf, SC, nbf, 0d0, 0d0, 0, 0, 1d-8, k, e, &
              ev, nbf, isuppz, work, lwork, iwork, liwork, i)
  deallocate(SC, isuppz, work, iwork)

  allocate(rootS(nbf,nbf), source=0d0)
  forall(i = 1:nbf) rootS(i,i) = DSQRT(DABS(e(i)))
  deallocate(e)
  allocate(SC(nbf,nbf), source=0d0) ! use SC to temporarily store US^1/2
  call dsymm('R', 'L', nbf, nbf, 1d0, rootS, nbf, ev, nbf, 0d0, SC, nbf)
  call dgemm('N', 'T', nbf, nbf, nbf, 1d0, SC, nbf, ev, nbf, 0d0, rootS, nbf)
  deallocate(ev, SC)
  allocate(SC(nbf,nif), source=0d0) ! use SC to temporarily store (S^-1/2)C
  call dsymm('L', 'L', nbf, nif, 1d0, rootS, nbf, coeff, nbf, 0d0, SC, nbf)
  deallocate(rootS)
  allocate(gross(natom,nif,nif))

  do i = 1, nif, 1
   do j = i, nif, 1
    do k = 1, natom, 1
     i1 = bfirst(k); i2 = bfirst(k+1)-1
     i3 = i2 - i1 + 1
     gross(k,j,i) = ddot(i3, SC(i1:i2,i), 1, SC(i1:i2,j), 1)
    end do ! for k
   end do ! for j
  end do ! for i

 case default
  write(6,'(/,A)') 'ERROR in subroutine pm: wrong population method provided.'
  write(6,'(A)') "Only 'mulliken' or 'lowdin' supported. But input pop="//pop
  stop
 end select

 forall(i=1:nif,j=1:nif,k=1:natom, (j>i)) gross(k,i,j) = gross(k,j,i)
 deallocate(bfirst, SC)

 call serial2by2(nbf, nif, new_coeff, natom, gross)
 deallocate(gross)

 t1 = time()
 write(6,'(A,I0)') 'Localization time(sec):', t1-t0
end subroutine pm

! perform serial 2-by-2 rotation on given MOs
subroutine serial2by2(nbf, nif, coeff, ncomp, mo_dipole)
 implicit none
 integer :: i, j, k, niter
 integer, intent(in) :: nbf, nif, ncomp
 integer, parameter :: niter_max = 9999
 ! niter_max: max number of iterations

 real(kind=8), parameter :: threshold1 = 1d-8, threshold2 = 1d-6
 ! threshold1: determine whether to rotate (and update MOs, dipole integrals)
 ! threshold2: determine whether rotation/localization converged
 ! if threshold2 is set to 1d-5, the localization is sometimes not thorough
 real(kind=8), parameter :: PI = 4d0*DATAN(1d0)

 real(kind=8), intent(inout) :: coeff(nbf,nif), mo_dipole(ncomp,nif,nif)
 ! for Boys, ncomp = 3
 ! for PM,   ncomp = natom, mo_dipole is actually population matrix

 real(kind=8) :: ddot, rtmp, sum_change
 real(kind=8) :: Aij, Bij, alpha, sin_4a, cos_a, sin_a
 real(kind=8) :: cc, ss, sin_2a, cos_2a
 real(kind=8), allocatable :: dipole(:,:,:), motmp(:,:), vtmp(:,:), vdiff(:)
 ! dipole: for Boys, store updated dipole integrals matrix
 !         for PM,   store updated population matrix
 ! motmp: store 2 MOs
 ! vtmp: array dipole(:,y,x)
 ! vdiff := vtmp(:,1) - vtmp(:,3)

 allocate(dipole(ncomp,nif,2), motmp(nbf,2), vtmp(ncomp,3), vdiff(ncomp))

 ! perform 2*2 rotation
 niter = 0
 do while(niter <= niter_max)
  sum_change = 0d0

  do i = 1, nif-1, 1
   do j = i+1, nif, 1
    vtmp(:,1) = mo_dipole(:,i,i)
    vtmp(:,2) = mo_dipole(:,j,i)
    vtmp(:,3) = mo_dipole(:,j,j)
    vdiff = vtmp(:,1) - vtmp(:,3)
    Aij = ddot(ncomp,vtmp(:,2),1,vtmp(:,2),1) - 0.25d0*ddot(ncomp,vdiff,1,vdiff,1)
    Bij = ddot(ncomp, vdiff, 1, vtmp(:,2), 1)
    rtmp = HYPOT(Aij, Bij)
    sin_4a = Bij/rtmp
    rtmp = rtmp + Aij
    if(rtmp < threshold1) cycle

    sum_change = sum_change + rtmp
    sin_4a = MAX(-1d0, MIN(sin_4a, 1d0)) ! in case of numerical error
    alpha = DASIN(sin_4a) ! [-PI/2,PI/2]
    if(Aij > 0d0) then
     alpha = PI - alpha
    else if(Aij<0d0 .and. Bij<0d0) then
     alpha = 2d0*PI + alpha
    end if
    alpha = 0.25d0*alpha
    cos_a = DCOS(alpha); sin_a = DSIN(alpha)

    ! update two orbitals
    motmp(:,1) = coeff(:,i)
    motmp(:,2) = coeff(:,j)
    coeff(:,i) = cos_a*motmp(:,1) + sin_a*motmp(:,2)
    coeff(:,j) = cos_a*motmp(:,2) - sin_a*motmp(:,1)

    ! update corresponding dipole integrals, only indices in range to be updated
    cc = cos_a*cos_a
    ss = sin_a*sin_a
    sin_2a = 2d0*sin_a*cos_a
    cos_2a = cc - ss
    dipole(:,i,1) = cc*vtmp(:,1) + ss*vtmp(:,3) + sin_2a*vtmp(:,2)
    dipole(:,j,2) = ss*vtmp(:,1) + cc*vtmp(:,3) - sin_2a*vtmp(:,2)
    dipole(:,j,1) = cos_2a*vtmp(:,2) - 0.5d0*sin_2a*vdiff
    dipole(:,i,2) = dipole(:,j,1)

    do k = 1, nif, 1
     if(k==i .or. k==j) cycle
     dipole(:,k,1) = cos_a*mo_dipole(:,k,i) + sin_a*mo_dipole(:,k,j)
     dipole(:,k,2) = cos_a*mo_dipole(:,k,j) - sin_a*mo_dipole(:,k,i)
    end do ! for k

    mo_dipole(:,:,i) = dipole(:,:,1)
    mo_dipole(:,:,j) = dipole(:,:,2)
    mo_dipole(:,i,:) = dipole(:,:,1)
    mo_dipole(:,j,:) = dipole(:,:,2)
   end do ! for j
  end do ! for i

  niter = niter + 1
  write(6,'(A,I5,A,F15.7)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < threshold2) exit
 end do ! for while

 deallocate(motmp, vdiff, vtmp, dipole)

 if(niter <= niter_max) then
  write(6,'(A)') 'Orbital localization converged successfully.'
 else
  write(6,'(/,A)') 'ERROR in subroutine serial2by2: niter_max exceeded.'
  write(6,'(A,I0)') 'niter_max=', niter_max
  stop
 end if
end subroutine serial2by2

! get the value of the modified Boys function
subroutine get_mboys(nif, ncore, npair, nopen, mo_dipole)
 implicit none
 integer :: i, j, nocc
 integer, intent(in) :: nif, ncore, npair, nopen
!f2py intent(in) :: nif, nocc, npair, nopen
 real(kind=8), intent(in) :: mo_dipole(3,nif,nif)
!f2py intent(in) :: mo_dipole
!f2py depend(nif) :: mo_dipole
 real(kind=8) :: ddot, fBoys, temp_dipole(3)

 nocc = ncore + npair + nopen

 fBoys = 0d0
 j = ncore + npair
 do i = 1, j, 1
  temp_dipole =  mo_dipole(1:3,i,i)
  fBoys = fBoys + ddot(3, temp_dipole, 1, temp_dipole, 1)
 end do
 fBoys = DSQRT(fBoys/DBLE(npair))
 write(6,'(A,F13.6)') 'In occ, Modified f(Boys)=', fBoys
 fBoys = 0d0
 j = nocc + npair
 do i = nocc+1, j, 1
  temp_dipole =  mo_dipole(1:3,i,i)
  fBoys = fBoys + ddot(3, temp_dipole, 1, temp_dipole, 1)
 end do
 fBoys = DSQRT(fBoys/DBLE(npair))
 write(6,'(A,F13.6)') 'In vir, Modified f(Boys)=', fBoys
end subroutine get_mboys

! perform immediate Boys localization by diagonalizing the DxDx+DyDy+DzDz
! This was an idea came up with me during 2017 October, and showed in group
!  meeting on 2017 Oct.12. However, I found that this subroutine works poorly
!  since it always converged to saddle points, not minima. This is because at
!  that time I didn't fully deduce the final solution, but put a restriction
!  and deduce this special solution. During 2022 May. 22~28, I deduce the
!  correct solution (see subroutine boys_noiter)
subroutine boys_diag(nbf, nmo, mo_coeff, mo_dipole, new_coeff)
 implicit none
 integer :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8) :: mo_coeff(nbf,nmo), mo_dipole(3,nmo,nmo), new_coeff(nbf,nmo)
!f2py intent(in) :: mo_coeff, mo_dipole
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nmo) :: mo_coeff, new_coeff
!f2py depend(nmo) :: mo_dipole
 real(kind=8), allocatable :: f(:,:), w(:)

 allocate(f(nmo,nmo), source=0d0)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(1,:,:),nmo, mo_dipole(1,:,:),nmo, 0d0,f,nmo)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(2,:,:),nmo, mo_dipole(2,:,:),nmo, 1d0,f,nmo)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(3,:,:),nmo, mo_dipole(3,:,:),nmo, 1d0,f,nmo)

 allocate(w(nmo), source=0d0)
 call diag_get_e_and_vec(nmo, f, w)
 deallocate(w)
 call dgemm('N','N',nbf,nmo,nmo, 1d0,mo_coeff,nbf, f,nmo, 0d0,new_coeff,nbf)
 deallocate(f)
end subroutine boys_diag

! written at the same time as subroutine boys_diag
subroutine solve_boys_lamda_matrix(nbf, nmo, coeff, lo_coeff, mo_dipole)
 implicit none
 integer :: i, j
 integer :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8) :: coeff(nbf,nmo), lo_coeff(nbf,nmo), mo_dipole(3,nmo,nmo)
!f2py intent(in) :: coeff, lo_coeff, mo_dipole
!f2py depend(nbf,nmo) :: coeff, lo_coeff
!f2py depend(nmo) ::  mo_dipole
 real(kind=8), allocatable :: f(:,:), U(:,:), fU(:,:), lamda(:,:)

 allocate(f(nmo,nmo), source=0d0)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(1,:,:),nmo, mo_dipole(1,:,:),nmo, 0d0,f,nmo)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(2,:,:),nmo, mo_dipole(2,:,:),nmo, 1d0,f,nmo)
 call dsymm('L','L',nmo,nmo, 1d0,mo_dipole(3,:,:),nmo, mo_dipole(3,:,:),nmo, 1d0,f,nmo)

 allocate(U(nmo,nmo))
 call get_u(nbf, nmo, coeff, lo_coeff, U) 
 allocate(fU(nmo,nmo), source=0d0)
 call dsymm('L','L',nmo,nmo, 1d0,f,nmo, U,nmo, 0d0,fU,nmo)
 deallocate(f)
 allocate(lamda(nmo,nmo), source=0d0)
 call dgemm('T', 'N', nmo,nmo,nmo, -4d0, U,nmo, fU,nmo, 0d0, lamda,nmo)
 deallocate(fU, U)

 do i = 1, nmo, 1
  write(6,'(20F14.4)') (lamda(j,i),j=1,i)
 end do ! for i
 deallocate(lamda)
end subroutine solve_boys_lamda_matrix

subroutine idx_map_2d(n, k, p, q)
 implicit none
 integer, intent(in) :: n, k
 integer, intent(out) :: p, q

 p = k/n
 if(k - n*p > 0) p = p+1
 q = k - (p-1)*n
end subroutine idx_map_2d

! A non-iterative Boys orbital localization solver
subroutine boys_noiter(nbf, nmo, mo_coeff, mo_dipole, new_coeff)
 implicit none
 integer :: i, j, k, m, p, q, nmo_s
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8) :: rtmp(3)
 real(kind=8), intent(in) :: mo_coeff(nbf,nmo), mo_dipole(3,nmo,nmo)
!f2py intent(in) :: mo_coeff, mo_dipole
!f2py depend(nbf,nmo) :: mo_coeff
!f2py depend(nmo) :: mo_dipole
 real(kind=8), intent(out) :: new_coeff(nbf,nmo)
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nmo) :: new_coeff
 real(kind=8), allocatable :: e(:), u(:,:), y(:,:), vtmp(:)

 if(nmo == 1) then
  new_coeff = mo_coeff
  return
 else
  new_coeff = 0d0 ! initialization
 end if

 nmo_s = nmo*nmo
 allocate(y(nmo_s,nmo_s), source=0d0)

 ! diagonal elements
 do i = 1, nmo_s, 1
  call idx_map_2d(nmo, i, k, m)
  rtmp = mo_dipole(:,m,k)
  y(i,i) = DOT_PRODUCT(rtmp, rtmp)
 end do ! for i

 ! non-diagonal elements
 do i = 1, nmo_s-1, 1
  call idx_map_2d(nmo, i, k, m)
  rtmp = mo_dipole(:,m,k)
  do j = i+1, nmo_s, 1
   call idx_map_2d(nmo, j, p, q)
   y(j,i) = DOT_PRODUCT(rtmp, mo_dipole(:,q,p))
   y(i,j) = y(j,i)
  end do ! for j
 end do ! for i

 allocate(e(nmo_s)) ! ((X')^T)Y(X') = (-1/2)(\epsilon)
 call diag_get_e_and_vec(nmo_s, y, e)
 write(6,'(5F18.12)') e
 stop
 ! interestingly, there ae only 3 non-zero eigenvalues in the array e, maybe
 ! because the dipole moment has only 3 components

 ! reverse the oerder of eigenvalues and eigenvectors
 allocate(vtmp(nmo_s), source=e)
 forall(i = 1:nmo_s) e(i) = vtmp(nmo_s-i+1)
 do i = 1, nmo_s/2, 1
  vtmp = y(:,i); y(:,i) = y(:,nmo_s-i+1); y(:,nmo_s-i+1) = vtmp
 end do ! for i
 deallocate(vtmp)

 ! solve the unitary matrix U
 allocate(u(nmo,nmo), source=0d0)
 do i = 1, nmo, 1
  do j = 1, nmo, 1
   rtmp(1) = y(j,(i-1)*nmo+i)
   if(rtmp(1) < -1d-5) then
    write(6,'(A)') 'ERROR in subroutine boys_noiter: too negative.'
    write(6,'(A,F16.8)') 'rtmp(1)=', rtmp(1)
    stop
   else if(rtmp(1) > 0d0) then
    u(j,i) = DABS(rtmp(1))
   end if
   ! for -1d-5 < rtmp(1) < 0d0, u(j,i) is set to zero (the initialized value)
  end do ! for j
 end do ! for i

 ! set u(:,1) > 0, determine the signs of other elements in matrix u
 do i = 2, nmo, 1
  do j = 1, nmo, 1
   if(y(j,i) < 0d0) u(j,i) = -u(j,i)
  end do ! for j
 end do ! for i

 deallocate(y)
 allocate(y(nmo,nmo), source=0d0)
 y = MATMUL(u, TRANSPOSE(u))
 forall(i = 1:nmo) y(i,i) = y(i,i) - 1d0
 y = DABS(y)
 write(6,'(A,F12.6)') 'max(U)=', MAXVAL(u)
 write(6,'(A,F12.6)') 'sum(U(U^T) - I)=', SUM(y)
 stop

! y = MATMUL(TRANSPOSE(u), u)
! forall(i = 1:nmo) y(i,i) = y(i,i) - 1d0
! y = DABS(y)
! write(6,'(A,F12.6)') 'max(U)=', MAXVAL(u)
! write(6,'(A,F12.6)') 'sum((U^T)U - I)=', SUM(y)

 deallocate(y)
 call dgemm('N','N',nbf,nmo,nmo, 1d0,mo_coeff,nbf, u,nmo, 0d0, new_coeff,nbf)
 deallocate(u)
end subroutine boys_noiter

