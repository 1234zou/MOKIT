! written by jxzou at 20190723
! updated by jxzou at 20200411: add Pipek-Mezey orbital localization (DOI: 10.1063/1.456588)
! updated by jxzou at 20200413: add Cholesky decomposition LMOs (DOI: 10.1063/1.2360264)
! updated by jxzou at 20220520: generate NO from a NSO .fch file

! Note: before PySCF-1.6.4, its dumped .molden file is wrong when using Cartesian functions.

! localize singly occupied orbitals in a .fch file
subroutine localize_singly_occ_orb(fchname, pm_loc)
 implicit none
 integer :: na, nb
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, intent(in) :: pm_loc
!f2py intent(in) :: pm_loc

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
!f2py intent(in) :: i1, i2
 character(len=240) :: lmofch, pyname, outname
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname
 logical, intent(in) :: pm_loc
!f2py intent(in) :: pm_loc

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
 write(fid,'(A)') 'from mokit.lib.gaussian import BOHR2ANG, load_mol_from_fch'
 write(fid,'(A)') 'from mokit.lib.rwwfn import read_nbf_and_nif_from_fch, \'
 write(fid,'(A)') ' read_eigenvalues_from_fch'
 write(fid,'(A)') 'from mokit.lib.fch2py import fch2py'
 write(fid,'(A)') 'from mokit.lib.py2fch import py2fch'
 if(.not. pm_loc) then
  write(fid,'(A)') 'from mokit.lib.gaussian import get_ao_dip'
 end if
 write(fid,'(A,/)') 'from mokit.lib.auto import loc_ini_guess, loc_driver'
 write(fid,'(A)') 'import os'

 write(fid,'(/,A)') "fchname = '"//TRIM(fchname)//"'"
 write(fid,'(A)') "lmofch = '"//TRIM(lmofch)//"'"
 write(fid,'(A)') 'mol = load_mol_from_fch(fchname)'
 write(fid,'(A)') 'nbf, nif = read_nbf_and_nif_from_fch(fchname)'
 write(fid,'(A)') "mo_coeff = fch2py(fchname, nbf, nif, 'a')"
 write(fid,'(2(A,I0),A)') 'idx = range(',i1-1,',',i2,')'
 write(fid,'(A)') 'nmo = len(idx)'
 write(fid,'(A)') 'nmo1, lmo_ini = loc_ini_guess(mol, mo_coeff[:,idx], nmo)'
 if(pm_loc) then
  write(fid,'(A)') "loc_orb = loc_driver(mol, lmo_ini, nmo, method='pm')"
 else
  write(fid,'(A)') 'center, ao_dip = get_ao_dip(mol)'
  write(fid,'(A)') 'ao_dip = ao_dip*BOHR2ANG'
  write(fid,'(A)') "loc_orb = loc_driver(mol, lmo_ini, nmo, method='boys', ao_d&
                   &ip=None)"
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
  write(6,'(/,A)') 'ERROR in subroutine gen_ao_dm_from_fch: itype out of range.'
  write(6,'(A,I0)') 'Only 0~4 are allowed. But got itype=', itype
  stop
 end select

 allocate(n(nif,nif),source=0d0)
 forall(i = 1:nif) n(i,i) = noon(i)
 deallocate(noon)

 allocate(mo(nbf,nif))
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
 call gen_no_from_dm_and_ao_ovlp(nbf, nif, dm, S, noon, mo)
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
 integer, intent(in) :: nbf, nif
!f2py intent(in) :: nbf, nif
 real(kind=8), intent(in) :: coeff(nbf,nif), S(nbf,nbf), P(nbf,nbf)
!f2py intent(in) :: coeff, S, P
!f2py depend(nbf,nif) :: coeff
!f2py depend(nbf) :: S, P
 real(kind=8), intent(out) ::dm(nif,nif)
!f2py intent(out) :: dm
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
!f2py intent(in) :: ndb, nopen
 ! ndb and nopen cannot be determined from .dat file, manual input required
 real(kind=8) :: a, a2, fac
 real(kind=8), allocatable :: coeff(:,:), ci_coeff(:,:), rtmp(:,:)
 ! coeff: MO coefficients
 ! ci_coeff: GVB CI coefficients, pair coefficients
 ! rtmp: temporary array to hold two orbitals
 character(len=240), intent(in) :: datname
!f2py intent(in) :: datname

 if(ndb<0 .or. nopen<0) then
  write(6,'(/,A)') 'ERROR in subroutine gen_cf_orb: ndb<0 or nopen<0 found.'
  write(6,'(A)') 'Correct values should be assigned to these two parameters.'
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
  write(6,'(/,A)') 'ERROR in subroutine gen_cf_orb: pair coefficients in file '&
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
  write(6,'(/,A)') 'ERROR in subroutine gen_cf_orb: probably wrong ndb or nopen.'
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

! perform Cholesky factorization with complete pivoting to find Cholesky LMOs
subroutine cholesky_dm(nbf, nif, dm, mo)
 implicit none
 integer :: i, j, rank
 integer, intent(in) :: nbf, nif
 integer, allocatable :: piv(:)
 real(kind=8), intent(inout) :: dm(nbf,nbf)
 real(kind=8), intent(out) :: mo(nbf,nif)
 real(kind=8), allocatable :: work(:), piv2(:,:)

 allocate(piv(nbf), source=0)
 allocate(work(2*nbf), source=0d0)
 call dpstrf('L', nbf, dm, nbf, piv, rank, -1d0, work, i)
 deallocate(work)
 ! Note that i>0 is possible since rank<<nbf usually
 if(i < 0) then
  write(6,'(/,A)') 'ERROR in subroutine cholesky_dm: dpstrf failed.'
  write(6,'(A,4I7)') 'nbf, nif, rank, info=', nbf, nif, rank, i
  stop
 end if
 if(rank < nif) then
  write(6,'(/,A)') 'ERROR in subroutine cholesky_dm: rank<nif.'
  write(6,'(A,3I7)') 'nbf, nif, rank=', nbf, nif, rank
  stop
 end if

!$omp parallel do schedule(dynamic) default(shared) private(i,j)
 do i = 2, nbf, 1
  do j = 1, i-1, 1
   dm(j,i) = 0d0
  end do ! for j
 end do ! for i
!$omp end parallel do

 if(rank < nbf) dm(:,rank+1:nbf) = 0d0
 allocate(piv2(nbf,nbf), source=0d0)
 forall(i = 1:nbf) piv2(piv(i),i) = 1d0
 deallocate(piv)

 mo = 0d0
 call dgemm('N','N', nbf,nif,nbf, 1d0,piv2,nbf, dm(:,1:nif),nbf, 0d0,mo,nbf)
 deallocate(piv2)
end subroutine cholesky_dm

! Use Cholesky factorization/decomposition of the density matrix to generate
! LMOs. These LMOs are usually less localized than Boys/PM localized ones, but
! they can be good initial guess of Boys/PM. A set of MOs and a unitary matrix
! will be returned, i.e. new_mo = mo*u
subroutine cholesky_mo(nbf, nif, mo)
 implicit none
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(inout) :: mo(nbf,nif)
 real(kind=8), allocatable :: dm(:,:)

 if(nif == 1) then
  write(6,'(/,A)') 'Warning from subroutine cholesky_mo: only 1 orbital. No rot&
                   &ation.'
  return
 end if

 allocate(dm(nbf,nbf)) ! P = C(C^T)
 call calc_cct(nbf, nif, mo, dm)

 ! generate Cholesky LMOs
 call cholesky_dm(nbf, nif, dm, mo)
 deallocate(dm)
end subroutine cholesky_mo

! Use Cholesky factorization/decomposition of the density matrix to generate
! LMOs. A slightly difference with subroutine cholesky_mo above is that this
! subroutine use MO coefficients expanded on the symmetrically/canonically
! orthogonalized MOs.
subroutine cholesky_mo2(nbf, nif, ao_ovlp, mo)
 implicit none
 integer :: k
 integer, intent(in) :: nbf, nif
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
 real(kind=8), intent(inout) :: mo(nbf,nif)
 real(kind=8), allocatable :: orth_mo(:,:), v(:,:)

 if(nif == 1) then
  write(6,'(/,A)') 'Warning from subroutine cholesky_mo2: only 1 orbital. No ro&
                   &tation.'
  return
 end if

 ! Generate a set of orthonormalized MOs from the AO overlap. Basis set linear
 ! dependency is possible, so k <= nbf.
 allocate(orth_mo(nbf,nbf))
 call gen_ortho_mo(nbf, ao_ovlp, .true., k, orth_mo)
 write(6,'(2(A,I0))') 'nbf=', nbf, ', k=', k

 ! calculate the unitary matrix V in equation `mo = orth_mo*V`
 ! V = (orth_mo^T)S(mo)
 allocate(v(k,nif))
 call calc_CTSCp2(nbf, k, nif, orth_mo(:,1:k), ao_ovlp, mo, v)

 call cholesky_mo(k, nif, v)

 ! Rotate v back to AO basis
 mo = 0d0
 call dgemm('N','N', nbf,nif,k, 1d0, orth_mo(:,1:k), nbf, v, k, 0d0, mo, nbf)
 deallocate(orth_mo, v)
end subroutine cholesky_mo2

! Make a set of given MOs resembles (some) atomic orbitals. Learning from atomic
!  initial guess in `pyscf/lo/boys.py`.
! According to jxzou's tests, the locality of obtained LMOs are usually atomic >
!  Cho2 > Cho (but not always), where Cho2 means Cholesky LMOs based on orthogonal
!  AO basis.
subroutine resemble_ao(nbf, nmo, ao_ovlp, pseudo_inv, mo)
 implicit none
 integer :: i, k
 integer, intent(in) :: nbf, nmo
 integer, allocatable :: idx(:)
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
 real(kind=8), intent(inout) :: mo(nbf,nmo)
 real(kind=8), allocatable :: orth_mo(:,:), u0(:,:), u1(:,:), u(:,:), vt(:,:), &
  uvt(:,:), norm(:), s(:), d(:)
 logical, intent(in) :: pseudo_inv

 if(nmo == 1) then
  write(6,'(A)') 'Warning from subroutine resemble_ao: only 1 orbital. No need &
                 &to change.'
  return
 end if

 allocate(orth_mo(nbf,nbf))
 call gen_ortho_mo(nbf, ao_ovlp, pseudo_inv, k, orth_mo)
 if(nmo > k) then
  write(6,'(/,A)') 'ERROR in subroutine resemble_ao: nmo>k. Required number of &
                   &MOs is larger'
  write(6,'(A)') 'than that of linear independent MOs.'
  stop
 end if

 ! calculate the unitary matrix U in equation `mo = orth_mo*U_0`
 ! U_0 = (orth_mo^T)S(mo)
 allocate(u0(k,nmo))
 call calc_CTSCp2(nbf, k, nmo, orth_mo(:,1:k), ao_ovlp, mo, u0)
 deallocate(orth_mo)

 ! find the AOs which has largest overlap with input MOs
 allocate(norm(k))
 forall(i = 1:k) norm(i) = DOT_PRODUCT(u0(i,:), u0(i,:))
 allocate(idx(k))
 call sort_dp_array(k, norm, .false., idx)
 deallocate(norm)
 allocate(u1(nmo,nmo))
 forall(i = 1:nmo) u1(i,:) = u0(idx(i),:)
 deallocate(idx, u0)

 ! compute the rotation matrix of original MOs
 allocate(u(nmo,nmo), vt(nmo,nmo), s(nmo))
 call do_svd(nmo, nmo, u1, u, vt, s)
 deallocate(u1)
 allocate(d(nmo))
 call calc_usut_diag_elem(nmo, s, u, d)
 deallocate(s)
 allocate(idx(nmo))
 call sort_dp_array(nmo, d, .false., idx)
 deallocate(d)

 ! U(V^T)
 allocate(uvt(nmo,nmo), source=0d0)
 call dgemm('N', 'N', nmo, nmo, nmo, 1d0, u, nmo, vt, nmo, 0d0, uvt, nmo)
 deallocate(u, vt)

 ! MO*(V(U^T))
 allocate(u0(nbf,nmo), source=0d0)
 call dgemm('N','T', nbf, nmo, nmo, 1d0, mo, nbf, uvt, nmo, 0d0, u0, nbf)
 deallocate(uvt)

 ! MOs in u0 are not necessarily in overlap descending order, sorting using
 ! idx is needed.
 forall(i = 1:nmo) mo(:,i) = u0(:,idx(i))
 deallocate(idx, u0)
end subroutine resemble_ao

! Combining AO-resemblance with Cholesky LMO techniques. The Cholesky LMO
! technique does not always make target orbitals spatially localized. This
! subroutine remains to be improved.
subroutine combined_ao_cholesky(nbf, nmo, ao_ovlp, pseudo_inv, mo)
 implicit none
 integer :: i, j, k, nlmo
 integer, intent(in) :: nbf, nmo
 integer, allocatable :: idx(:)
 real(kind=8), parameter :: ovlp_thres = 0.56d0
 ! threshold for selecting MOs transformed to Cholesky LMOs
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
 real(kind=8), intent(inout) :: mo(nbf,nmo)
 real(kind=8), allocatable :: orth_mo(:,:), u0(:,:), u1(:,:), u(:,:), vt(:,:),&
  uvt(:,:), norm(:), s(:), d(:)
 logical, intent(in) :: pseudo_inv

 if(nmo == 1) then
  write(6,'(A)') 'Warning from subroutine combined_ao_cholesky: only 1 orbital.&
                 & No need to change.'
  return
 end if

 allocate(orth_mo(nbf,nbf))
 call gen_ortho_mo(nbf, ao_ovlp, pseudo_inv, k, orth_mo)
 if(nmo > k) then
  write(6,'(/,A)') 'ERROR in subroutine combined_ao_cholesky: nmo>k. Required n&
                   &umber of MOs'
  write(6,'(A)') 'is larger than that of linear independent MOs.'
  stop
 end if

 ! calculate the unitary matrix U in equation `mo = orth_mo*U_0`
 ! U_0 = (orth_mo^T)S(mo)
 allocate(u0(k,nmo))
 call calc_CTSCp2(nbf, k, nmo, orth_mo(:,1:k), ao_ovlp, mo, u0)

 ! find the AOs which has largest overlap with input MOs
 allocate(norm(k))
 forall(i = 1:k) norm(i) = DOT_PRODUCT(u0(i,:), u0(i,:))
 allocate(idx(k))
 call sort_dp_array(k, norm, .false., idx)
 deallocate(norm)
 allocate(u1(nmo,nmo))
 forall(i = 1:nmo) u1(i,:) = u0(idx(i),:)
 deallocate(idx)

 ! compute the rotation matrix of original MOs
 allocate(u(nmo,nmo), vt(nmo,nmo), s(nmo))
 call do_svd(nmo, nmo, u1, u, vt, s)
 deallocate(u1)
 allocate(d(nmo))
 call calc_usut_diag_elem(nmo, s, u, d)
 deallocate(s)
 allocate(idx(nmo))
 call sort_dp_array(nmo, d, .false., idx)
 nlmo = COUNT(d > ovlp_thres)
 deallocate(d)
 write(6,'(A,I0)') 'nlmo=', nlmo
 ! If nlmo=nmo, all orbitals match well with AOs.
 ! If 0< nlmo <nmo, some orbitals match well with AOs, for remaining orbitals
 ! we generate their Cholesky LMOs.
 ! If nlmo=0, no orbital match well with AOs, generate Cholesky LMOs for all
 ! orbitals.

 if(nlmo == 0) then
  deallocate(u, vt)
  call cholesky_mo(k, nmo, u0)
  mo = 0d0
  call dgemm('N','N', nbf,nmo,k, 1d0,orth_mo(:,1:k),nbf, u0,k, 0d0,mo,nbf)
  deallocate(u0)
 else
  ! U(V^T)
  deallocate(u0)
  allocate(uvt(nmo,nmo), source=0d0)
  call dgemm('N', 'N', nmo, nmo, nmo, 1d0, u, nmo, vt, nmo, 0d0, uvt, nmo)
  deallocate(u, vt)
  ! MO*(V(U^T))
  allocate(u1(nbf,nmo), source=0d0)
  call dgemm('N','T', nbf, nmo, nmo, 1d0, mo, nbf, uvt, nmo, 0d0, u1, nbf)
  deallocate(uvt)
  ! MOs in u1 are not necessarily in overlap descending order, sorting using
  ! idx is needed.
  forall(i = 1:nmo) mo(:,i) = u1(:,idx(i))
  deallocate(idx, u1)
  if(nlmo < nmo-1) then ! generate Cholesky LMOs for remaining orbitals
   j = nmo - nlmo
   allocate(u0(k,j))
   allocate(u1(nbf,j), source=mo(:,nlmo+1:nmo))
   call calc_CTSCp2(nbf, k, j, orth_mo(:,1:k), ao_ovlp, u1, u0)
   call cholesky_mo(k, j, u0)
   u1 = 0d0
   call dgemm('N','N', nbf,j,k, 1d0,orth_mo(:,1:k),nbf, u0,k, 0d0, u1, nbf)
   mo(:,nlmo+1:nmo) = u1
   deallocate(u0, u1)
  end if
 end if

 deallocate(orth_mo)
end subroutine combined_ao_cholesky

! Generate initial guess of orbital localization. If all atoms are chosen, the
! subroutine resemble_ao will be caled. Otherwise AO-like LMOs centered on
! specified atoms are constructed.
! new_mo is the unitary transformation of mo.
subroutine gen_loc_ini_guess(natom, nbf, nif, chosen, bfirst, ao_ovlp, mo, &
                             nif1, new_mo)
 implicit none
 integer :: i, j, k, nbf1
 integer, intent(in) :: natom, nbf, nif
!f2py intent(in) :: natom, nbf, nif
 integer, intent(in) :: bfirst(natom+1)
!f2py intent(in) :: bfirst
!f2py depend(natom) :: bfirst
 integer, intent(out) :: nif1
!f2py intent(out) :: nif1
 integer, allocatable :: idx(:)
 real(kind=8), parameter :: sv_thres = 0.92d0
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
!f2py intent(in) :: ao_ovlp
!f2py depend(nbf,nbf) :: ao_ovlp
 real(kind=8), intent(in) :: mo(nbf,nif)
!f2py intent(in) :: mo
!f2py depend(nbf,nif) :: mo
 real(kind=8), intent(out) :: new_mo(nbf,nif)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nif) :: new_mo
 real(kind=8), allocatable :: ao_ovlp1(:,:), ao_ovlp2(:,:), orth_mo(:,:), &
  mo_ovlp(:,:), proj_mo(:,:), u(:,:), vt(:,:), s(:)
 logical, intent(in) :: chosen(natom)
!f2py intent(in) :: chosen
!f2py depend(natom) :: chosen

 write(6,'(/,A)') 'Construct AO-like LMOs as initial guess...'

 if(ALL(chosen .eqv. .true.)) then
  nif1 = nif
  new_mo = mo
  call resemble_ao(nbf, nif, ao_ovlp, .true., new_mo)
  write(6,'(A)') 'Done construction.'
  return
 end if

 if(nif == 1) then
  write(6,'(A)') 'Warning from subroutine gen_loc_ini_guess: only 1 orbital. No&
                 & need to change.'
  return
 end if

 write(6,'(A)') 'Not all atoms are chosen. Only localized orbitals centered on &
                &specified atoms to'
 write(6,'(A)') 'be solved...'
 allocate(idx(nbf), source=0)
 nbf1 = 0

 do i = 1, natom, 1
  if(chosen(i)) then
   k = bfirst(i+1) - bfirst(i)
   forall(j = 1:k) idx(nbf1+j) = bfirst(i)+j-1
   nbf1 = nbf1 + k
  end if
 end do ! for i

 allocate(ao_ovlp1(nbf1,nbf1), ao_ovlp2(nbf,nbf1))

!$omp parallel do schedule(dynamic) default(shared) private(i,j,k)
 do i = 1, nbf1, 1
  k = idx(i)
  do j = 1, i, 1
   ao_ovlp1(j,i) = ao_ovlp(idx(j),k)
  end do ! for j
  ao_ovlp2(:,i) = ao_ovlp(:,k)
 end do ! for i
!$omp end parallel do

 call symmetrize_dmat(nbf1, ao_ovlp1)
 deallocate(idx)

 ! generate OAOs based on basis functions of specified atoms
 allocate(orth_mo(nbf1,nbf1))
 call gen_ortho_mo(nbf1, ao_ovlp1, .true., k, orth_mo)

 ! calculate matrix products (A^T)BC, i.e. MO cross overlap matrix
 ! ((C')^T)SC
 allocate(mo_ovlp(nif,k))
 call calc_atbc(nif, nbf, nbf1, k, mo, ao_ovlp2, orth_mo(:,1:k), mo_ovlp)

 ! perform SVD on the cross overlap matrix
 allocate(u(nif,nif), vt(k,k), s(nif))
 call do_svd(nif, k, mo_ovlp, u, vt, s)
 nif1 = COUNT(s > sv_thres)
 write(6,'(A)') 'Singular values:'
 write(6,'(5(1X,ES15.8))') s
 deallocate(mo_ovlp, u, s)

 ! solve reference MOs proj_mo(:,1:nif1)
 allocate(proj_mo(nbf1,k), source=0d0)
 call dgemm('N','T',nbf1,k,k,1d0,orth_mo(:,1:k),nbf1,vt,k,0d0,proj_mo,nbf1)
 deallocate(orth_mo, vt)
 call resemble_ao(nbf1, nif1, ao_ovlp1, .true., proj_mo(:,1:nif1))
 deallocate(ao_ovlp1)

 ! Make new_mo(:,1:nif1) resembles proj_mo(:,1:nif1).
 ! new_mo(:,nif1+1:) are remaining orbitals with non-zero MO coefficients, and
 ! will probably be unused in subsequent orbital localization.
 call orb_resemble_ref1(nbf, nif, mo, nbf1, nif1, proj_mo, ao_ovlp2, new_mo)
 deallocate(proj_mo, ao_ovlp2)
 write(6,'(A)') 'Done construction.'
end subroutine gen_loc_ini_guess

subroutine classify_lmo(nmo, natom, gross, conn)
 implicit none
 integer :: i, j, n, npair
 integer, intent(in) :: nmo, natom
!f2py intent(in) :: nmo, natom
 integer, allocatable :: ijmap(:,:)
 integer, intent(out) :: conn(nmo,nmo)
!f2py intent(out) :: conn
!f2py depend(nmo) :: conn
 real(kind=8), parameter :: thres = 1d-6
 real(kind=8) :: ddot, rtmp, Aij, Bij
 real(kind=8), intent(in) :: gross(natom,nmo,nmo)
!f2py intent(in) :: gross
!f2py depend(natom,nmo) :: gross
 real(kind=8), allocatable :: vtmp(:,:), vdiff(:)

 conn = 0
 npair = nmo*(nmo-1)/2
 allocate(ijmap(2,npair))
 call get_triu_idx1(nmo, ijmap)
 allocate(vtmp(natom,3), vdiff(natom))

 do n = 1, npair, 1
  i = ijmap(1,n); j = ijmap(2,n)
  vtmp(:,1) = gross(:,i,i)
  vtmp(:,2) = gross(:,i,j)
  vtmp(:,3) = gross(:,j,j)
  vdiff = vtmp(:,1) - vtmp(:,3)
  Aij = ddot(natom,vtmp(:,2),1,vtmp(:,2),1) - 0.25d0*ddot(natom,vdiff,1,vdiff,1)
  Bij = ddot(natom, vdiff, 1, vtmp(:,2), 1)
  rtmp = HYPOT(Aij, Bij) + Aij
  if(rtmp > thres) then
   conn(j,i) = 1; conn(i,j) = 1
  end if
 end do ! for n

 deallocate(ijmap)
end subroutine classify_lmo

! Perform Berry orbital localization (Jacobian 2*2 rotations) on a set of MOs.
! This subroutine is used for PBC gamma-point orbital localization, and cannot
! be used for an isolated molecule.
! The input ao_zdip must be in unit Angstrom.
subroutine berry(natom, nbf, nmo, maxcyc, lo_diis, bfirst, dis, mo, ao_zdip, &
                 dis_tol, conv_tol, new_mo)
 use lo_info, only: max_niter, dis_thres, conv_thres, npair, ijmap, diis, u, &
  cayley_k, k_old, k_diis, nfile, binfile
 implicit none
 integer :: k
 integer, intent(in) :: natom, nbf, nmo, maxcyc
!f2py intent(in) :: natom, nbf, nmo, maxcyc
 integer, intent(in) :: bfirst(natom+1)
!f2py intent(in) :: bfirst
!f2py depend(natom) :: bfirst
 real(kind=8), intent(in) :: dis(natom,natom), mo(nbf,nmo), dis_tol, conv_tol
!f2py intent(in) :: dis, mo, dis_tol, conv_tol
!f2py depend(natom) :: dis
!f2py depend(nbf,nmo) :: mo
 real(kind=8), intent(out) :: new_mo(nbf,nmo)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nmo) :: new_mo
 complex(kind=8), intent(in) :: ao_zdip(3,nbf,nbf) ! complex symmetric
!f2py intent(in) :: ao_zdip
!f2py depend(nbf) :: ao_zdip
 complex(kind=8), allocatable :: mo_zdip(:,:,:)
 character(len=240) :: proname
 logical, intent(in) :: lo_diis
!f2py intent(in) :: lo_diis

 write(6,'(/,A)') 'Berry orbital localization begins:'
 new_mo = mo; dis_thres = dis_tol; conv_thres = conv_tol
 diis = lo_diis; max_niter = maxcyc

 if(nmo == 1) then
  write(6,'(A)') 'Warning from subroutine berry: only 1 orbital. No rotation.'
  return
 end if

 write(6,'(/,A)') 'Transform AO complex dipole integrals to MO ones...'
 allocate(mo_zdip(3,nmo,nmo))
 call ao2mo_zdip(nbf, nmo, new_mo, ao_zdip, mo_zdip)
 write(6,'(A)') 'Done update.'

 npair = nmo*(nmo-1)/2
 allocate(ijmap(2,npair))
 call get_triu_idx1(nmo, ijmap)

 if(diis) then
  call get_a_random_int(k)
  write(proname,'(A,I0)') 'berry', k
  call init_lo_diis(proname, nmo)
 end if

 if(nmo < 10) then
  call serial2by2_cmplx(nbf, nmo, new_mo, mo_zdip)
 else
  call serial22berry(natom, nbf, nmo, bfirst, dis, new_mo, mo_zdip)
 end if
 !else if(nmo < 500) then
 ! call serial22berry(natom, nbf, nmo, bfirst, dis, new_mo, mo_zdip)
 !else
 ! call para22berry(natom, nbf, nmo, bfirst, dis, new_mo, mo_zdip)
 !end if

 ! job accomplished, no need to update the lower triangle part of mo_zdip
 deallocate(mo_zdip, ijmap)
 if(diis) then
  call delete_files(nfile, binfile)
  deallocate(u, cayley_k, k_old, k_diis, binfile)
 end if
end subroutine berry

! Perform Boys orbital localization (Jacobian 2*2 rotations) on a set of MOs.
! The input ao_dip must be in unit Angstrom.
subroutine boys(natom, nbf, nmo, bfirst, dis, mo, ao_dip, dis_tol, conv_tol, &
                new_mo)
 use lo_info, only: dis_thres, conv_thres, npair, ijmap
 implicit none
 integer, intent(in) :: natom, nbf, nmo
!f2py intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
!f2py intent(in) :: bfirst
!f2py depend(natom) :: bfirst
 real(kind=8), intent(in) :: dis(natom,natom), mo(nbf,nmo), ao_dip(3,nbf,nbf), &
                             dis_tol, conv_tol
!f2py intent(in) :: dis, mo, ao_dip, dis_tol, conv_tol
!f2py depend(natom) :: dis
!f2py depend(nbf,nmo) :: mo
!f2py depend(nbf) :: ao_dip
 real(kind=8), intent(out) :: new_mo(nbf,nmo)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nmo) :: new_mo
 real(kind=8), allocatable :: mo_dip(:,:,:) !, u(:,:), save_mo(:,:)

 write(6,'(/,A)') 'Boys orbital localization begins:'
 new_mo = mo; dis_thres = dis_tol; conv_thres = conv_tol

 if(nmo == 1) then
  write(6,'(A)') 'Warning from subroutine boys: only 1 orbital. No rotation.'
  return
 end if

 npair = nmo*(nmo-1)/2
 allocate(ijmap(2,npair))
 call get_triu_idx1(nmo, ijmap)

 write(6,'(/,A)') 'Transform AO dipole integrals to MO ones...'
 allocate(mo_dip(3,nmo,nmo))
 call ao2mo_dip(nbf, nmo, new_mo, ao_dip, mo_dip)
 write(6,'(A)') 'Done update.'

 if(nmo < 10) then
  call serial2by2(nbf, nmo, 3, new_mo, mo_dip)
 else if(nmo < 500) then
  call serial22boys(natom, nbf, nmo, bfirst, dis, new_mo, mo_dip)
 else
  call para22boys(natom, nbf, nmo, bfirst, dis, new_mo, mo_dip)
 end if

 !call ao2mo_dip(nbf, nmo, new_mo, ao_dip, mo_dip)
 !allocate(u(nmo,nmo))
 !call boys_polar(nmo, mo_dip, u)
 !allocate(save_mo(nbf,nmo), source=new_mo)
 !new_mo = 0d0
 !call dgemm('N','N',nbf,nmo,nmo,1d0,save_mo,nbf,u,nmo,0d0,new_mo,nbf)
 !deallocate(u, save_mo)

 ! job accomplished, no need to update the lower triangle part of mo_dip
 deallocate(mo_dip, ijmap)
end subroutine boys

! perform Pipek-Mezey orbital localization (Jacobian 2*2 rotations) on a set of MOs
subroutine pm(natom, nbf, nmo, bfirst, dis, mo, ao_ovlp, popm, dis_tol, &
              conv_tol, new_mo)
 use lo_info, only: dis_thres, conv_thres, npair, ijmap
 implicit none
 integer, intent(in) :: natom, nbf, nmo
!f2py intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
!f2py intent(in) :: bfirst
!f2py depend(natom) :: bfirst
 real(kind=8), intent(in) :: dis(natom,natom), mo(nbf,nmo), ao_ovlp(nbf,nbf),&
                             dis_tol, conv_tol
!f2py intent(in) :: dis, mo, ao_ovlp, dis_tol, conv_tol
!f2py depend(natom) :: dis
!f2py depend(nbf) :: ao_ovlp
!f2py depend(nbf,nmo) :: mo
 real(kind=8), intent(out) :: new_mo(nbf,nmo)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nmo) :: new_mo
 real(kind=8), allocatable :: gross(:,:,:) ! size (natom,nmo,nmo)
 character(len=8), intent(in) :: popm ! 'mulliken'/'lowdin'
!f2py intent(in) :: popm

 write(6,'(/,A)') 'PM orbital localization begins: using '//TRIM(popm)//&
                  ' population'
 new_mo = mo; dis_thres = dis_tol; conv_thres = conv_tol

 if(nmo == 1) then
  write(6,'(A)') 'Warning from subroutine pm: only 1 orbital. No rotation.'
  return
 end if

 npair = nmo*(nmo-1)/2
 allocate(ijmap(2,npair))
 call get_triu_idx1(nmo, ijmap)

 write(6,'(/,A)') 'Construct gross matrix...'
 allocate(gross(natom,nmo,nmo))
 call calc_gross_pop(natom, nbf, nmo, bfirst, ao_ovlp, new_mo, popm, gross)
 write(6,'(A)') 'Done construction.'

 if(nmo < 10) then
  call serial2by2(nbf, nmo, natom, new_mo, gross)
 else if(nmo < 500) then
  call serial22pm(natom, nbf, nmo, dis, new_mo, gross)
 else
  call para22pm(natom, nbf, nmo, dis, new_mo, gross)
 end if

 ! job accomplished, no need to update the lower triangle part of gross
 deallocate(gross, ijmap)
end subroutine pm

! perform serial 2-by-2 rotation on given MOs
subroutine serial2by2_cmplx(nbf, nmo, coeff, mo_dipole)
 use lo_info, only: npair, ijmap, max_niter, QPI, HPI, upd_thres, conv_thres
 implicit none
 integer :: i, j, k, m, niter
 integer, intent(in) :: nbf, nmo
 real(kind=8) :: rtmp, Aij, Bij, alpha, sin_4a, cos_a, sin_a, cc, ss, sc, &
  sin_2a, cos_2a, tot_change, sum_change
 real(kind=8), intent(inout) :: coeff(nbf,nmo)
 real(kind=8), external :: cmplx_v3_square, cmplx_v3_dot_real
 complex(kind=8) :: vtmp(3,4), vdiff(3)
 complex(kind=8), allocatable :: dipole(:,:,:)
 complex(kind=8), intent(inout) :: mo_dipole(3,nmo,nmo)

 write(6,'(/,A)') 'Perform serial 2*2 rotation...'
 allocate(dipole(3,nmo,2))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  sum_change = 0d0

  do m = 1, npair, 1
   i = ijmap(1,m); j = ijmap(2,m)
   vtmp(:,1) = mo_dipole(:,i,i)
   vtmp(:,2) = mo_dipole(:,i,j)
   vtmp(:,3) = mo_dipole(:,j,j)
   vdiff = vtmp(:,1) - vtmp(:,3)
   Aij = cmplx_v3_square(vtmp(:,2)) - 0.25d0*cmplx_v3_square(vdiff)
   Bij = cmplx_v3_dot_real(vdiff, vtmp(:,2))
   rtmp = HYPOT(Aij, Bij)
   sin_4a = Bij/rtmp
   rtmp = rtmp + Aij
   if(rtmp < upd_thres) cycle

   sum_change = sum_change + rtmp
   alpha = 0.25d0*DASIN(MAX(-1d0, MIN(sin_4a, 1d0)))
   if(Aij > 0d0) then
    alpha = QPI - alpha
   else if(Aij<0d0 .and. Bij<0d0) then
    alpha = HPI + alpha
   end if
   ! if alpha>PI/4, make an equivalent rotation so that these two orbital do
   ! not change significantly
   if(alpha > QPI) alpha = alpha - HPI
   cos_a = DCOS(alpha); sin_a = DSIN(alpha)

   ! update two orbitals
   call rotate_mo_ij(cos_a, sin_a, nbf, coeff(:,i), coeff(:,j))

   ! update corresponding dipole integrals, only indices in range to be updated
   cc = cos_a*cos_a; ss = sin_a*sin_a; sc = sin_a*cos_a
   cos_2a = 2d0*cc-1d0; sin_2a = 2d0*sc
   vtmp(:,4) = sin_2a*vtmp(:,2)
   dipole(:,i,1) = cc*vtmp(:,1) + ss*vtmp(:,3) + vtmp(:,4)
   dipole(:,j,2) = ss*vtmp(:,1) + cc*vtmp(:,3) - vtmp(:,4)
   dipole(:,j,1) = cos_2a*vtmp(:,2) - sc*vdiff
   dipole(:,i,2) = dipole(:,j,1)

   ! It seems that OpenMP makes this loop slower
   do k = 1, nmo, 1
    if(k==i .or. k==j) cycle
    dipole(:,k,1) = cos_a*mo_dipole(:,k,i) + sin_a*mo_dipole(:,k,j)
    dipole(:,k,2) = cos_a*mo_dipole(:,k,j) - sin_a*mo_dipole(:,k,i)
   end do ! for k

   mo_dipole(:,:,i) = dipole(:,:,1)
   mo_dipole(:,:,j) = dipole(:,:,2)
   mo_dipole(:,i,:) = dipole(:,:,1)
   mo_dipole(:,j,:) = dipole(:,:,2)
  end do ! for m

  tot_change = tot_change + sum_change
  niter = niter + 1
  write(6,'(A,I4,A,F16.8)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for while

 deallocate(dipole)
 write(6,'(A,F20.8)') 'tot_change=', tot_change
 call prt_loc_conv_remark(niter, max_niter, 'serial2by2_cmplx')
end subroutine serial2by2_cmplx

! perform serial 2-by-2 rotation on given MOs
subroutine serial2by2(nbf, nmo, ncomp, coeff, mo_dipole)
 use lo_info, only: npair, ijmap, max_niter, QPI, HPI, upd_thres, conv_thres
 implicit none
 integer :: i, j, k, m, niter
 integer, intent(in) :: nbf, nmo, ncomp
 real(kind=8) :: ddot, rtmp, tot_change, sum_change, Aij, Bij, alpha, sin_4a, &
  cos_a, sin_a, cc, ss, sc, sin_2a, cos_2a
 real(kind=8), intent(inout) :: coeff(nbf,nmo), mo_dipole(ncomp,nmo,nmo)
 ! for Boys, ncomp = 3
 ! for PM,   ncomp = natom, mo_dipole is actually the population matrix
 real(kind=8), allocatable :: dipole(:,:,:), vtmp(:,:), vdiff(:)
 ! dipole: for Boys, store updated dipole integrals matrix
 !         for PM,   store updated population matrix

 write(6,'(/,A)') 'Perform serial 2*2 rotation...'
 allocate(dipole(ncomp,nmo,2), vtmp(ncomp,4), vdiff(ncomp))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  sum_change = 0d0

  do m = 1, npair, 1
   i = ijmap(1,m); j = ijmap(2,m)
   vtmp(:,1) = mo_dipole(:,i,i)
   vtmp(:,2) = mo_dipole(:,i,j)
   vtmp(:,3) = mo_dipole(:,j,j)
   vdiff = vtmp(:,1) - vtmp(:,3)
   Aij = ddot(ncomp,vtmp(:,2),1,vtmp(:,2),1) - 0.25d0*ddot(ncomp,vdiff,1,vdiff,1)
   Bij = ddot(ncomp, vdiff, 1, vtmp(:,2), 1)
   rtmp = HYPOT(Aij, Bij)
   sin_4a = Bij/rtmp
   rtmp = rtmp + Aij
   if(rtmp < upd_thres) cycle

   sum_change = sum_change + rtmp
   alpha = 0.25d0*DASIN(MAX(-1d0, MIN(sin_4a, 1d0)))
   if(Aij > 0d0) then
    alpha = QPI - alpha
   else if(Aij<0d0 .and. Bij<0d0) then
    alpha = HPI + alpha
   end if
   ! if alpha>PI/4, make an equivalent rotation so that these two orbital do
   ! not change significantly
   if(alpha > QPI) alpha = alpha - HPI
   cos_a = DCOS(alpha); sin_a = DSIN(alpha)

   ! update two orbitals
   call rotate_mo_ij(cos_a, sin_a, nbf, coeff(:,i), coeff(:,j))

   ! update corresponding dipole integrals, only indices in range to be updated
   cc = cos_a*cos_a; ss = sin_a*sin_a; sc = sin_a*cos_a
   cos_2a = 2d0*cc-1d0; sin_2a = 2d0*sc
   vtmp(:,4) = sin_2a*vtmp(:,2)
   dipole(:,i,1) = cc*vtmp(:,1) + ss*vtmp(:,3) + vtmp(:,4)
   dipole(:,j,2) = ss*vtmp(:,1) + cc*vtmp(:,3) - vtmp(:,4)
   dipole(:,j,1) = cos_2a*vtmp(:,2) - sc*vdiff
   dipole(:,i,2) = dipole(:,j,1)

   ! It seems that OpenMP makes this loop slower
   do k = 1, nmo, 1
    if(k==i .or. k==j) cycle
    dipole(:,k,1) = cos_a*mo_dipole(:,k,i) + sin_a*mo_dipole(:,k,j)
    dipole(:,k,2) = cos_a*mo_dipole(:,k,j) - sin_a*mo_dipole(:,k,i)
   end do ! for k

   mo_dipole(:,:,i) = dipole(:,:,1)
   mo_dipole(:,:,j) = dipole(:,:,2)
   mo_dipole(:,i,:) = dipole(:,:,1)
   mo_dipole(:,j,:) = dipole(:,:,2)
  end do ! for m

  tot_change = tot_change + sum_change
  niter = niter + 1
  write(6,'(A,I4,A,F16.8)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for while

 deallocate(vdiff, vtmp, dipole)
 write(6,'(A,F20.8)') 'tot_change=', tot_change
 call prt_loc_conv_remark(niter, max_niter, 'serial2by2')
end subroutine serial2by2

subroutine serial22berry_kernel(nbf, nmo, mo, mo_zdip, change)
 use lo_info, only: QPI, HPI, upd_thres, diis, eff_npair, eff_ijmap, u
 implicit none
 integer :: i, j, m
 integer, intent(in) :: nbf, nmo
 real(kind=8) :: rtmp, Aij, Bij, alpha, sin_4a, cos_a, sin_a, cc, ss, sc, &
  sin_2a, cos_2a
 real(kind=8), intent(inout) :: mo(nbf,nmo)
 real(kind=8), intent(out) :: change
 real(kind=8), external :: cmplx_v3_square, cmplx_v3_dot_real
 complex(kind=8) :: vdiff(3), vtmp(3,4)
 complex(kind=8), intent(inout) :: mo_zdip(3,nmo,nmo)

 change = 0d0

 do m = 1, eff_npair, 1
  i = eff_ijmap(1,m); j = eff_ijmap(2,m)
  vtmp(:,1) = mo_zdip(:,i,i)
  vtmp(:,2) = mo_zdip(:,i,j)
  vtmp(:,3) = mo_zdip(:,j,j)
  vdiff = vtmp(:,1) - vtmp(:,3)
  Aij = cmplx_v3_square(vtmp(:,2)) - 0.25d0*cmplx_v3_square(vdiff)
  Bij = cmplx_v3_dot_real(vdiff, vtmp(:,2))
  rtmp = HYPOT(Aij, Bij)
  sin_4a = Bij/rtmp
  rtmp = rtmp + Aij
  if(rtmp < upd_thres) cycle

  change = change + rtmp
  alpha = 0.25d0*DASIN(MAX(-1d0, MIN(sin_4a, 1d0)))
  if(Aij > 0d0) then
   alpha = QPI - alpha
  else if(Aij<0d0 .and. Bij<0d0) then
   alpha = HPI + alpha
  end if
  if(alpha > QPI) alpha = alpha - HPI
  cos_a = DCOS(alpha); sin_a = DSIN(alpha)

  ! update two orbitals
  if(diis) call rotate_mo_ij(cos_a, sin_a, nmo, u(:,i), u(:,j))
  call rotate_mo_ij(cos_a, sin_a, nbf, mo(:,i), mo(:,j))

  ! update corresponding dipole integrals
  cc = cos_a*cos_a; ss = sin_a*sin_a; sc = sin_a*cos_a
  cos_2a = 2d0*cc-1d0; sin_2a = 2d0*sc
  vtmp(:,4) = sin_2a*vtmp(:,2)
  mo_zdip(:,i,i) = cc*vtmp(:,1) + ss*vtmp(:,3) + vtmp(:,4)
  mo_zdip(:,j,j) = ss*vtmp(:,1) + cc*vtmp(:,3) - vtmp(:,4)
  mo_zdip(:,i,j) = cos_2a*vtmp(:,2) - sc*vdiff
  call update_up_tri_ij_zdip(nmo, i, j, cos_a, sin_a, mo_zdip)
 end do ! for m

 deallocate(eff_ijmap)
end subroutine serial22berry_kernel

! perform serial Berry 2-by-2 Jacobi rotations on given MOs with distance considered
subroutine serial22berry(natom, nbf, nmo, bfirst, dis, mo, mo_zdip)
 use lo_info, only: npair, max_niter, dis_thres, conv_thres, diis, ndiis, &
  ijmap, u, k_old, cayley_k, k_diis, binfile, mo_cen, find_eff_ijmap
 implicit none
 integer :: i
 integer, intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
 real(kind=8) :: fberry0, fberry, sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: mo(nbf,nmo)
 real(kind=8), allocatable :: u_tmp(:,:), mo_dis(:)
 complex(kind=8), intent(inout) :: mo_zdip(3,nmo,nmo)
 logical :: solve_new_mo

 write(6,'(/,A)') 'Perform serial Berry 2*2 rotation...'
 write(6,'(A,L1,A,F7.2,A,ES9.2,A,I0)') 'DIIS=',diis,', dis_thres=',dis_thres, &
  ', conv_thres=', conv_thres, ', nmo=', nmo
 allocate(mo_cen(nmo), mo_dis(npair))
 tot_change = 0d0

 do i = 1, max_niter, 1
  solve_new_mo = .true.

  if(diis) then
   if((i>1 .and. i<ndiis+3) .or. (i>ndiis+2 .and. MOD(i,2)==1)) then
    call save_cayley_k_old(nmo, ndiis, k_old, binfile)
   end if
   if(i>ndiis .and. MOD(i,2)==0) then
    allocate(u_tmp(nmo,nmo))
    call get_mo_and_u_from_cayley_k_diis(nbf, nmo, ndiis, ijmap, binfile, &
                                     k_diis, u, mo, k_old, cayley_k, u_tmp)
    call get_fberry(nmo, mo_zdip, fberry0)
    call symmetrize_mo_zdip(nmo, mo_zdip) ! remember to symmetrize mo_zdip
    call update_mo_zdip_by_u(nmo, u_tmp, mo_zdip)
    deallocate(u_tmp)
    call get_fberry(nmo, mo_zdip, fberry)
    sum_change = fberry - fberry0
    solve_new_mo = .false.
   end if
  end if

  if(solve_new_mo) then
   if(diis) k_old = cayley_k
   call get_mo_center_by_scpa(natom, nbf, nmo, bfirst, mo, mo_cen)
   call atm_dis2mo_dis(natom, nmo, npair, dis, mo_cen, mo_dis)
   call find_eff_ijmap(nmo, mo_dis)
   call serial22berry_kernel(nbf, nmo, mo, mo_zdip, sum_change)
   if(diis) then
    call cayley_trans(nmo, u, cayley_k)
    call save_cayley_k_diff(nmo, ndiis, ijmap, binfile, k_old, cayley_k, k_diis)
   end if
  end if

  write(6,'(A,I4,A,F16.8)') 'niter=', i, ', sum_change=', sum_change
  tot_change = tot_change + sum_change
  if( (diis.and.solve_new_mo) .or. (.not.diis) ) then
   if(sum_change < conv_thres) exit
  end if
 end do ! for i

 deallocate(mo_cen, mo_dis)
 write(6,'(A,F20.8)') 'tot_change=', tot_change
 call prt_loc_conv_remark(i, max_niter, 'serial22berry')
end subroutine serial22berry

subroutine para22berry(natom, nbf, nmo, bfirst, dis, mo, mo_zdip)
 use lo_info, only: np, npair, nsweep, max_niter, dis_thres, conv_thres, diis, &
  ndiis, u, k_old, cayley_k, k_diis, binfile, ijmap, rrmap, mo_cen, rot_idx, &
  screen_jacobi_idx_by_dis, para22berry_kernel
 implicit none
 integer :: i
 integer, intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
 real(kind=8) :: fberry0, fberry, sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: mo(nbf,nmo)
 real(kind=8), allocatable :: u_tmp(:,:), mo_dis(:)
 complex(kind=8), intent(inout) :: mo_zdip(3,nmo,nmo)
 logical :: solve_new_mo

 if(nmo < 4) then
  write(6,'(/,A)') 'ERROR in subroutine para22berry: nmo<4. Too few orbitals. T&
                   &he subroutine'
  write(6,'(A)') 'serial22berry is supposed to be called in this case.'
  write(6,'(A,3I7)') 'natom, nbf, nmo=', natom, nbf, nmo
  stop
 end if

 write(6,'(/,A)') 'Perform parallel Berry 2*2 rotation...'
 write(6,'(A,L1,A,F7.2,A,ES9.2,A,I0)') 'DIIS=',diis,', dis_thres=',dis_thres, &
  ', conv_thres=', conv_thres, ', nmo=', nmo
 np = (nmo+1)/2
 nsweep = 4*np - 2
 allocate(rrmap(2,np,nsweep))
 call init_ring_jacobi_idx(nmo, np, rrmap)

 allocate(mo_cen(nmo), mo_dis(npair))
 tot_change = 0d0

 do i = 1, max_niter, 1
  solve_new_mo = .true.

  if(diis) then
   if((i>1 .and. i<ndiis+3) .or. (i>ndiis+2 .and. MOD(i,2)==1)) then
    call save_cayley_k_old(nmo, ndiis, k_old, binfile)
   end if
   if(i>ndiis .and. MOD(i,2)==0) then
    allocate(u_tmp(nmo,nmo))
    call get_mo_and_u_from_cayley_k_diis(nbf, nmo, ndiis, ijmap, binfile, &
                                     k_diis, u, mo, k_old, cayley_k, u_tmp)
    call get_fberry(nmo, mo_zdip, fberry0)
    call symmetrize_mo_zdip(nmo, mo_zdip) ! remember to symmetrize mo_zdip
    call update_mo_zdip_by_u(nmo, u_tmp, mo_zdip)
    deallocate(u_tmp)
    call get_fberry(nmo, mo_zdip, fberry)
    sum_change = fberry - fberry0
    solve_new_mo = .false.
   end if
  end if

  if(solve_new_mo) then
   if(diis) k_old = cayley_k
   call get_mo_center_by_scpa(natom, nbf, nmo, bfirst, mo, mo_cen)
   call atm_dis2mo_dis(natom, nmo, npair, dis, mo_cen, mo_dis)
   call screen_jacobi_idx_by_dis(nmo, mo_dis)
   call para22berry_kernel(nbf, nmo, mo, mo_zdip, sum_change)
   if(diis) then
    call cayley_trans(nmo, u, cayley_k)
    call save_cayley_k_diff(nmo, ndiis, ijmap, binfile, k_old, cayley_k, k_diis)
   end if
  end if

  write(6,'(A,I4,A,F16.8)') 'niter=', i, ', sum_change=', sum_change
  tot_change = tot_change + sum_change
  if( (diis.and.solve_new_mo) .or. (.not.diis) ) then
   if(sum_change < conv_thres) exit
  end if
 end do ! for i

 deallocate(mo_cen, mo_dis, rrmap, rot_idx)
 write(6,'(A,F20.8)') 'tot_change=', tot_change
 call prt_loc_conv_remark(i, max_niter, 'para22berry')
end subroutine para22berry

subroutine serial22boys_kernel(nbf, nmo, coeff, mo_dip, change)
 use lo_info, only: QPI, HPI, upd_thres, eff_npair, eff_ijmap
 implicit none
 integer :: i, j, m
 integer, intent(in) :: nbf, nmo
 real(kind=8) :: ddot, rtmp, Aij, Bij, alpha, sin_4a, cos_a, sin_a, cc, ss, sc,&
  sin_2a, cos_2a, vdiff(3), vtmp(3,4)
 real(kind=8), intent(inout) :: coeff(nbf,nmo), mo_dip(3,nmo,nmo)
 real(kind=8), intent(out) :: change

 change = 0d0

 do m = 1, eff_npair, 1
  i = eff_ijmap(1,m); j = eff_ijmap(2,m)
  vtmp(:,1) = mo_dip(:,i,i)
  vtmp(:,2) = mo_dip(:,i,j)
  vtmp(:,3) = mo_dip(:,j,j)
  vdiff = vtmp(:,1) - vtmp(:,3)
  Aij = ddot(3,vtmp(:,2),1,vtmp(:,2),1) - 0.25d0*ddot(3,vdiff,1,vdiff,1)
  Bij = ddot(3, vdiff, 1, vtmp(:,2), 1)
  rtmp = HYPOT(Aij, Bij)
  sin_4a = Bij/rtmp
  rtmp = rtmp + Aij
  if(rtmp < upd_thres) cycle

  change = change + rtmp
  alpha = 0.25d0*DASIN(MAX(-1d0, MIN(sin_4a, 1d0)))
  if(Aij > 0d0) then
   alpha = QPI - alpha
  else if(Aij<0d0 .and. Bij<0d0) then
   alpha = HPI + alpha
  end if
  if(alpha > QPI) alpha = alpha - HPI
  cos_a = DCOS(alpha); sin_a = DSIN(alpha)

  ! update two orbitals
  call rotate_mo_ij(cos_a, sin_a, nbf, coeff(:,i), coeff(:,j))

  ! update corresponding dipole integrals
  cc = cos_a*cos_a; ss = sin_a*sin_a; sc = sin_a*cos_a
  cos_2a = 2d0*cc-1d0; sin_2a = 2d0*sc
  vtmp(:,4) = sin_2a*vtmp(:,2)
  mo_dip(:,i,i) = cc*vtmp(:,1) + ss*vtmp(:,3) + vtmp(:,4)
  mo_dip(:,j,j) = ss*vtmp(:,1) + cc*vtmp(:,3) - vtmp(:,4)
  mo_dip(:,i,j) = cos_2a*vtmp(:,2) - sc*vdiff
  call update_up_tri_ij_dip(nmo, i, j, cos_a, sin_a, mo_dip)
 end do ! for m

 deallocate(eff_ijmap)
end subroutine serial22boys_kernel

! perform serial Boys 2-by-2 Jacobi rotations on given MOs with distance considered
subroutine serial22boys(natom, nbf, nmo, bfirst, dis, mo, mo_dip)
 use lo_info, only: npair, max_niter, dis_thres, conv_thres, mo_cen, &
  find_eff_ijmap
 implicit none
 integer :: i
 integer, intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: mo(nbf,nmo), mo_dip(3,nmo,nmo)
 real(kind=8), allocatable :: mo_dis(:)

 write(6,'(/,A)') 'Perform serial Boys 2*2 rotation...'
 write(6,'(A,F8.2,A,ES9.2,A,I0)') 'dis_thres=', dis_thres, ', conv_thres=', &
                                   conv_thres, ', nmo=', nmo
 allocate(mo_cen(nmo), mo_dis(npair))
 tot_change = 0d0

 do i = 1, max_niter, 1
  call get_mo_center_by_scpa(natom, nbf, nmo, bfirst, mo, mo_cen)
  call atm_dis2mo_dis(natom, nmo, npair, dis, mo_cen, mo_dis)
  call find_eff_ijmap(nmo, mo_dis)
  call serial22boys_kernel(nbf, nmo, mo, mo_dip, sum_change)
  tot_change = tot_change + sum_change
  write(6,'(A,I4,A,F16.8)') 'niter=', i, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for i

 deallocate(mo_cen, mo_dis)
 write(6,'(A,F20.8)') 'tot_change=', tot_change
 call prt_loc_conv_remark(i, max_niter, 'serial22boys')
end subroutine serial22boys

! perform parallel Foster-Boys 2-by-2 Jacobi rotation on given MOs
subroutine para22boys(natom, nbf, nmo, bfirst, dis, coeff, mo_dip)
 use lo_info, only: np, npair, nsweep, max_niter, dis_thres, conv_thres, &
  mo_cen, rrmap, rot_idx, screen_jacobi_idx_by_dis, para22boys_kernel
 implicit none
 integer :: i
 integer, intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: coeff(nbf,nmo), mo_dip(3,nmo,nmo)
 real(kind=8), allocatable :: mo_dis(:)

 if(nmo < 4) then
  write(6,'(/,A)') 'ERROR in subroutine para22boys: nmo<4. Too few orbitals. Th&
                   &e subroutine'
  write(6,'(A)') 'serial22boys is supposed to be called in this case.'
  write(6,'(A,3I7)') 'natom, nbf, nmo=', natom, nbf, nmo
  stop
 end if

 write(6,'(/,A)') 'Perform parallel Boys 2*2 rotation...'
 write(6,'(A,F8.2,A,ES9.2,A,I0)') 'dis_thres=', dis_thres, ', conv_thres=', &
                                   conv_thres, ', nmo=', nmo
 np = (nmo+1)/2
 nsweep = 4*np - 2
 allocate(rrmap(2,np,nsweep))
 call init_ring_jacobi_idx(nmo, np, rrmap)

 allocate(mo_cen(nmo), mo_dis(npair))
 tot_change = 0d0

 do i = 1, max_niter, 1
  call get_mo_center_by_scpa(natom, nbf, nmo, bfirst, coeff, mo_cen)
  call atm_dis2mo_dis(natom, nmo, npair, dis, mo_cen, mo_dis)
  call screen_jacobi_idx_by_dis(nmo, mo_dis)
  call para22boys_kernel(nbf, nmo, coeff, mo_dip, sum_change)
  tot_change = tot_change + sum_change
  write(6,'(A,I4,A,F16.8)') 'niter=', i, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for i

 deallocate(mo_cen, mo_dis, rrmap, rot_idx)
 write(6,'(A,F20.8)') 'tot_change=', tot_change
 call prt_loc_conv_remark(i, max_niter, 'para22boys')
end subroutine para22boys

subroutine serial22pm_kernel(nbf, nmo, natom, coeff, gross, change)
 use lo_info, only: QPI, HPI, upd_thres, eff_npair, eff_ijmap
 implicit none
 integer :: i, j, m
 integer, intent(in) :: nbf, nmo, natom
 real(kind=8) :: ddot, rtmp, Aij, Bij, alpha, sin_4a, cos_a, sin_a
 real(kind=8), intent(inout) :: coeff(nbf,nmo), gross(natom,nmo,nmo)
 real(kind=8), intent(out) :: change
 real(kind=8), allocatable :: vdiff(:), vtmp(:,:), tmp_mo(:)

 change = 0d0
 allocate(vdiff(natom), vtmp(natom,3), tmp_mo(nbf))

 do m = 1, eff_npair, 1
  i = eff_ijmap(1,m); j = eff_ijmap(2,m)
  call dcopy(natom, gross(:,i,i), 1, vtmp(:,1), 1)
  call dcopy(natom, gross(:,i,j), 1, vtmp(:,2), 1)
  call dcopy(natom, gross(:,j,j), 1, vtmp(:,3), 1)
  call dcopy(natom, vtmp(:,1)-vtmp(:,3), 1, vdiff, 1)
  Aij = ddot(natom,vtmp(:,2),1,vtmp(:,2),1) - 0.25d0*ddot(natom,vdiff,1,vdiff,1)
  Bij = ddot(natom, vdiff, 1, vtmp(:,2), 1)
  rtmp = HYPOT(Aij, Bij)
  sin_4a = Bij/rtmp
  rtmp = rtmp + Aij
  if(rtmp < upd_thres) cycle

  change = change + rtmp
  alpha = 0.25d0*DASIN(MAX(-1d0, MIN(sin_4a, 1d0)))
  if(Aij > 0d0) then
   alpha = QPI - alpha
  else if(Aij<0d0 .and. Bij<0d0) then
   alpha = HPI + alpha
  end if
  if(alpha > QPI) alpha = alpha - HPI
  cos_a = DCOS(alpha); sin_a = DSIN(alpha)

  call rotate_mo_ij(cos_a, sin_a, nbf, coeff(:,i), coeff(:,j))
  call update_gross_ii_jj_ji(cos_a, sin_a, natom, vtmp, vdiff, gross(:,i,i), &
                             gross(:,j,j), gross(:,i,j))
  call update_up_tri_ij_gross(natom, nmo, i, j, cos_a, sin_a, gross)
 end do ! for m

 deallocate(vdiff, vtmp, tmp_mo, eff_ijmap)
end subroutine serial22pm_kernel

! perform serial Pipek-Mezey 2-by-2 Jacobi rotation on given MOs with distance considered
subroutine serial22pm(natom, nbf, nmo, dis, coeff, gross)
 use lo_info, only: npair, max_niter, dis_thres, conv_thres, mo_cen, &
  find_eff_ijmap
 implicit none
 integer :: i
 integer, intent(in) :: natom, nbf, nmo
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: coeff(nbf,nmo), gross(natom,nmo,nmo)
 real(kind=8), allocatable :: mo_dis(:)

 write(6,'(/,A)') 'Perform serial PM 2*2 rotation...'
 write(6,'(A,F8.2,A,ES9.2,A,I0)') 'dis_thres=', dis_thres, ', conv_thres=', &
                                   conv_thres, ', nmo=', nmo
 allocate(mo_cen(nmo), mo_dis(npair))
 tot_change = 0d0

 do i = 1, max_niter, 1
  call get_mo_center_from_gross(natom, nmo, gross, mo_cen)
  call atm_dis2mo_dis(natom, nmo, npair, dis, mo_cen, mo_dis)
  call find_eff_ijmap(nmo, mo_dis)
  call serial22pm_kernel(nbf, nmo, natom, coeff, gross, sum_change)
  tot_change = tot_change + sum_change
  write(6,'(A,I4,A,F16.8)') 'niter=', i, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for i

 deallocate(mo_cen, mo_dis)
 write(6,'(A,F20.8)') 'tot_change=', tot_change
 call prt_loc_conv_remark(i, max_niter, 'serial22pm')
end subroutine serial22pm

! perform parallel Pipek-Mezey 2-by-2 Jacobi rotation on given MOs
subroutine para22pm(natom, nbf, nmo, dis, coeff, gross)
 use lo_info, only: np, npair, nsweep, max_niter, dis_thres, conv_thres, &
  mo_cen, rrmap, rot_idx, screen_jacobi_idx_by_dis, para22pm_kernel
 implicit none
 integer :: i
 integer, intent(in) :: natom, nbf, nmo
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: coeff(nbf,nmo), gross(natom,nmo,nmo)
 real(kind=8), allocatable :: mo_dis(:)

 if(nmo < 4) then
  write(6,'(/,A)') 'ERROR in subroutine para22pm: nmo<4. Too few orbitals. The &
                   &subroutine'
  write(6,'(A)') 'serial22pm is supposed to be called in this case.'
  write(6,'(A,3I7)') 'natom, nbf, nmo=', natom, nbf, nmo
  stop
 end if

 write(6,'(/,A)') 'Perform parallel PM 2*2 rotation...'
 write(6,'(A,F8.2,A,ES9.2,A,I0)') 'dis_thres=', dis_thres, ', conv_thres=', &
                                   conv_thres, ', nmo=', nmo
 np = (nmo+1)/2
 nsweep = 4*np - 2
 allocate(rrmap(2,np,nsweep))
 call init_ring_jacobi_idx(nmo, np, rrmap)

 allocate(mo_cen(nmo), mo_dis(npair))
 tot_change = 0d0

 do i = 1, max_niter, 1
  call get_mo_center_from_gross(natom, nmo, gross, mo_cen)
  call atm_dis2mo_dis(natom, nmo, npair, dis, mo_cen, mo_dis)
  call screen_jacobi_idx_by_dis(nmo, mo_dis)
  call para22pm_kernel(nbf, nmo, natom, coeff, gross, sum_change)
  tot_change = tot_change + sum_change
  write(6,'(A,I4,A,F16.8)') 'niter=', i, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for i

 deallocate(mo_cen, mo_dis, rrmap, rot_idx)
 write(6,'(A,F20.8)') 'tot_change=', tot_change
 call prt_loc_conv_remark(i, max_niter, 'para22pm')
end subroutine para22pm

!subroutine boys_polar(nmo, mo_dip, u)
! implicit none
! integer :: i, j, n
! integer, intent(in) :: nmo
! integer, parameter :: maxit = 500
! real(kind=8) :: ddot
! real(kind=8), intent(in) :: mo_dip(3,nmo,nmo)
! real(kind=8), intent(out) :: u(nmo,nmo)
! real(kind=8), allocatable :: dc(:,:,:), u0(:,:), x(:,:), a(:,:), a_u(:,:), &
!  a_vt(:,:), a_s(:)
!
! call init_identity_mat(nmo, u)
! allocate(u0(nmo,nmo))
!
! do n = 1, maxit, 1
!  u0 = u
!  allocate(dc(nmo,nmo,3), source=0d0)
!  call dsymm('N','N',nmo,nmo,nmo,1d0,mo_dip(1,:,:),nmo,u,nmo,0d0,dc(:,:,1),nmo)
!  call dsymm('N','N',nmo,nmo,nmo,1d0,mo_dip(2,:,:),nmo,u,nmo,0d0,dc(:,:,2),nmo)
!  call dsymm('N','N',nmo,nmo,nmo,1d0,mo_dip(3,:,:),nmo,u,nmo,0d0,dc(:,:,3),nmo)
!  allocate(x(nmo,3))
!  !$omp parallel do schedule(dynamic) default(shared) private(i,j) collapse(2)
!  do i = 1, 3
!   do j = 1, nmo, 1
!    x(j,i) = ddot(nmo, u(:,j), 1, dc(:,j,i), 1)
!   end do ! for j
!  end do ! for i
!  !$omp end parallel do
!  allocate(a(nmo,nmo))
!  !$omp parallel do schedule(dynamic) default(shared) private(i,j) collapse(2)
!  do i = 1, nmo, 1
!   do j = 1, nmo, 1
!    a(j,i) = dc(j,i,1)*x(i,1) + dc(j,i,2)*x(i,2) + dc(j,i,3)*x(i,3)
!   end do ! for j
!  end do ! for i
!  !$omp end parallel do
!  deallocate(dc, x)
!  allocate(a_u(nmo,nmo), a_vt(nmo,nmo), a_s(nmo))
!  call do_svd(nmo, nmo, a, a_u, a_vt, a_s)
!  deallocate(a, a_s)
!  u = 0d0
!  call dgemm('N', 'N', nmo, nmo, nmo, 1d0, a_u, nmo, a_vt, nmo, 0d0, u, nmo)
!  deallocate(a_u, a_vt)
! end do ! for n
!
! deallocate(u0)
! write(6,'(A,I0)') 'niter=', n
!end subroutine boys_polar

