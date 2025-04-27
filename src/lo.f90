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
 write(fid,'(A)') 'from pyscf import gto'
 write(fid,'(A)') 'from shutil import copyfile'
 write(fid,'(A)') 'from mokit.lib.gaussian import BOHR2ANG, load_mol_from_fch'
 write(fid,'(A)') 'from mokit.lib.rwwfn import read_nbf_and_nif_from_fch, \'
 write(fid,'(A)') ' read_eigenvalues_from_fch'
 write(fid,'(A)') 'from mokit.lib.fch2py import fch2py'
 write(fid,'(A)') 'from mokit.lib.py2fch import py2fch'
 if(pm_loc) then
  write(fid,'(A)') 'from mokit.lib.lo import gen_loc_ini_guess, pm'
 else
  write(fid,'(A)') 'from mokit.lib.gaussian import ao_dipole_int'
  write(fid,'(A)') 'from mokit.lib.lo import gen_loc_ini_guess, boys'
 end if
 write(fid,'(A)') 'import numpy as np'
 write(fid,'(A)') 'import os'

 write(fid,'(/,A)') "fchname = '"//TRIM(fchname)//"'"
 write(fid,'(A)') "lmofch = '"//TRIM(lmofch)//"'"
 write(fid,'(A)') 'mol = load_mol_from_fch(fchname)'
 write(fid,'(A)') 'nbf, nif = read_nbf_and_nif_from_fch(fchname)'
 write(fid,'(A)') "mo_coeff = fch2py(fchname, nbf, nif, 'a')"
 write(fid,'(2(A,I0),A)') 'idx = range(',i1-1,',',i2,')'
 write(fid,'(A)') 'nmo = len(idx)'
 write(fid,'(A)') 'natom = mol.natm'
 write(fid,'(A)') 'dis = BOHR2ANG*gto.inter_distance(mol)'
 write(fid,'(A)') 'bfirst = np.ones(natom+1, dtype=np.int32)'
 write(fid,'(A)') 'bfirst[1:] = mol.aoslice_by_atom()[:,3] + 1'
 write(fid,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
 write(fid,'(A)') 'chosen = np.ones(natom, dtype=bool)'
 write(fid,'(A)') 'nmo1, lmo_ini = gen_loc_ini_guess(natom,nbf,nmo,chosen,bfirs&
                  &t,S,mo_coeff[:,idx])'
 if(pm_loc) then
  write(fid,'(A)') "loc_orb = pm(natom,nbf,nmo,bfirst,dis,lmo_ini,S,'mulliken',&
                   &17.0,1e-5)"
 else
  write(fid,'(A)') 'center, ao_dip = ao_dipole_int(mol)'
  write(fid,'(A)') 'ao_dip = ao_dip*BOHR2ANG'
  write(fid,'(A)') 'loc_orb = boys(natom,nbf,nmo,bfirst,dis,lmo_ini,ao_dip,17.0&
                   &,1e-5)'
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
  write(6,'(A)') 'ERROR in subroutine gen_cf_orb: ndb<0 or nopen<0 found.'
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
 call gen_ortho_mo(nbf, ao_ovlp, k, orth_mo)
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
! According to jxzou's tests, the locality of obtained LMOs are usually
!  atomic > Cho2 > Cho, where Cho2 means Cholesky LMOs based on orthogonal AO
!  basis.
subroutine resemble_ao(nbf, nif, ao_ovlp, mo)
 implicit none
 integer :: i, k
 integer, intent(in) :: nbf, nif
 integer, allocatable :: idx(:)
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
 real(kind=8), intent(inout) :: mo(nbf,nif)
 real(kind=8), allocatable :: orth_mo(:,:), u0(:,:), u1(:,:), u(:,:), vt(:,:), &
  uvt(:,:), norm(:), s(:), d(:)

 allocate(orth_mo(nbf,nbf))
 call gen_ortho_mo(nbf, ao_ovlp, k, orth_mo)
 if(k < nbf) then
  write(6,'(/,A)') 'ERROR in subroutine resemble_ao: linear dependency detecte&
                   &d in AO overlap'
  write(6,'(A)') 'used. Not tested so far. You can discard this geometry and u&
                 &se another one,'
  write(6,'(A)') 'or simply contact MOKIT developers.'
  stop
 end if

 ! calculate the unitary matrix U in equation `mo = orth_mo*U_0`
 ! U_0 = (orth_mo^T)S(mo)
 allocate(u0(k,nif))
 call calc_CTSCp2(nbf, k, nif, orth_mo(:,1:k), ao_ovlp, mo, u0)
 deallocate(orth_mo)

 ! find the AOs which has largest overlap with input MOs
 allocate(norm(k))
 forall(i = 1:k) norm(i) = DOT_PRODUCT(u0(i,:), u0(i,:))
 allocate(idx(k))
 call sort_dp_array(k, norm, .false., idx)
 deallocate(norm)
 allocate(u1(nif,nif))
 forall(i = 1:nif) u1(i,:) = u0(idx(i),:)
 deallocate(idx, u0)

 ! compute the rotation matrix of original MOs
 allocate(u(nif,nif), vt(nif,nif), s(nif))
 call do_svd(nif, nif, u1, u, vt, s)
 deallocate(u1)
 allocate(d(nif))
 call calc_usut_diag_elem(nif, s, u, d)
 deallocate(s)
 allocate(idx(nif))
 call sort_dp_array(nif, d, .false., idx)
 deallocate(d)

 ! U(V^T)
 allocate(uvt(nif,nif), source=0d0)
 call dgemm('N', 'N', nif, nif, nif, 1d0, u, nif, vt, nif, 0d0, uvt, nif)
 deallocate(u, vt)

 ! MO*(V(U^T))
 allocate(u0(nbf,nif), source=0d0)
 call dgemm('N','T', nbf, nif, nif, 1d0, mo, nbf, uvt, nif, 0d0, u0, nbf)
 deallocate(uvt)

 ! MOs in u0 are not necessarily in overlap descending order, sorting using
 ! idx is needed.
 forall(i = 1:nif) mo(:,i) = u0(:,idx(i))
 deallocate(idx, u0)
end subroutine resemble_ao

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
  call resemble_ao(nbf, nif, ao_ovlp, new_mo)
  write(6,'(A)') 'Done construction.'
  return
 end if

 write(6,'(A)') 'Not all atoms are chosen. Only orbitals centered on specified &
                &atoms to be solved...'
 allocate(idx(nbf), source=0)
 nbf1 = 0

 do i = 1, natom, 1
  if(chosen(i)) then
   k = bfirst(i+1) - bfirst(i)
   forall(j=1:k) idx(nbf1+j) = bfirst(i)+j-1
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
 call gen_ortho_mo(nbf1, ao_ovlp1, k, orth_mo)

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
 call resemble_ao(nbf1, nif1, ao_ovlp1, proj_mo(:,1:nif1))
 deallocate(ao_ovlp1)

 ! Make new_mo(:,1:nif1) resembles proj_mo(:,1:nif1).
 ! new_mo(:,nif1+1:) are remaining orbitals with non-zero MO coefficients, but
 ! will probably be discarded in subsequent orbital localization.
 call orb_resemble_ref1(nbf, nif, mo, nbf1, nif1, proj_mo, ao_ovlp2, new_mo)
 deallocate(proj_mo, ao_ovlp2)
 write(6,'(A)') 'Done construction.'
end subroutine gen_loc_ini_guess

! Combining AO-resemblance with Cholesky LMO techniques. The Cholesky LMO
! technique does not always make target orbitals spatially localized. This
! subroutine remains to be improved.
subroutine combined_ao_cholesky(nbf, nif, ao_ovlp, mo)
 implicit none
 integer :: i, j, k, nlmo
 integer, intent(in) :: nbf, nif
 integer, allocatable :: idx(:)
 real(kind=8), parameter :: ovlp_thres = 0.56d0
 ! threshold for selecting MOs transformed to Cholesky LMOs
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf)
 real(kind=8), intent(inout) :: mo(nbf,nif)
 real(kind=8), allocatable :: orth_mo(:,:), u0(:,:), u1(:,:), u(:,:), vt(:,:),&
  uvt(:,:), norm(:), s(:), d(:)

 allocate(orth_mo(nbf,nbf))
 call gen_ortho_mo(nbf, ao_ovlp, k, orth_mo)

 if(k < nbf) then
  write(6,'(/,A)') 'Warning from subroutine combined_ao_cholesky: linear depend&
                   &ency detected in'
  write(6,'(A)') 'given AO overlap. Initial LMOs not constructed.'
  return
 end if

 ! calculate the unitary matrix U in equation `mo = orth_mo*U_0`
 ! U_0 = (orth_mo^T)S(mo)
 allocate(u0(k,nif))
 call calc_CTSCp2(nbf, k, nif, orth_mo(:,1:k), ao_ovlp, mo, u0)

 ! find the AOs which has largest overlap with input MOs
 allocate(norm(k))
 forall(i = 1:k) norm(i) = DOT_PRODUCT(u0(i,:), u0(i,:))
 allocate(idx(k))
 call sort_dp_array(k, norm, .false., idx)
 deallocate(norm)
 allocate(u1(nif,nif))
 forall(i = 1:nif) u1(i,:) = u0(idx(i),:)
 deallocate(idx)

 ! compute the rotation matrix of original MOs
 allocate(u(nif,nif), vt(nif,nif), s(nif))
 call do_svd(nif, nif, u1, u, vt, s)
 deallocate(u1)
 allocate(d(nif))
 call calc_usut_diag_elem(nif, s, u, d)
 deallocate(s)
 allocate(idx(nif))
 call sort_dp_array(nif, d, .false., idx)
 nlmo = COUNT(d > ovlp_thres)
 deallocate(d)
 write(6,'(A,I0)') 'nlmo=', nlmo
 ! If nlmo=nif, all orbitals match well with AOs.
 ! If 0< nlmo <nif, some orbitals match well with AOs, for remaining orbitals
 ! we generate their Cholesky LMOs.
 ! If nlmo=0, no orbital match well with AOs, generate Cholesky LMOs for all
 ! orbitals.

 if(nlmo == 0) then
  deallocate(u, vt)
  call cholesky_mo(k, nif, u0)
  mo = 0d0
  call dgemm('N','N', nbf,nif,k, 1d0,orth_mo(:,1:k),nbf, u0,k, 0d0,mo,nbf)
  deallocate(u0)
 else
  ! U(V^T)
  deallocate(u0)
  allocate(uvt(nif,nif), source=0d0)
  call dgemm('N', 'N', nif, nif, nif, 1d0, u, nif, vt, nif, 0d0, uvt, nif)
  deallocate(u, vt)
  ! MO*(V(U^T))
  allocate(u1(nbf,nif), source=0d0)
  call dgemm('N','T', nbf, nif, nif, 1d0, mo, nbf, uvt, nif, 0d0, u1, nbf)
  deallocate(uvt)
  ! MOs in u1 are not necessarily in overlap descending order, sorting using
  ! idx is needed.
  forall(i = 1:nif) mo(:,i) = u1(:,idx(i))
  deallocate(idx, u1)
  if(nlmo < nif-1) then ! generate Cholesky LMOs for remaining orbitals
   j = nif - nlmo
   allocate(u0(k,j))
   allocate(u1(nbf,j), source=mo(:,nlmo+1:nif))
   call calc_CTSCp2(nbf, k, j, orth_mo(:,1:k), ao_ovlp, u1, u0)
   call cholesky_mo(k, j, u0)
   u1 = 0d0
   call dgemm('N','N', nbf,j,k, 1d0,orth_mo(:,1:k),nbf, u0,k, 0d0, u1, nbf)
   mo(:,nlmo+1:nif) = u1
   deallocate(u0, u1)
  end if
 end if

 deallocate(orth_mo)
end subroutine combined_ao_cholesky

! Calculate/get the integer array bfirst, which stores the AO index range of
! each atom.
!subroutine get_bfirst_from_shl(cart, nshl, natom, shl2atm, ang, ibas, nbf, bfirst)
! implicit none
! integer :: i, j
! integer, intent(in) :: nshl, natom
!!f2py intent(in) :: nshl, natom
! integer, intent(in) :: shl2atm(nshl), ang(nshl), ibas(nshl)
!!f2py intent(in) :: shl2atm, ang, ibas
!!f2py depend(nshl) :: shl2atm, ang, ibas
! integer, intent(out) :: nbf, bfirst(natom+1)
!!f2py intent(out) :: nbf, bfirst
!!f2py depend(natom) :: bfirst
! integer, allocatable :: ang1(:)
! logical, intent(in) :: cart
!!f2py intent(in) :: cart
!
! if(ANY(ang < 0)) then
!  write(6,'(/,A)') 'ERROR in subroutine get_bfirst_from_shl: there exists some &
!                   &ang(i)<0.'
!  write(6,'(A,L1,A,I0)') 'cart=', cart, ', nshl=', nshl
!  write(6,'(A)') 'ang='
!  write(6,'(20I4)') ang
!  stop
! end if
!
! allocate(ang1(nshl))
! if(cart) then
!  forall(i = 1:nshl) ang1(i) = (ang(i)+1)*(ang(i)+2)/2
! else
!  forall(i = 1:nshl) ang1(i) = 2*ang(i) + 1
! end if
!
! nbf = DOT_PRODUCT(ang1, ibas)
! ! return nbf, so that one is able to check whether the number of basis functions
! ! is equal to nbf here.
! write(6,'(A,I0)') 'natom=', natom
! bfirst = 0; bfirst(1) = 1
!
! do i = 1, nshl, 1
!  j = shl2atm(i) + 2
!  bfirst(j) = bfirst(j) + ang1(i)*ibas(i)
! end do ! for i
! deallocate(ang1)
!
! do i = 2, natom+1, 1
!  bfirst(i) = bfirst(i) + bfirst(i-1)
! end do ! for i
!end subroutine get_bfirst_from_shl

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
  vtmp(:,2) = gross(:,j,i)
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
subroutine berry(natom, nbf, nmo, bfirst, dis, mo, ao_zdip, dis_tol, conv_tol, &
                 new_mo)
 use lo_info, only: dis_thres, conv_thres, npair, ijmap
 implicit none
 integer, intent(in) :: natom, nbf, nmo
!f2py intent(in) :: natom, nbf, nmo
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

 write(6,'(/,A)') 'Berry orbital localization begins:'
 new_mo = mo; dis_thres = dis_tol; conv_thres = conv_tol

 if(nmo == 1) then
  write(6,'(A)') 'Warning in subroutine berry: only 1 orbital. No rotation.'
  return
 end if

 write(6,'(/,A)') 'Transform AO complex dipole integrals to MO ones...'
 allocate(mo_zdip(3,nmo,nmo))
 call ao2mo_zdip(nbf, nmo, new_mo, ao_zdip, mo_zdip)
 write(6,'(A)') 'Done update.'

 npair = nmo*(nmo-1)/2
 allocate(ijmap(2,npair))
 call get_triu_idx1(nmo, ijmap)

 !if(nmo < 10) then
  !call serial2by2_cmplx(nbf, nmo, new_mo, mo_zdip)
 !else if(nmo < 500) then
  call serial22berry(natom, nbf, nmo, bfirst, dis, new_mo, mo_zdip)
 !else
 ! call para22berry(natom, nbf, nmo, bfirst, dis, new_mo, mo_zdip)
 !end if

 ! job accomplished, no need to update the lower triangle part of mo_zdip
 deallocate(mo_zdip, ijmap)
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
 real(kind=8), intent(in) :: dis(natom,natom), mo(nbf,nmo), dis_tol, conv_tol
!f2py intent(in) :: dis, mo, dis_tol, conv_tol
!f2py depend(natom) :: dis
!f2py depend(nbf,nmo) :: mo
 real(kind=8), intent(out) :: new_mo(nbf,nmo)
!f2py intent(out) :: new_mo
!f2py depend(nbf,nmo) :: new_mo
 real(kind=8), intent(in) :: ao_dip(3,nbf,nbf)
!f2py intent(in) :: ao_dip
!f2py depend(nbf) :: ao_dip
 real(kind=8), allocatable :: mo_dip(:,:,:)

 write(6,'(/,A)') 'Boys orbital localization begins:'
 new_mo = mo; dis_thres = dis_tol; conv_thres = conv_tol

 if(nmo == 1) then
  write(6,'(A)') 'Warning in subroutine boys: only 1 orbital. No rotation.'
  return
 end if

 npair = nmo*(nmo-1)/2
 allocate(ijmap(2,npair))
 call get_triu_idx1(nmo, ijmap)

 write(6,'(/,A)') 'Transform AO dipole integrals to MO ones...'
 allocate(mo_dip(3,nmo,nmo))
 call ao2mo_dip(nbf, nmo, ao_dip, new_mo, mo_dip)
 write(6,'(A)') 'Done update.'

 if(nmo < 10) then
  call serial2by2(nbf, nmo, 3, new_mo, mo_dip)
 else if(nmo < 500) then
  call serial22boys(natom, nbf, nmo, bfirst, dis, new_mo, mo_dip)
 else
  call para22boys(natom, nbf, nmo, bfirst, dis, new_mo, mo_dip)
 end if

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
 character(len=*), intent(in) :: popm ! 'mulliken'/'lowdin'
!f2py intent(in) :: popm

 write(6,'(/,A)') 'PM orbital localization begins: using '//TRIM(popm)//&
                  ' population'
 new_mo = mo; dis_thres = dis_tol; conv_thres = conv_tol

 if(nmo == 1) then
  write(6,'(A)') 'Warning in subroutine pm: only 1 orbital. No rotation.'
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
 real(kind=8) :: rtmp, Aij, Bij, alpha, sin_4a, cos_a, sin_a, cc, ss, sin_2a, &
  cos_2a, tot_change, sum_change
 real(kind=8), intent(inout) :: coeff(nbf,nmo)
 real(kind=8), allocatable :: tmp_mo(:)
 real(kind=8), external :: cmplx_v3_square, cmplx_v3_dot_real
 complex(kind=8) :: vtmp(3,3), vdiff(3)
 complex(kind=8), allocatable :: dipole(:,:,:)
 complex(kind=8), intent(inout) :: mo_dipole(3,nmo,nmo)

 write(6,'(/,A)') 'Perform serial 2*2 rotation...'
 allocate(tmp_mo(nbf), dipole(3,nmo,2))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  sum_change = 0d0

  do m = 1, npair, 1
   i = ijmap(1,m); j = ijmap(2,m)
   vtmp(:,1) = mo_dipole(:,i,i)
   vtmp(:,2) = mo_dipole(:,j,i)
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
!dir$ ivdep
   tmp_mo = coeff(:,i)
!dir$ ivdep
   coeff(:,i) = cos_a*tmp_mo + sin_a*coeff(:,j)
!dir$ ivdep
   coeff(:,j) = cos_a*coeff(:,j) - sin_a*tmp_mo

   ! update corresponding dipole integrals, only indices in range to be updated
   cc = cos_a*cos_a
   ss = sin_a*sin_a
   sin_2a = 2d0*sin_a*cos_a
   cos_2a = cc - ss
   dipole(:,i,1) = cc*vtmp(:,1) + ss*vtmp(:,3) + sin_2a*vtmp(:,2)
   dipole(:,j,2) = ss*vtmp(:,1) + cc*vtmp(:,3) - sin_2a*vtmp(:,2)
   dipole(:,j,1) = cos_2a*vtmp(:,2) - 0.5d0*sin_2a*vdiff
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

 deallocate(tmp_mo, dipole)
 write(6,'(A,F20.8)') 'tot_change=', tot_change

 if(niter <= max_niter) then
  write(6,'(A)') 'Orbital localization converged successfully.'
 else
  write(6,'(/,A)') 'ERROR in subroutine serial2by2_cmplx: max_niter exceeded.'
  write(6,'(A,I0)') 'max_niter=', max_niter
  stop
 end if
end subroutine serial2by2_cmplx

! perform serial 2-by-2 rotation on given MOs
subroutine serial2by2(nbf, nmo, ncomp, coeff, mo_dipole)
 use lo_info, only: npair, ijmap, max_niter, QPI, HPI, upd_thres, conv_thres
 implicit none
 integer :: i, j, k, m, niter
 integer, intent(in) :: nbf, nmo, ncomp
 real(kind=8) :: ddot, rtmp, tot_change, sum_change
 real(kind=8) :: Aij, Bij, alpha, sin_4a, cos_a, sin_a
 real(kind=8) :: cc, ss, sin_2a, cos_2a
 real(kind=8), intent(inout) :: coeff(nbf,nmo), mo_dipole(ncomp,nmo,nmo)
 ! for Boys, ncomp = 3
 ! for PM,   ncomp = natom, mo_dipole is actually the population matrix
 real(kind=8), allocatable :: dipole(:,:,:), tmp_mo(:), vtmp(:,:), vdiff(:)
 ! dipole: for Boys, store updated dipole integrals matrix
 !         for PM,   store updated population matrix
 ! tmp_mo: store one MO

 write(6,'(/,A)') 'Perform serial 2*2 rotation...'
 allocate(dipole(ncomp,nmo,2), tmp_mo(nbf), vtmp(ncomp,3), vdiff(ncomp))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  sum_change = 0d0

  do m = 1, npair, 1
   i = ijmap(1,m); j = ijmap(2,m)
   vtmp(:,1) = mo_dipole(:,i,i)
   vtmp(:,2) = mo_dipole(:,j,i)
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
!dir$ ivdep
   tmp_mo = coeff(:,i)
!dir$ ivdep
   coeff(:,i) = cos_a*tmp_mo + sin_a*coeff(:,j)
!dir$ ivdep
   coeff(:,j) = cos_a*coeff(:,j) - sin_a*tmp_mo

   ! update corresponding dipole integrals, only indices in range to be updated
   cc = cos_a*cos_a
   ss = sin_a*sin_a
   sin_2a = 2d0*sin_a*cos_a
   cos_2a = cc - ss
   dipole(:,i,1) = cc*vtmp(:,1) + ss*vtmp(:,3) + sin_2a*vtmp(:,2)
   dipole(:,j,2) = ss*vtmp(:,1) + cc*vtmp(:,3) - sin_2a*vtmp(:,2)
   dipole(:,j,1) = cos_2a*vtmp(:,2) - 0.5d0*sin_2a*vdiff
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

 deallocate(tmp_mo, vdiff, vtmp, dipole)
 write(6,'(A,F20.8)') 'tot_change=', tot_change

 if(niter <= max_niter) then
  write(6,'(A)') 'Orbital localization converged successfully.'
 else
  write(6,'(/,A)') 'ERROR in subroutine serial2by2: max_niter exceeded.'
  write(6,'(A,I0)') 'max_niter=', max_niter
  stop
 end if
end subroutine serial2by2

! perform serial Berry 2-by-2 Jacobi rotations on given MOs with distance considered
subroutine serial22berry(natom, nbf, nmo, bfirst, dis, mo, mo_zdip)
 use lo_info, only: npair, max_niter, max_ncenter, dis_thres, conv_thres, &
  get_mo_center_by_scpa, atm_dis2mo_dis, find_eff_ijmap, serial22berry_kernel
 implicit none
 integer :: ncenter, niter
 integer, intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
 integer, allocatable :: mo_center(:,:)
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: mo(nbf,nmo)
 real(kind=8), allocatable :: mo_dis(:)
 complex(kind=8), intent(inout) :: mo_zdip(3,nmo,nmo)

 write(6,'(/,A)') 'Perform serial Berry 2*2 rotation...'
 write(6,'(A,F8.2)') 'dis_thres=', dis_thres
 ncenter = MIN(natom, max_ncenter)
 allocate(mo_center(0:ncenter,nmo), mo_dis(npair))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  call get_mo_center_by_scpa(natom, nbf, nmo, ncenter, bfirst, mo, mo_center)
  call atm_dis2mo_dis(natom, nmo, ncenter, dis, mo_center, mo_dis)
  call find_eff_ijmap(nmo, mo_dis)
  call serial22berry_kernel(nbf, nmo, mo, mo_zdip, sum_change)
  tot_change = tot_change + sum_change
  niter = niter + 1
  write(6,'(A,I4,A,F16.8)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for while

 write(6,'(A,F20.8)') 'tot_change=', tot_change
 deallocate(mo_center, mo_dis)

 if(niter <= max_niter) then
  write(6,'(A)') 'Orbital localization converged successfully.'
 else
  write(6,'(/,A)') 'ERROR in subroutine serial22berry: max_niter exceeded.'
  write(6,'(A,I0)') 'max_niter=', max_niter
  stop
 end if
end subroutine serial22berry

subroutine para22berry(natom, nbf, nmo, bfirst, dis, coeff, mo_zdip)
 use lo_info, only: np, npair, nsweep, max_niter, max_ncenter, dis_thres, &
  conv_thres, rrmap, rot_idx, get_mo_center_by_scpa, atm_dis2mo_dis, &
  screen_jacobi_idx_by_dis, para22berry_kernel
 implicit none
 integer :: ncenter, niter
 integer, intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
 integer, allocatable :: mo_center(:,:)
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: coeff(nbf,nmo)
 real(kind=8), allocatable :: mo_dis(:)
 complex(kind=8), intent(inout) :: mo_zdip(3,nmo,nmo)

 if(nmo < 4) then
  write(6,'(/,A)') 'ERROR in subroutine para22berry: nmo<4. Too few orbitals.'
  write(6,'(A,I0)') 'nmo=', nmo
  stop
 end if

 write(6,'(/,A)') 'Perform parallel Berry 2*2 rotation...'
 write(6,'(A,F8.2)') 'dis_thres=', dis_thres

 np = (nmo+1)/2
 nsweep = 4*np - 2
 allocate(rrmap(2,np,nsweep))
 call init_ring_jacobi_idx(nmo, np, rrmap)

 ncenter = MIN(natom, max_ncenter)
 allocate(mo_center(0:ncenter,nmo), mo_dis(npair))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  call get_mo_center_by_scpa(natom, nbf, nmo, ncenter, bfirst, coeff, mo_center)
  call atm_dis2mo_dis(natom, nmo, ncenter, dis, mo_center, mo_dis)
  call screen_jacobi_idx_by_dis(nmo, mo_dis)
  call para22berry_kernel(nbf, nmo, coeff, mo_zdip, sum_change)
  tot_change = tot_change + sum_change
  niter = niter + 1
  write(6,'(A,I4,A,F16.8)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for while

 deallocate(mo_center, mo_dis, rrmap, rot_idx)
 write(6,'(A,F20.8)') 'tot_change=', tot_change

 if(niter <= max_niter) then
  write(6,'(A)') 'Orbital localization converged successfully.'
 else
  write(6,'(/,A)') 'ERROR in subroutine para22berry: max_niter exceeded.'
  write(6,'(A,I0)') 'max_niter=', max_niter
  stop
 end if
end subroutine para22berry

! perform serial Boys 2-by-2 Jacobi rotations on given MOs with distance considered
subroutine serial22boys(natom, nbf, nmo, bfirst, dis, mo, mo_dip)
 use lo_info, only: npair, max_niter, max_ncenter, dis_thres, conv_thres, &
  get_mo_center_by_scpa, atm_dis2mo_dis, find_eff_ijmap, serial22boys_kernel
 implicit none
 integer :: ncenter, niter
 integer, intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
 integer, allocatable :: mo_center(:,:)
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: mo(nbf,nmo), mo_dip(3,nmo,nmo)
 real(kind=8), allocatable :: mo_dis(:)

 write(6,'(/,A)') 'Perform serial Boys 2*2 rotation...'
 write(6,'(A,F8.2)') 'dis_thres=', dis_thres
 ncenter = MIN(natom, max_ncenter)
 allocate(mo_center(0:ncenter,nmo), mo_dis(npair))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  call get_mo_center_by_scpa(natom, nbf, nmo, ncenter, bfirst, mo, mo_center)
  call atm_dis2mo_dis(natom, nmo, ncenter, dis, mo_center, mo_dis)
  call find_eff_ijmap(nmo, mo_dis)
  call serial22boys_kernel(nbf, nmo, mo, mo_dip, sum_change)
  tot_change = tot_change + sum_change
  niter = niter + 1
  write(6,'(A,I4,A,F16.8)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for while

 write(6,'(A,F20.8)') 'tot_change=', tot_change
 deallocate(mo_center, mo_dis)

 if(niter <= max_niter) then
  write(6,'(A)') 'Orbital localization converged successfully.'
 else
  write(6,'(/,A)') 'ERROR in subroutine serial22boys: max_niter exceeded.'
  write(6,'(A,I0)') 'max_niter=', max_niter
  stop
 end if
end subroutine serial22boys

! perform parallel Foster-Boys 2-by-2 Jacobi rotation on given MOs
subroutine para22boys(natom, nbf, nmo, bfirst, dis, coeff, mo_dip)
 use lo_info, only: np, npair, nsweep, max_niter, max_ncenter, dis_thres, &
  conv_thres, rrmap, rot_idx, get_mo_center_by_scpa, atm_dis2mo_dis, &
  screen_jacobi_idx_by_dis, para22boys_kernel
 implicit none
 integer :: ncenter, niter
 integer, intent(in) :: natom, nbf, nmo
 integer, intent(in) :: bfirst(natom+1)
 integer, allocatable :: mo_center(:,:)
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: coeff(nbf,nmo), mo_dip(3,nmo,nmo)
 real(kind=8), allocatable :: mo_dis(:)

 if(nmo < 4) then
  write(6,'(/,A)') 'ERROR in subroutine para22boys: nmo<4. Too few orbitals.'
  write(6,'(A,I0)') 'nmo=', nmo
  stop
 end if

 write(6,'(/,A)') 'Perform parallel Boys 2*2 rotation...'
 write(6,'(A,F8.2)') 'dis_thres=', dis_thres

 np = (nmo+1)/2
 nsweep = 4*np - 2
 allocate(rrmap(2,np,nsweep))
 call init_ring_jacobi_idx(nmo, np, rrmap)

 ncenter = MIN(natom, max_ncenter)
 allocate(mo_center(0:ncenter,nmo), mo_dis(npair))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  call get_mo_center_by_scpa(natom, nbf, nmo, ncenter, bfirst, coeff, mo_center)
  call atm_dis2mo_dis(natom, nmo, ncenter, dis, mo_center, mo_dis)
  call screen_jacobi_idx_by_dis(nmo, mo_dis)
  call para22boys_kernel(nbf, nmo, coeff, mo_dip, sum_change)
  tot_change = tot_change + sum_change
  niter = niter + 1
  write(6,'(A,I4,A,F16.8)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for while

 deallocate(mo_center, mo_dis, rrmap, rot_idx)
 write(6,'(A,F20.8)') 'tot_change=', tot_change

 if(niter <= max_niter) then
  write(6,'(A)') 'Orbital localization converged successfully.'
 else
  write(6,'(/,A)') 'ERROR in subroutine para22boys: max_niter exceeded.'
  write(6,'(A,I0)') 'max_niter=', max_niter
  stop
 end if
end subroutine para22boys

! perform serial Pipek-Mezey 2-by-2 Jacobi rotation on given MOs with distance considered
subroutine serial22pm(natom, nbf, nmo, dis, coeff, gross)
 use lo_info, only: npair, max_niter, max_ncenter, dis_thres, conv_thres, &
  get_mo_center_from_gross, atm_dis2mo_dis, find_eff_ijmap, serial22pm_kernel
 implicit none
 integer :: ncenter, niter
 integer, intent(in) :: natom, nbf, nmo
 integer, allocatable :: mo_center(:,:)
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: coeff(nbf,nmo), gross(natom,nmo,nmo)
 real(kind=8), allocatable :: mo_dis(:)

 write(6,'(/,A)') 'Perform serial PM 2*2 rotation...'
 write(6,'(A,F8.2)') 'dis_thres=', dis_thres
 ncenter = MIN(natom, max_ncenter)
 allocate(mo_center(0:ncenter,nmo), mo_dis(npair))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  call get_mo_center_from_gross(natom, nmo, ncenter, gross, mo_center)
  call atm_dis2mo_dis(natom, nmo, ncenter, dis, mo_center, mo_dis)
  call find_eff_ijmap(nmo, mo_dis)
  call serial22pm_kernel(nbf, nmo, natom, coeff, gross, sum_change)
  tot_change = tot_change + sum_change
  niter = niter + 1
  write(6,'(A,I4,A,F16.8)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for while

 deallocate(mo_center, mo_dis)
 write(6,'(A,F20.8)') 'tot_change=', tot_change

 if(niter <= max_niter) then
  write(6,'(A)') 'Orbital localization converged successfully.'
 else
  write(6,'(/,A)') 'ERROR in subroutine serial22pm: max_niter exceeded.'
  write(6,'(A,I0)') 'max_niter=', max_niter
  stop
 end if
end subroutine serial22pm

! perform parallel Pipek-Mezey 2-by-2 Jacobi rotation on given MOs
subroutine para22pm(natom, nbf, nmo, dis, coeff, gross)
 use lo_info, only: np, npair, nsweep, max_niter, max_ncenter, dis_thres, &
  conv_thres, rrmap, rot_idx, get_mo_center_from_gross, atm_dis2mo_dis, &
  screen_jacobi_idx_by_dis, para22pm_kernel
 implicit none
 integer :: ncenter, niter
 integer, intent(in) :: natom, nbf, nmo
 integer, allocatable :: mo_center(:,:)
 real(kind=8) :: sum_change, tot_change
 real(kind=8), intent(in) :: dis(natom,natom)
 real(kind=8), intent(inout) :: coeff(nbf,nmo), gross(natom,nmo,nmo)
 real(kind=8), allocatable :: mo_dis(:)

 if(nmo < 4) then
  write(6,'(/,A)') 'ERROR in subroutine para22pm: nmo<4. Too small number of or&
                   &bitals.'
  write(6,'(A,I0)') 'nmo=', nmo
  stop
 end if

 write(6,'(/,A)') 'Perform parallel PM 2*2 rotation...'
 write(6,'(A,F8.2)') 'dis_thres=', dis_thres

 np = (nmo+1)/2
 nsweep = 4*np - 2
 allocate(rrmap(2,np,nsweep))
 call init_ring_jacobi_idx(nmo, np, rrmap)

 ncenter = MIN(natom, max_ncenter)
 allocate(mo_center(0:ncenter,nmo), mo_dis(npair))
 tot_change = 0d0; niter = 0

 do while(niter <= max_niter)
  call get_mo_center_from_gross(natom, nmo, ncenter, gross, mo_center)
  call atm_dis2mo_dis(natom, nmo, ncenter, dis, mo_center, mo_dis)
  call screen_jacobi_idx_by_dis(nmo, mo_dis)
  call para22pm_kernel(nbf, nmo, natom, coeff, gross, sum_change)
  tot_change = tot_change + sum_change
  niter = niter + 1
  write(6,'(A,I4,A,F16.8)') 'niter=', niter, ', sum_change=', sum_change
  if(sum_change < conv_thres) exit
 end do ! for while

 deallocate(mo_center, mo_dis, rrmap, rot_idx)
 write(6,'(A,F20.8)') 'tot_change=', tot_change

 if(niter <= max_niter) then
  write(6,'(A)') 'Orbital localization converged successfully.'
 else
  write(6,'(/,A)') 'ERROR in subroutine para22pm: max_niter exceeded.'
  write(6,'(A,I0)') 'max_niter=', max_niter
  stop
 end if
end subroutine para22pm

! get the value of the modified Boys function
subroutine get_mboys(nif, mo_dipole)
 implicit none
 integer :: i
 integer, intent(in) :: nif
!f2py intent(in) :: nif
 real(kind=8), intent(in) :: mo_dipole(3,nif,nif)
!f2py intent(in) :: mo_dipole
!f2py depend(nif) :: mo_dipole
 real(kind=8) :: ddot, fBoys, tmp_dip(3)

 fBoys = 0d0
 do i = 1, nif, 1
  tmp_dip = mo_dipole(:,i,i)
  fBoys = fBoys + ddot(3, tmp_dip, 1, tmp_dip, 1)
 end do ! for i

 fBoys = DSQRT(fBoys/DBLE(nif))
 write(6,'(A,F20.7)') 'Modified f(Boys)=', fBoys
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
 real(kind=8), intent(in) :: mo_coeff(nbf,nmo), mo_dipole(3,nmo,nmo)
!f2py intent(in) :: mo_coeff, mo_dipole
!f2py depend(nbf,nmo) :: mo_coeff
!f2py depend(nmo) :: mo_dipole
 real(kind=8), intent(out) :: new_coeff(nbf,nmo)
!f2py intent(out) :: new_coeff
!f2py depend(nbf,nmo) :: new_coeff
 real(kind=8), allocatable :: f(:,:), w(:)

 allocate(f(nmo,nmo), source=0d0)
 call dsymm('L','L',nmo,nmo,1d0,mo_dipole(1,:,:),nmo,mo_dipole(1,:,:),nmo,0d0,f,nmo)
 call dsymm('L','L',nmo,nmo,1d0,mo_dipole(2,:,:),nmo,mo_dipole(2,:,:),nmo,1d0,f,nmo)
 call dsymm('L','L',nmo,nmo,1d0,mo_dipole(3,:,:),nmo,mo_dipole(3,:,:),nmo,1d0,f,nmo)

 allocate(w(nmo), source=0d0)
 call diag_get_e_and_vec(nmo, f, w)
 write(6,'(5(1X,ES15.8))') w
 stop
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
subroutine boys_noiter(nbf, nmo, ao_ovlp, mo, lmo, old_dip)
 implicit none
 integer :: i, j, k, m, p, q, nmo_s
 integer, intent(in) :: nbf, nmo
!f2py intent(in) :: nbf, nmo
 real(kind=8) :: rtmp(3)
 real(kind=8), intent(in) :: ao_ovlp(nbf,nbf), mo(nbf,nmo), lmo(nbf,nmo), &
  old_dip(3,nmo,nmo)
!f2py intent(in) :: ao_ovlp, mo, lmo, old_dip
!f2py depend(nbf,nmo) :: mo, lmo
!f2py depend(nbf) :: ao_ovlp
!f2py depend(nmo) :: old_dip
 real(kind=8), allocatable :: y(:,:), u(:,:), x(:,:)

 nmo_s = nmo*nmo
 allocate(y(nmo_s,nmo_s), source=0d0)

 ! diagonal elements
 do i = 1, nmo_s, 1
  call idx_map_2d(nmo, i, k, m)
  rtmp = old_dip(:,m,k)
  y(i,i) = DOT_PRODUCT(rtmp, rtmp)
 end do ! for i

 ! non-diagonal elements
 do i = 1, nmo_s-1, 1
  call idx_map_2d(nmo, i, k, m)
  rtmp = old_dip(:,m,k)
  do j = i+1, nmo_s, 1
   call idx_map_2d(nmo, j, p, q)
   y(j,i) = DOT_PRODUCT(rtmp, old_dip(:,q,p))
   y(i,j) = y(j,i)
  end do ! for j
 end do ! for i

 allocate(u(nmo,nmo))
 call calc_CTSCp(nbf, nmo, mo, ao_ovlp, lmo, u)
 allocate(x(nmo_s,nmo))
 do i = 1, nmo_s, 1
  call idx_map_2d(nmo, i, k, m)
  forall(j = 1:nmo) x(i,j) = u(k,j)*u(m,j)
 end do ! for i
 write(6,'(/,A)') 'X_iik'
 do i = 1, nmo, 1
  do j = 1, nmo, 1
   k = (j-1)*nmo + j
   write(6,'(2I5,F20.8)') j,i,x(k,i)
  end do ! for j
 end do ! for i

 call calc_CTSC(nmo_s, nmo, x, y, u)
 deallocate(x, y)

 write(6,'(/,A)') '(X^T)YX='
 do i = 1, nmo, 1
  write(6,'(I5,F20.8)') i, u(i,i)
 end do
 do i = 1, nmo-1, 1
  do j = i+1, nmo, 1
   write(6,'(2I5,2F20.8)') i, j, u(j,i), u(i,j)
  end do
 end do
 deallocate(u)
end subroutine boys_noiter

