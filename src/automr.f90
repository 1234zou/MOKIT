! written by jxzou at 20200420: black-box multireference calculations
! updated by jxzou at 20200511: framework of the program
! updated by jxzou at 20201213: %casscf+noiter -> %mrci for correct CASCI NOONs in ORCA
! updated by jxzou at 20201222: import grad for CASSCF force in PySCF; read CASSCF force for Molpro
! updated by jxzou at 20210111: add subroutine do_mrpt3
! updated by jxzou at 20210119: add OpenMolcas-QCMaquis, OpenMolcas-CheMPS2 interface

! The input file is just like Gaussian .gjf format. MOKIT keywords should be
!  specified in the Title Card line like 'mokit{}'.

program main
 use mr_keyword, only: read_program_path
 implicit none
 integer :: i, j
 character(len=240) :: fname = ' '

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine automr: wrong command line argument!'
  write(6,'(A)')   " Example: automr h2o.gjf >h2o.out 2>&1 &"
  write(6,'(A,/)') ' See help: automr -h'
  stop
 end if

 call getarg(1, fname)

 select case(TRIM(fname))
 case('-v', '-V', '--version')
  write(6,'(A)') 'AutoMR 1.2.6rc30 :: MOKIT, release date: 2024-May-8'
  stop
 case('-h','-help','--help')
  write(6,'(/,A)') 'Usage: automr [gjfname] > [outname]'
  write(6,'(A)')   "  Example: automr h2o.gjf >h2o.out 2>&1 &"
  write(6,'(/,A)') 'Options:'
  write(6,'(A)')   '  -h, -help, --help: Print this message and exit.'
  write(6,'(A)')   '  -v, -V, --version: Print the version number of automr and exit.'
  write(6,'(A)')   '  -t, --testprog: Print the path of programs detected by automr and exit.'
  write(6,'(/,A)') 'Methods(#p ...):'
  write(6,'(A)')   '  GVB, CASCI, CASSCF, DMRGCI, DMRGSCF, NEVPT2, NEVPT3,&
                   & CASPT2, CASPT2-K,'
  write(6,'(A)')   '  CASPT3, MRMP2, OVBMP2, SDSPT2, MRCISD, MRCISDT, MCPDFT,&
                   & FICMRCCSD,'
  write(6,'(A)')   '  MkMRCCSD, MkMRCCSD(T), BWMRCCSD, BWMRCCSD(T), BCCC2b, BCCC3b'
  write(6,'(/,A)') 'Frequently used keywords in MOKIT{}:'
  write(6,'(A)')   '      HF_prog=Gaussian/PySCF/ORCA/PSI4'
  write(6,'(A)')   '     GVB_prog=GAMESS/Gaussian/QChem'
  write(6,'(A)')   '  CASSCF_prog=PySCF/OpenMolcas/ORCA/Molpro/GAMESS/Gaussian/BDF/PSI4/Dalton'
  write(6,'(A)')   '   CASCI_prog=PySCF/OpenMolcas/ORCA/Molpro/GAMESS/Gaussian/BDF/PSI4/Dalton'
  write(6,'(A)')   '  NEVPT2_prog=PySCF/OpenMolcas/ORCA/Molpro/BDF'
  write(6,'(A)')   '  CASPT2_prog=OpenMolcas/Molpro/ORCA'
  write(6,'(A)')   '  MCPDFT_prog=OpenMolcas/PySCF/GAMESS'
  write(6,'(A)')   '  MRCISD_prog=OpenMolcas/Molpro/ORCA/Gaussian/GAMESS/PSI4/Dalton'
  write(6,'(A)')   '      CtrType=1/2/3 for uc-/ic-/FIC-MRCISD'
  write(6,'(A,/)') '    MRCC_prog=ORCA'
  stop
 case('-t','--testprog')
  call check_mokit_root()
  call read_program_path()
  stop
 end select

 i = INDEX(fname, '.gjf', back=.true.)
 j = INDEX(fname, '.fch', back=.true.)
 if(i>0 .and. j>0) then
  write(6,'(/,A)') "ERROR in subroutine automr: both '.gjf' and '.fch' keys det&
                   &ected in"
  write(6,'(A)') 'filename '//TRIM(fname)//'.'
  write(6,'(A)') "Better to use a filename only with suffix '.gjf'."
  stop
 else if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine automr: '.gjf' key not found in filenam&
                   &e "//TRIM(fname)
  stop
 end if

 call require_file_exist(fname)
 call automr(fname)
end program main

! automatically do multireference calculations in a block-box way
subroutine automr(fname)
 use mr_keyword, only: gjfname, read_program_path, parse_keyword, check_kywd_compatible
 implicit none
 character(len=24) :: data_string
 character(len=240), intent(in) :: fname

 gjfname = fname
 ! read paths of various programs from environment variables
 call check_mokit_root()
 call read_program_path()
 call parse_keyword()
 call check_kywd_compatible()

 call do_hf(.true.)          ! RHF and/or UHF
 call do_minimal_basis_gvb() ! GVB/STO-6G, only valid for ist=6
 call get_paired_LMO()
 call do_gvb()        ! GVB
 call do_cas(.false.) ! CASCI/DMRG-CASCI
 call do_cas(.true.)  ! CASSCF/DMRG-CASSCF, including SS-CASSCF
 call do_mrpt2()      ! CASPT2/NEVPT2/SDSPT2/MRMP2
 call do_mrpt3()      ! CASPT3/NEVPT3
 call do_mrcisd()     ! uncontracted/ic-/FIC- MRCISD
 call do_mrcisdt()    ! uncontracted MRCISDT
 call do_mcpdft()     ! MC-PDFT
 call do_mrcc()       ! MRCC

 call do_cis()        ! CIS/TDHF
 call do_sa_cas()     ! SA-CASSCF
 call do_pes_scan()   ! PES scan

 call fdate(data_string)
 write(6,'(/,A)') 'Normal termination of AutoMR at '//TRIM(data_string)
end subroutine automr

! check if MOKIT_ROOT present. Set it if not present.
subroutine check_mokit_root()
 integer :: i, SYSTEM
 character(len=240) :: mokit_root

 write(6,'(/,A)') 'Checking MOKIT_ROOT ... '
 mokit_root = ' '
 call getenv('MOKIT_ROOT', mokit_root)
 if (len_trim(mokit_root) < 1) then
   ! assume we are under conda install
   i = SYSTEM("echo `get_mokit_loc.py` > ~/.mokitrc ")
   !i = SYSTEM("echo `get_mokit_loc.py` ")
   !i = SYSTEM('echo $MOKIT_ROOT')
   if(i /= 0) then
    write(6,'(/,A)') 'ERROR in subroutine check_mokit_root: '
    write(6,'(A)') '    fail to set MOKIT_ROOT for conda-installed MOKIT.'
    write(6,'(A)') 'If MOKIT is installed via conda, please report this issue.'
    write(6,'(A)') 'Otherwise, it means your MOKIT_ROOT is not properly set.'
    stop
   else
    write(6,'(A)') 'reset MOKIT_ROOT for conda-installed version'
   end if
 endif
end subroutine check_mokit_root

! generate PySCF input file .py from Gaussian .fch(k) file, and get paired LMOs
subroutine get_paired_LMO()
 use mr_keyword, only: eist, mo_rhf, ist, hf_fch, bgchg, chgname, nskip_uno
 use mol, only: nbf, nif, ndb, nacte, nacto, nacta, nactb, npair, npair0, nopen,&
  lin_dep, chem_core, ecp_core
 use util_wrapper, only: bas_fch2py_wrap
 implicit none
 integer :: i, SYSTEM
 real(kind=8) :: unpaired_e
 character(len=24) :: data_string = ' '
 character(len=240) :: proname, pyname, outname, fchname

 if(eist == 1) return ! excited state calculation
 if(ist == 5) return ! no need for this subroutine
 write(6,'(//,A)') 'Enter subroutine get_paired_LMO...'

 if(ist == 0) then
  write(6,'(A)') 'ERROR in subroutine get_paired_LMO: ist=0. It should be non-&
                 &zero before this subroutine.'
  stop
 end if

 ! calculate the number of core orbitals from array core_orb
 call calc_ncore(hf_fch, chem_core, ecp_core)
 write(6,'(3(A,I0))') 'chem_core=', chem_core, ', ecp_core=', ecp_core, &
                      ', Nskip_UNO=', nskip_uno

 i = INDEX(hf_fch, '.fch', back=.true.)
 proname = hf_fch(1:i-1)

 if(mo_rhf) then
  if(ist == 6) then
   write(6,'(A)') 'One set of MOs: GVB/STO-6G -> MO projection -> Rotate has be&
                  &en invoked.'
  else
   write(6,'(A)') 'One set of MOs: invoke RHF virtual MO projection -> localiza&
                  &tion -> paring.'
   pyname = TRIM(proname)//'_proj_loc_pair.py'
   call prt_rhf_proj_script_into_py(pyname)
   call prt_auto_pair_script_into_py(pyname)
   call submit_pyscf_job(pyname)
  end if

 else
  if(ist == 1) then
   write(6,'(A)') 'Two sets of MOs, ist=1, invoke UNO associated rotation.'
  else if(ist == 2) then
   write(6,'(A)') 'Two sets of MOs, ist=2, invoke UNO generation.'
  end if

  fchname = hf_fch(1:i-1)//'_uno.fch'
  if(ist == 1) then
   pyname = TRIM(proname)//'_uno_asrot.py'
   outname = TRIM(proname)//'_uno_asrot.out'
  else
   pyname = TRIM(proname)//'_uno.py'
   outname = TRIM(proname)//'_uno.out'
  end if
  call bas_fch2py_wrap(hf_fch, .false., pyname)
  call prt_uno_script_into_py(pyname)
  if(ist == 1) call prt_assoc_rot_script_into_py(pyname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(pyname))
  call submit_pyscf_job(pyname)
  call calc_unpaired_from_fch(fchname, 1, .false., unpaired_e)

  ! when ist=2, GVB will not be performed, so we need to read variables before CASCI
  if(ist == 2) then
   call read_npair_from_uno_out(nbf, nif, ndb, npair, nopen, lin_dep)
   ! find npair0: the number of active orbitals (NOON <- [0.02,0.98])
   call find_npair0_from_fch(fchname, nopen, npair0)
   ! determine the number of orbitals/electrons in following CAS/DMRG computations
   nacta = npair0 + nopen
   nactb = npair0
   nacte = nacta + nactb
   nacto = nacte
  end if
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine get_paired_LMO at '//TRIM(data_string)
end subroutine get_paired_LMO

! print RHF virtual MOs projection scheme into a given .py file
subroutine prt_rhf_proj_script_into_py(pyname)
 use mr_keyword, only: hf_prog, hf_fch
 use mol, only: nuc, beyond_xe
 implicit none
 integer :: i, fid
 character(len=240) :: proj_fch, proname, minbas_fch
 character(len=240), intent(in) :: pyname
 logical :: beyond_kr

 call find_specified_suffix(hf_fch, '.fch', i)
 proj_fch = hf_fch(1:i-1)//'_proj.fch'
 proname = hf_fch(1:i-1)//'_minb'
 minbas_fch = hf_fch(1:i-1)//'_minb.fch'

 if(ANY(nuc > 54)) beyond_xe = .true.
 beyond_kr = .false.
 if(ANY(nuc > 36)) beyond_kr = .true.

 ! Perform a HF/STO-6G calculation using Gaussian/PySCF.
 ! Note: STO-6G in PySCF ranges from H~Kr, but in Gaussian ranges from H~Xe,
 !  so Gaussian is preferred.
 select case(TRIM(hf_prog))
 case('gaussian')
  if(beyond_xe) then
   return
  else
   call do_min_bas_hf_gau(proname)
  end if
 case('pyscf')
  if(beyond_kr) then
   beyond_xe = .true.
   return
  else
   call do_min_bas_hf_pyscf(proname)
  end if
 end select

 open(newunit=fid,file=TRIM(pyname),status='replace')
 write(fid,'(A)') 'from pyscf import gto'
 write(fid,'(A)') 'from mokit.lib.gaussian import load_mol_from_fch'
 write(fid,'(A)') 'from mokit.lib.rwwfn import read_nbf_and_nif_from_fch, read_&
                  &na_and_nb_from_fch, get_occ_from_na_nb'
 write(fid,'(A)') 'from mokit.lib.mo_svd import mo_svd'
 write(fid,'(A)') 'from mokit.lib.fch2py import fch2py'
 write(fid,'(A)') 'from mokit.lib.py2fch import py2fch'
 write(fid,'(A,/)') 'from shutil import copyfile'

 write(fid,'(A)') "hf_fch = '"//TRIM(hf_fch)//"'"
 write(fid,'(A)') 'mol = load_mol_from_fch(hf_fch)'
 write(fid,'(A)') 'nbf, nif = read_nbf_and_nif_from_fch(hf_fch)'
 write(fid,'(A,/)') "mo = fch2py(hf_fch, nbf, nif, 'a')"

 write(fid,'(A)') "minbas_fch = '"//TRIM(minbas_fch)//"'"
 write(fid,'(A)') 'mol2 = load_mol_from_fch(minbas_fch)'
 write(fid,'(A)') 'nbf2, nif2 = read_nbf_and_nif_from_fch(minbas_fch)'
 write(fid,'(A)') "mo2 = fch2py(minbas_fch, nbf2, nif2, 'a')"
 write(fid,'(/,A)') '# project virtual MOs onto those of R(O)HF/STO-6G'
 write(fid,'(A)') 'na, nb = read_na_and_nb_from_fch(hf_fch)'
 write(fid,'(A)') 'na2, nb2 = read_na_and_nb_from_fch(minbas_fch)'
 write(fid,'(A)') 'npair = nif2 - na2'
 write(fid,'(A)') 'coeff1 = mo[:,range(na,nif)]'
 write(fid,'(A)') 'coeff2 = mo2[:,range(na2,nif2)]'
 write(fid,'(A)') "cross_S = gto.intor_cross('int1e_ovlp', mol, mol2)"
 write(fid,'(A)') 'mo_svd(nbf, nif-na, nbf2, npair, coeff1, coeff2, cross_S, False)'
 write(fid,'(A)') 'mo[:,range(na,nif)] = coeff1.copy()'
 write(fid,'(A)') '# project done'
 write(fid,'(/,A)') '# save projected MOs into a new .fch file'
 write(fid,'(A)') "proj_fch = '"//TRIM(proj_fch)//"'"
 write(fid,'(A)') "copyfile(hf_fch, proj_fch)"
 write(fid,'(A)') 'mo_occ = get_occ_from_na_nb(nif, na, nb)'
 write(fid,'(A)') "py2fch(proj_fch, nbf, nif, mo, 'a', mo_occ, False, False)"
 write(fid,'(A)') '# save done'
 close(fid)
end subroutine prt_rhf_proj_script_into_py

! print localization and automatically pairing information into a given .py file
subroutine prt_auto_pair_script_into_py(pyname)
 use mr_keyword, only: localm, hf_fch
 use mol, only: beyond_xe, chem_core, ecp_core
 implicit none
 integer :: i, ncore, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1, loc_fch, lmo_fch
 character(len=240), intent(in) :: pyname

 buf = ' '
 i = INDEX(pyname, '.py', back=.true.)
 pyname1 = pyname(1:i-1)//'.t'
 i = INDEX(hf_fch, '.fch', back=.true.)
 loc_fch = hf_fch(1:i-1)//'_proj_loc_pair.fch'
 lmo_fch = hf_fch(1:i-1)//'_LMO.fch'
 ncore = chem_core - ecp_core

 if(beyond_xe) then ! some element >Xe
  open(newunit=fid1,file=TRIM(pyname),status='replace')
  write(fid1,'(A)') 'from pyscf.lo.boys import dipole_integral'
  write(fid1,'(A)') 'from mokit.lib.gaussian import load_mol_from_fch, loc'
  write(fid1,'(A)') 'from mokit.lib.rwwfn import read_nbf_and_nif_from_fch, rea&
                    &d_na_and_nb_from_fch, \'
  write(fid1,'(A)') ' read_eigenvalues_from_fch, get_occ_from_na_nb, get_core_&
                     &valence_sep_idx'
  write(fid1,'(A)') 'from mokit.lib.rwwfn import get_1e_exp_and_sort_pair as sor&
                    &t_pair'
  write(fid1,'(A)') 'from mokit.lib.assoc_rot import assoc_loc'
  write(fid1,'(A)') 'from mokit.lib.fch2py import fch2py'
  write(fid1,'(A)') 'from mokit.lib.py2fch import py2fch'
  write(fid1,'(A)') 'from shutil import copyfile'
  write(fid1,'(/,A)') "hf_fch = '"//TRIM(hf_fch)//"'"
  write(fid1,'(A)') "lmo_fch = '"//TRIM(lmo_fch)//"'"
  write(fid1,'(A,I0)') 'ncore = ', ncore
  write(fid1,'(A)') 'na, nb = read_na_and_nb_from_fch(hf_fch)'
  write(fid1,'(A)') 'nopen = na - nb'
  write(fid1,'(A)',advance='no') 'loc(fchname=hf_fch, idx=range(ncore,nb)'
  if(TRIM(localm) == 'boys') write(fid1,'(A)',advance='no') ", method='boys'"
  write(fid1,'(A)') ')'
  write(fid1,'(A)') 'sort_pair(lmo_fch, hf_fch, nb-ncore)'
  write(fid1,'(/,A)') 'mol = load_mol_from_fch(lmo_fch)'
  write(fid1,'(A)') 'nbf, nif = read_nbf_and_nif_from_fch(lmo_fch)'
  write(fid1,'(A)') "mo = fch2py(lmo_fch, nbf, nif, 'a')"
  write(fid1,'(A)') 'mo_dipole = dipole_integral(mol, mo)'
  write(fid1,'(A)') 'ref1 = get_core_valence_sep_idx(lmo_fch) - 1'
  write(fid1,'(A)') 'npair = nb - ref1'
  write(fid1,'(A)') 'new_mo = assoc_loc(nbf, nif, ref1, nb, na, nif, mo, mo_dipole)'
  write(fid1,'(/,A)') "loc_fch = '"//TRIM(loc_fch)//"'"
  write(fid1,'(A)') 'copyfile(lmo_fch, loc_fch)'
  write(fid1,'(A)') 'mo_occ = get_occ_from_na_nb(nif, na, nb)'
  write(fid1,'(A)') "py2fch(loc_fch, nbf, nif, new_mo, 'a', mo_occ, False, False)"
 else               ! all elements <= Xe
  open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
  open(newunit=fid2,file=TRIM(pyname1),status='replace')
  do while(.true.)
   read(fid1,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(LEN_TRIM(buf) == 0) exit
   write(fid2,'(A)') TRIM(buf)
  end do

  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine prt_auto_pair_script_into_py: end-of-fi&
                    &le detected.'
   write(6,'(A)') 'File may be incomplete: '//TRIM(pyname)
   close(fid1)
   close(fid2,status='delete')
   stop
  end if

  write(fid2,'(A)') 'from mokit.lib.rwwfn import get_1e_exp_and_sort_pair as sor&
                    &t_pair'
  if(TRIM(localm) == 'pm') then ! Pipek-Mezey localization
   write(fid2,'(A)') 'from mokit.lib.lo import pm'
  else                          ! Foster-Boys localization
   write(fid2,'(A)') 'from mokit.lib.lo import boys'
  end if
  write(fid2,'(A)') 'from pyscf.lo.boys import dipole_integral'
  write(fid2,'(A,/)') 'from mokit.lib.auto_pair import pair_by_tdm'

  do while(.true.)
   read(fid1,'(A)',iostat=i) buf
   if(i /= 0) exit
   write(fid2,'(A)') TRIM(buf)
  end do
  close(fid1,status='delete')

  write(fid2,'(/,A)') '# localize the projected MOs of larger basis'
  write(fid2,'(A,I0)') 'ncore = ', ncore
  write(fid2,'(A)') 'nopen = na - nb'
  write(fid2,'(A)') 'nvir_lmo = npair # backup'
  write(fid2,'(A)') 'npair = min(npair, nb-ncore)'
  write(fid2,'(A)') 'occ_idx = range(ncore,nb)'
  write(fid2,'(A)') 'vir_idx = range(na,nif2)'

  if(TRIM(localm) == 'pm') then ! Pipek-Mezey localization
   write(fid2,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
   write(fid2,'(A)') 'occ_loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],mol.&
                     &_bas[:,3],mol.cart,nbf,'
   write(fid2,'(17X,A)') "nb-ncore,mo[:,occ_idx],S,'mulliken')"
   write(fid2,'(A)') 'vir_loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],mol.&
                     &_bas[:,3],mol.cart,nbf,'
   write(fid2,'(17X,A)') "nif2-na,mo[:,vir_idx],S,'mulliken')"
  else                          ! Foster-Boys localization
   write(fid2,'(A)') 'mo_dipole = dipole_integral(mol, mo[:,occ_idx])'
   write(fid2,'(A)') 'occ_loc_orb = boys(nbf, nb-ncore, mo[:,occ_idx], mo_dipole)'
   write(fid2,'(A)') 'mo_dipole = dipole_integral(mol, mo[:,vir_idx])'
   write(fid2,'(A)') 'vir_loc_orb = boys(nbf, nif2-na, mo[:,vir_idx], mo_dipole)'
  end if

  write(fid2,'(A)') 'mo[:,occ_idx] = occ_loc_orb.copy()'
  write(fid2,'(A)') 'mo[:,vir_idx] = vir_loc_orb.copy()'
  write(fid2,'(A)') '# localization done'
  write(fid2,'(/,A)') '# pair the active orbitals'
  write(fid2,'(A)') 'mo_dipole = dipole_integral(mol, mo)'
  write(fid2,'(A)') 'alpha_coeff = pair_by_tdm(ncore,npair,nopen,na,nvir_lmo,nb&
                    &f,nif,mo,mo_dipole)'
  write(fid2,'(A)') 'mo = alpha_coeff.copy()'
  write(fid2,'(A)') '# pair done'
  write(fid2,'(/,A)') '# save the paired LMO into .fch file'
  write(fid2,'(A)') "loc_fch = '"//TRIM(loc_fch)//"'"
  write(fid2,'(A)') "copyfile(hf_fch, loc_fch)"
  write(fid2,'(A)') "py2fch(loc_fch, nbf, nif, mo, 'a', mo_occ, False, False)"
  write(fid2,'(A)') "sort_pair(loc_fch, hf_fch, npair)"
  write(fid2,'(A)') '# save done'
  close(fid2)
  i = RENAME(TRIM(pyname1), TRIM(pyname))
  open(newunit=fid1,file=TRIM(pyname),status='old',position='append')
 end if

 write(fid1,'(/,A)') "f = open('uno.out', 'w+')"
 write(fid1,'(A)') "f.write('nbf=%i\n' %nbf)"
 write(fid1,'(A)') "f.write('nif=%i\n\n' %nif)"
 write(fid1,'(A)') 'idx1 = nb - npair'
 write(fid1,'(A)') "f.write('ndb=%i\n\n' %idx1)"
 write(fid1,'(A)') "f.write('idx=%i %i %i' %(idx1+1,na+npair+1,nopen))"
 write(fid1,'(A)') 'f.close()'
 close(fid1)
end subroutine prt_auto_pair_script_into_py

! print UNO script into a given .py file
subroutine prt_uno_script_into_py(pyname)
 use mr_keyword, only: nproc, hf_fch, uno_thres
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1, uno_fch
 character(len=240), intent(in) :: pyname

 buf = ' '
 pyname1 = TRIM(pyname)//'.t'
 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine prt_uno_script_into_py: end-of-file det&
                   &ected.'
  write(6,'(A)') 'File may be problematic: '//TRIM(pyname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A)') 'from pyscf import lib'
 write(fid2,'(A)') 'import numpy as np'
 write(fid2,'(A)') 'import os'
 write(fid2,'(A)') 'from mokit.lib.py2fch import py2fch'
 write(fid2,'(A)') 'from mokit.lib.uno import uno'
 write(fid2,'(A)') 'from mokit.lib.gaussian import mo_fch2py'
 write(fid2,'(A)') 'from mokit.lib.rwwfn import read_nbf_and_nif_from_fch, \'
 write(fid2,'(A)') '                            read_na_and_nb_from_fch'
 write(fid2,'(A,/)') 'from mokit.lib.construct_vir import construct_vir'
 write(fid2,'(A,I0)') 'nproc = ', nproc
 write(fid2,'(A,/)') 'lib.num_threads(nproc)'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:9) == 'mf.kernel') exit
  write(fid2,'(A)') TRIM(buf)
 end do

 write(fid2,'(/,A)') '# read MOs from .fch(k) file'
 write(fid2,'(A)') "hf_fch = '"//TRIM(hf_fch)//"'"
 write(fid2,'(A)') 'mf.mo_coeff = mo_fch2py(hf_fch)'
 write(fid2,'(A)') 'nbf, nif = read_nbf_and_nif_from_fch(hf_fch)'
 write(fid2,'(A)') '# read done'
 write(fid2,'(/,A)') '# check if input MOs are orthonormal'
 write(fid2,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
 write(fid2,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff[0], S)'
 write(fid2,'(A)') 'check_orthonormal(nbf, nif, mf.mo_coeff[1], S)'

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
 i = INDEX(hf_fch, '.fch', back=.true.)
 uno_fch = hf_fch(1:i-1)//'_uno.fch'

 open(newunit=fid1,file=TRIM(pyname),status='old',position='append')
 write(fid1,'(/,A)') "uno_fch = '"//TRIM(uno_fch)//"'"
 write(fid1,'(/,A)') '# transform UHF canonical orbitals to UNO'
 write(fid1,'(A)') 'na, nb = read_na_and_nb_from_fch(hf_fch)'
 write(fid1,'(A,E12.5,A)') 'idx, noon, alpha_coeff = uno(nbf,nif,na,nb,mf.mo_coeff[0],&
                           &mf.mo_coeff[1],S,',uno_thres,')'
 write(fid1,'(A)') 'alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)'
 write(fid1,'(A)') 'mf.mo_coeff = (alpha_coeff, alpha_coeff)'
 write(fid1,'(A)') '# done transform'

 write(fid1,'(/,A)') '# save the UNO into .fch file'
 !write(fid1,'(A)') "os.system('fch_u2r "//TRIM(hf_fch)//' '//TRIM(uno_fch)//"')"
 !write(fid1,'(A)') 'sleep(1) # in some node, py2fch begins when fch_u2r unfinished'
 ! Thanks to the suggestion of Kalinite. The problem in the above two lines
 ! does not exist now.
 write(fid1,'(A)') "with os.popen('fch_u2r '+hf_fch+' '+uno_fch) as run:"
 write(fid1,'(A)') '  null = run.read()'
 write(fid1,'(A)') "py2fch(uno_fch,nbf,nif,mf.mo_coeff[0],'a',noon,True,True)"
 write(fid1,'(A)') '# save done'
 close(fid1)
end subroutine prt_uno_script_into_py

! print associated rotation into a given .py file
subroutine prt_assoc_rot_script_into_py(pyname)
 use mol, only: chem_core, ecp_core
! use mol, only: natom, elem, nuc, chem_core, ecp_core
 use mr_keyword, only : localm, hf_fch, npair_wish, nskip_uno
 implicit none
 integer :: i, ncore, fid1, fid2, RENAME
! integer, allocatable :: ntimes(:)
 character(len=240) :: buf, pyname1, assoc_fch
 character(len=240), intent(in) :: pyname

 ncore = chem_core - ecp_core
 pyname1 = TRIM(pyname)//'.t'
 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine prt_assoc_rot_script_into_py: end-of-fi&
                   &le detected.'
  write(6,'(A)') 'File may be incomplete: '//TRIM(pyname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 if(localm == 'pm') then
  write(fid2,'(A)') 'from mokit.lib.lo import pm, boys'
 else
  write(fid2,'(A)') 'from mokit.lib.lo import boys'
 end if
 write(fid2,'(A)') 'from pyscf.lo.boys import dipole_integral'
 write(fid2,'(A)') 'from mokit.lib.auto_pair import pair_by_tdm'
 write(fid2,'(A)') 'from mokit.lib.assoc_rot import assoc_rot'
 write(fid2,'(A)') 'from mokit.lib.mo_svd import proj_occ_get_act_vir'
 write(fid2,'(A)') 'from mokit.lib.rwwfn import get_1e_exp_and_sort_pair as sort_pair'
 write(fid2,'(A,/)') 'from shutil import copyfile'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')
 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 i = INDEX(hf_fch, '.fch', back=.true.)
 assoc_fch = hf_fch(1:i-1)//'_uno_asrot.fch'

 open(newunit=fid1,file=TRIM(pyname),status='old',position='append')
 write(fid1,'(/,A)') "assoc_fch = '"//TRIM(assoc_fch)//"'"
 write(fid1,'(/,A)') '# associated rotation'
 write(fid1,'(A)') 'npair = np.int64((idx[1]-idx[0]-idx[2])/2)'
 write(fid1,'(A)') 'if(npair > 0):'
 write(fid1,'(A)') '  idx2 = idx[0] + npair - 1'
 write(fid1,'(A)') '  idx3 = idx2 + idx[2]'
 if(npair_wish > 0) write(fid1,'(A,I0,A1)') '  npair = min(npair,',npair_wish,')'
 write(fid1,'(A)') '  idx1 = idx2 - npair'
 write(fid1,'(A)') '  idx4 = idx3 + npair'
 write(fid1,'(A,I0,A)') '  i = ',nskip_uno,' # pair(s) of UNO to be skipped'
 write(fid1,'(A)') '  occ_idx = range(idx1,idx2-i)'
 write(fid1,'(A)') '  vir_idx = range(idx3+i,idx4)'

 if(localm == 'pm') then ! Pipek-Mezey localization
  write(fid1,'(A)') "  occ_loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],mo&
                    &l._bas[:,3],mol.cart,nbf,"
  write(fid1,'(19X,A)') "npair-i,mf.mo_coeff[0][:,occ_idx],S,'mulliken')"
 else ! Boys localization
  write(fid1,'(A)') '  mo_dipole = dipole_integral(mol, mf.mo_coeff[0][:,occ_id&
                      &x])'
  write(fid1,'(A)') '  occ_loc_orb = boys(nbf, npair-i, mf.mo_coeff[0][:,occ_id&
                      &x], mo_dipole)'
 end if

 write(fid1,'(A)') '  vir_loc_orb = assoc_rot(nbf, npair-i, mf.mo_coeff[0][:,oc&
                     &c_idx], occ_loc_orb,'
 write(fid1,'(26X,A)') 'mf.mo_coeff[0][:,vir_idx])'
 write(fid1,'(A)') '  mf.mo_coeff[0][:,occ_idx] = occ_loc_orb.copy()'
 write(fid1,'(A)') '  mf.mo_coeff[0][:,vir_idx] = vir_loc_orb.copy()'
 write(fid1,'(A)') '# localization done'

! write(fid1,'(/,A)') '# using projection to supplement pairs (if needed)'

 write(fid1,'(/,A)') '# save associated rotation MOs into .fch(k) file'
 write(fid1,'(A)') 'copyfile(uno_fch, assoc_fch)'
 write(fid1,'(A)') 'noon = np.zeros(nif)'
 write(fid1,'(A)') "py2fch(assoc_fch, nbf, nif, mf.mo_coeff[0], 'a', noon, Fals&
                   &e, False)"
 write(fid1,'(A)') 'sort_pair(assoc_fch, uno_fch, npair)'
 close(fid1)
end subroutine prt_assoc_rot_script_into_py

! perform GVB/STO-6G, only valid for ist=6
subroutine do_minimal_basis_gvb()
 use mol, only: mult, nbf, nif, nopen, ndb, npair
 use mr_keyword, only: nproc, ist, npair_wish, gjfname, localm, hf_fch, mo_rhf,&
  nskip_uno, bgchg, inherit
 implicit none
 integer :: i, fid, SYSTEM
 real(kind=8) :: e(3), uhf_s2 ! RHF/UHF/GVB energies and UHF spin mult
 real(kind=8) :: gvb_mult     ! GVB spin mult
 real(kind=8), allocatable :: noon(:)
 character(len=24) :: data_string
 character(len=240) :: buf, proname, mbgjf, gvb_nofch, mbout, pyname, outname
 ! mbgjf: minimal basis gjf

 if(ist /= 6) return
 i = (mult - 1)/2
 gvb_mult = DBLE(i*(i+1))

 i = INDEX(gjfname, '.gjf', back=.true.)
 proname = gjfname(1:i-1)
 mbgjf = gjfname(1:i-1)//'_mb.gjf'
 mbout = gjfname(1:i-1)//'_mb.out'
 hf_fch = gjfname(1:i-1)//'_proj_rem.fch'
 pyname = gjfname(1:i-1)//'_proj_rem.py'
 outname = gjfname(1:i-1)//'_proj_rem.out'

 call prt_automr_mb_gvb_gjf(gjfname, mbgjf, npair_wish, nskip_uno, localm, &
                            bgchg, .true., inherit)
 call submit_automr_job(mbgjf)
 call read_hf_and_gvb_e_from_automr_out(mbout, e, uhf_s2)
 call read_ndb_npair_nopen_from_automr_out(mbout, ndb, npair, nopen)

 write(6,'(/)',advance='no')
 if(mult == 1) then
  write(6,'(A,F18.8,1X,A)') 'GVB/STO-6G E(RHF) = ',e(1),'a.u., <S**2>=  0.000'
 end if
 write(6,'(A,F18.8,1X,A,F7.3)') 'GVB/STO-6G E(UHF) = ',e(2),'a.u., <S**2>=',uhf_s2
 write(6,'(A,F18.8,1X,A,F7.3)') 'GVB/STO-6G E(GVB) = ',e(3),'a.u., <S**2>=',gvb_mult

 call find_mb_gvb_nofch(proname, gvb_nofch)

 ! compress these minimal basis set files
 buf = 'tar -zcf '//TRIM(proname)//'_minb.tar.gz '//TRIM(proname)//&
       '_mb* --remove-files'
 write(6,'(A)') '$'//TRIM(buf)
 i = SYSTEM(TRIM(buf))

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine do_minimal_basis_gvb: failed to compress&
                & minimal basis related files.'
  stop
 end if

 ! decompress the gvb_nofch file
 buf = 'tar -zxf '//TRIM(proname)//'_minb.tar.gz '//TRIM(gvb_nofch)
 write(6,'(A)') '$'//TRIM(buf)
 i = SYSTEM(TRIM(buf))

 write(6,'(A)') 'GVB/STO-6G finished. Rotate MOs at target basis to resemble&
               & GVB/STO-6G orbitals...'

 call gen_fch_from_gjf(gjfname, hf_fch)
 call prt_orb_resemble_py_script(nproc, hf_fch, gvb_nofch, pyname)
 call submit_pyscf_job(pyname)

 ! Read GVB NOON of minimal basis set. The GVB/STO-6G may have many pairs, but
 ! we are only interested in moderate/strong-correlated pairs.
 call read_nbf_and_nif_from_fch(gvb_nofch, nbf, nif)
 allocate(noon(nif))
 call read_eigenvalues_from_fch(gvb_nofch, nif, 'a', noon)
 i = COUNT(noon(1:ndb+npair) < 1.98d0) ! 0.02~0.98
 deallocate(noon)

 ndb = ndb + npair - i ! update ndb and npair
 npair = i

 if(i == 0) then
  if(nopen == 0) then
   write(6,'(/,A)') 'ERROR in subroutine do_minimal_basis_gvb: GVB/STO-6G &
                    &results shows'
   write(6,'(A)') 'that this is a closed-shell singlet molecule.'
   stop
  else ! nopen > 0
   write(6,'(/,A)') 'Warning from subroutine do_minimal_basis_gvb: GVB/STO-6G s&
                   &hows that this'
   write(6,'(A)') 'molecule has npair=0. Calculation will be proceeded, but GVB&
                  & and/or CASSCF will'
   write(6,'(A)') 'be identical to ROHF.'
  end if
 end if

 call delete_file(TRIM(gvb_nofch))
 mo_rhf = .true. ! set to .True., actually ist=6, but mimicking ist=3

 call read_nbf_and_nif_from_fch(hf_fch, nbf, nif)
 open(newunit=fid,file='uno.out',status='replace')
 write(fid,'(A,I0,/,A,I0)') 'nbf=', nbf, 'nif=', nif
 write(fid,'(/,A,I0,/)') 'ndb=',ndb
 write(fid,'(A,I0,2(1X,I0))') 'idx=', ndb+1, ndb+2*npair+nopen+1, nopen
 close(fid)

 write(6,'(A)') 'Rotation done. HF_fch='//TRIM(hf_fch)
 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_minimal_basis_gvb at '//TRIM(data_string)
end subroutine do_minimal_basis_gvb

! print/create a GVB/STO-6G MOKIT automr .gjf file
subroutine prt_automr_mb_gvb_gjf(gjfname, mbgjf, npair, nskip_uno, localm, &
                                 bgchg, fcgvb, inherit)
 implicit none
 integer :: i, j, fid1, fid2
 integer, intent(in) :: npair, nskip_uno
 character(len=240) :: buf
 character(len=4), intent(in) :: localm
 character(len=240), intent(in) :: gjfname, mbgjf
 logical, intent(in) :: bgchg, fcgvb, inherit

 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(mbgjf),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:1) == '#') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 write(fid2,'(A)',advance='no') '#p GVB'

 if(inherit) then
  if(npair > 0) write(fid2,'(A,I0,A)',advance='no') '(',npair,')'
  i = INDEX(buf,'/')
  j = INDEX(buf(i+1:),' ')
  buf = '/STO-6G'//buf(i+j:)
  write(fid2,'(A,//,A)',advance='no') TRIM(buf),'mokit{LocalM='//TRIM(localm)
  if(nskip_uno > 0) then
   write(fid2,'(A,I0)',advance='no') ',skip_uno=', nskip_uno
  end if
  if(fcgvb) write(fid2,'(A,I0)',advance='no') ',FcGVB'
 else
  write(fid2,'(A,//,A)',advance='no') '/STO-6G','mokit{LocalM=PM'
 end if

 if(bgchg) write(fid2,'(A)',advance='no') ',charge'
 write(fid2,'(A,/)') ',GVB_conv=5d-4}'

 read(fid1,'(A)') buf
 read(fid1,'(A)') buf
 read(fid1,'(A)') buf

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1)
 close(fid2)
end subroutine prt_automr_mb_gvb_gjf

! search the _mb_*s.fch filename
subroutine find_mb_gvb_nofch(proname, gvb_nofch)
 implicit none
 integer :: i, fid, SYSTEM
 character(len=240), intent(in) :: proname
 character(len=240), intent(out) :: gvb_nofch

 i = SYSTEM('ls '//TRIM(proname)//'_mb_*s.fch >mb.txt')
 open(newunit=fid,file='mb.txt',status='old',position='rewind')
 read(fid,'(A)') gvb_nofch
 close(fid,status='delete')
end subroutine find_mb_gvb_nofch

! read RHF, UHF, GVB energies and UHF <S**2> from automr output file
subroutine read_hf_and_gvb_e_from_automr_out(outname, e, uhf_s2)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e(3), uhf_s2
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 e = 0d0 ; uhf_s2 = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6)=='E(RHF)' .or. buf(1:6)=='E(UHF)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_hf_and_gvb_e_from_automr_out: fail&
                   &ed to read HF energies'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 if(buf(1:6) == 'E(RHF)') then
  i = INDEX(buf, '=')
  read(buf(i+1:),*) e(1)
  read(fid,'(A)') buf ! read E(UHF) line
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) e(2)
 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) uhf_s2

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'E(GVB)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_hf_and_gvb_e_from_automr_out: fail&
                   &ed to read GVB energies'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) e(3)
 close(fid)
end subroutine read_hf_and_gvb_e_from_automr_out

! read ndb, npair, and nopen from a AutoMR output file
subroutine read_ndb_npair_nopen_from_automr_out(mbout, ndb, npair, nopen)
 implicit none
 integer :: i, fid
 integer, intent(out) :: ndb, npair, nopen
 character(len=240) :: buf
 character(len=240), intent(in) :: mbout

 ndb = 0; npair = 0; nopen = 0
 open(newunit=fid,file=TRIM(mbout),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:23) == 'Enter subroutine do_gvb') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_ndb_npair_nopen_from_automr_out: &
                   &no 'Enter subroutine do_gvb' found in"
  write(6,'(A)') 'file '//TRIM(mbout)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 close(fid)

 call read_int_in_buf(buf, 'doubly_occ=', ndb)
 call read_int_in_buf(buf, 'npair=', npair)
 call read_int_in_buf(buf, 'nopen=', nopen)
end subroutine read_ndb_npair_nopen_from_automr_out

! read an integer in buf according to the key, e.g.
!  read_int_in_buf('doubly_occ=10,xxx','doubly_occ=',k) -> k=10
subroutine read_int_in_buf(buf, key, k)
 implicit none
 integer :: i, j, m, n
 integer, intent(out) :: k
 character(len=*), intent(in) :: buf, key

 k = 0
 m = LEN(key)
 i = INDEX(buf, key)
 j = INDEX(buf(i+m:), ',')
 read(buf(i+m:i+m+j-2),*,iostat=n) k

 if(n /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_int_in_buf: failed to read the int&
                   &eger from'
  write(6,'(A)') "buf='"//buf//"'"
  stop
 end if
end subroutine read_int_in_buf

! call Gaussian to generate fch file from a given gjf file
subroutine gen_fch_from_gjf(gjfname, hf_fch)
 use util_wrapper, only: formchk
 use mr_keyword, only: gau_path
 implicit none
 integer :: i, j, mult, fid1, fid2, SYSTEM
 character(len=4) :: method
 character(len=240) :: buf, tmpchk, tmpgjf, tmpout
 character(len=240), intent(in) :: gjfname, hf_fch

 call read_mult_from_gjf(gjfname, mult)
 method = 'RHF'
 if(mult > 1) method = 'ROHF'

 call get_a_random_int(i)
 write(tmpchk,'(I0,A)') i,'.chk'
 write(tmpgjf,'(I0,A)') i,'.gjf'
 write(tmpout,'(I0,A)') i,'.log'
 open(newunit=fid1,file=TRIM(gjfname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(tmpgjf),status='replace')
 write(fid2,'(A)') '%chk='//TRIM(tmpchk)

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:1) == '#') exit
  if(buf(1:4) == '%chk') cycle
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 i = INDEX(buf,' ')
 j = INDEX(buf,'/')
 if(i > j) then
  write(6,'(A)') 'ERROR in subroutine gen_fch_from_gjf: wrong syntax in file '&
                  //TRIM(gjfname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 buf = buf(1:i)//TRIM(method)//TRIM(buf(j:))//' guess(only,save) nosymm 5D 7F&
                                              & int=nobasistransform'
 write(fid2,'(A)') TRIM(buf)

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid1)
 close(fid2)

 i = SYSTEM(TRIM(gau_path)//' '//TRIM(tmpgjf))
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine gen_fch_from_gjf: Gaussian 'ONLY'-type&
                  & job failed."
  write(6,'(A)') 'You can open file '//TRIM(tmpout)//' and check why.'
  stop
 end if

 call formchk(tmpchk, hf_fch)
 call delete_files(3, [tmpchk, tmpgjf, tmpout])
end subroutine gen_fch_from_gjf

! read spin multiplicity from Gaussian gjf file
subroutine read_mult_from_gjf(gjfname, mult)
 implicit none
 integer :: charge, nblank, fid
 integer, intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 mult = 1
 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 nblank = 0

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 read(fid,*) charge, mult
 close(fid)
end subroutine read_mult_from_gjf

! create/print a PySCF .py file for rotating MOs in fchname1 to resemble MOs
! in fchname2
subroutine prt_orb_resemble_py_script(nproc, fchname1, fchname2, pyname)
 use util_wrapper, only: bas_fch2py_wrap
 implicit none
 integer :: i, nbf1, nif1, nbf2, nif2, fid1, fid2, fid3, RENAME
 integer, intent(in) :: nproc
 character(len=240) :: buf, pyname1, pyname2
 character(len=240), intent(in) :: fchname1, fchname2
 character(len=240), intent(out) :: pyname

 call bas_fch2py_wrap(fchname1, .false.)
 call bas_fch2py_wrap(fchname2, .false.)
 call read_nbf_and_nif_from_fch(fchname1, nbf1, nif1)
 call read_nbf_and_nif_from_fch(fchname2, nbf2, nif2)

 i = INDEX(fchname1, '.fch', back=.true.)
 pyname = fchname1(1:i-1)//'.py0'
 pyname1 = fchname1(1:i-1)//'.py'

 i = INDEX(fchname2, '.fch', back=.true.)
 pyname2 = fchname2(1:i-1)//'.py'

 open(newunit=fid1,file=TRIM(pyname1),status='old',position='rewind')
 open(newunit=fid3,file=TRIM(pyname),status='replace')
 write(fid3,'(A)') 'from mokit.lib.mo_svd import orb_resemble'
 write(fid3,'(A)') 'from mokit.lib.py2fch import py2fch'
 write(fid3,'(A)') 'from mokit.lib.rwwfn import read_nbf_and_nif_from_fch'
 write(fid3,'(A)') 'import numpy as np'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:17) == 'from pyscf import') exit
  write(fid3,'(A)') TRIM(buf)
 end do ! for while
 buf= TRIM(buf)//', lib'
 write(fid3,'(A)') TRIM(buf)

 do while(.true.)
  read(fid1,'(A)') buf
  write(fid3,'(A)') TRIM(buf)
  if(buf(1:6) == 'mol.bu') exit
 end do ! for while

 close(fid1,status='delete')
 write(fid3,'(/,A)') '# copy this molecule'
 write(fid3,'(A)') 'mol2 = mol.copy()'

 open(newunit=fid2,file=TRIM(pyname2),status='old',position='rewind')
 do while(.true.)
  read(fid2,'(A)') buf
  if(buf(1:6) == 'mol.ba') exit
 end do ! for while

 write(fid3,'(A)') 'mol2.basis = {'
 do while(.true.)
  read(fid2,'(A)') buf
  write(fid3,'(A)') TRIM(buf)
  if(buf(1:5) == "''')}") exit
 end do ! for while

 close(fid2,status='delete')
 write(fid3,'(A)') 'mol2.build(parse_arg=False)'
 write(fid3,'(/,A)') "nbf1, nif1 = read_nbf_and_nif_from_fch('"//TRIM(fchname1)//"')"
 write(fid3,'(A)') "nbf2, nif2 = read_nbf_and_nif_from_fch('"//TRIM(fchname2)//"')"

 write(fid3,'(/,A,I0)') 'nproc = ', nproc
 write(fid3,'(A)') 'lib.num_threads(nproc)'

 write(fid3,'(/,A)') '# rotate MOs at target basis to resemble known orbitals'
 write(fid3,'(A)') "cross_S = gto.intor_cross('int1e_ovlp', mol, mol2)"
 write(fid3,'(A)') "mo1 = fch2py('"//TRIM(fchname1)//"', nbf1, nif1, 'a')"
 write(fid3,'(A)') "mo2 = fch2py('"//TRIM(fchname2)//"', nbf2, nif2, 'a')"
 write(fid3,'(A)') "mo3 = orb_resemble(nbf1, nif1, mo1, nbf2, nif2, mo2, cross_S)"
 write(fid3,'(A)') 'noon = np.zeros(nif1)'
 write(fid3,'(A)') "py2fch('"//TRIM(fchname1)//"', nbf1, nif1, mo3, 'a', noon, &
                   &False, False)"
 i = RENAME(TRIM(pyname), TRIM(pyname1))
 pyname = pyname1
end subroutine prt_orb_resemble_py_script

! find the number of active UNO pairs from a given .fch(k) file
! Note that UNO are in pairs naturally, so npair0 from occupied space
!  must be equal to that from unoccupied space
subroutine find_npair0_from_fch(fchname, nopen, npair0)
 use mr_keyword, only: on_thres
 implicit none
 integer :: i, fid, nif
 integer, intent(in) :: nopen
 integer, intent(out) :: npair0
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 real(kind=8), allocatable :: noon(:)

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'Alpha O') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine find_npair0_from_fch: keyword 'Alpha O'&
                 & not found in file "//TRIM(fchname)
  stop
 end if

 BACKSPACE(fid)
 read(fid,'(A49,2X,I10)') buf, nif
 allocate(noon(nif), source=0d0)
 read(fid,'(5(1X,ES15.8))') (noon(i),i=1,nif)
 close(fid)

 npair0 = 0
 do i = 1, nif, 1
  if(noon(i)>on_thres .and. noon(i)<(2d0-on_thres)) npair0 = npair0 + 1
 end do ! for i
 deallocate(noon)

 if(MOD(npair0-nopen,2) /= 0) then
  write(6,'(A)') 'ERROR in subroutine find_npair0_from_fch: npair0-nopen is not&
                & an even integer.'
  write(6,'(A)') "This is probably because UNO occupation numbers in 'Alpha O'&
                & are probably incorrect."
  stop
 end if

 npair0 = (npair0 - nopen)/2
end subroutine find_npair0_from_fch

! perform a RHF/STO-6G calculation, or a UHF/STO-6G calculation plus UNO
subroutine do_min_bas_hf_gau(proname)
 use mr_keyword, only: mem, nproc, gau_path
 use mol, only: charge, mult, natom, elem, coor
 use util_wrapper, only: formchk, fch_u2r_wrap
 implicit none
 integer :: i, fid
 character(len=240) :: gjfname, chkname, fchname, logname
 character(len=240), intent(in) :: proname

 gjfname = TRIM(proname)//'.gjf'
 chkname = TRIM(proname)//'.chk'
 fchname = TRIM(proname)//'.fch'
#ifdef _WIN32
 logname = TRIM(proname)//'.out'
#else
 logname = TRIM(proname)//'.log'
#endif

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p '
 if(mult == 1) then
  write(fid,'(A)',advance='no') 'R'
 else
  write(fid,'(A)',advance='no') 'U'
 end if
 write(fid,'(A)',advance='no') 'HF/STO-6G nosymm int=nobasistransform scf(xqc,m&
                               &axcycle=500)'
 if(mult > 1) write(fid,'(A)',advance='no') ' stable=opt'

 write(fid,'(//,A,//,I0,1X,I0)') 'Title', charge, mult
 do i = 1, natom, 1
  write(fid,'(A2,3(1X,F18.8))') elem(i), coor(:,i)
 end do ! for i

 ! We will not use ROHF here since it is sometimes not easy to be converged and
 ! `scf=qc` cannot be used for ROHF in Gaussian. Instead, we will use UNO.
 if(mult > 1) then
  write(fid,'(//,A)') '--Link1--'
  write(fid,'(A)') '%chk='//TRIM(chkname)
  write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
  write(fid,'(A,I0)') '%nprocshared=', nproc
  write(fid,'(A)') '#p chkbasis nosymm int=nobasistransform guess(read,only,sav&
                   &e,NaturalOrbitals) geom=allcheck'
 end if

 write(fid,'(/)',advance='no')
 close(fid)

 call submit_gau_job(gau_path, gjfname, .false.)
 call formchk(chkname)
 if(mult > 1) call fch_u2r_wrap(fchname, fchname)
 call delete_files(3, [gjfname, chkname, logname])
end subroutine do_min_bas_hf_gau

subroutine do_min_bas_hf_pyscf(proname)
 use mr_keyword, only: mem, nproc, dkh2_or_x2c
 use mol, only: charge, mult, natom, elem, coor
 implicit none
 integer :: i, fid
 character(len=240) :: pyname, outname, fchname, r_fch
 character(len=240), intent(in) :: proname

 pyname = TRIM(proname)//'.py'
 outname = TRIM(proname)//'.out'
 fchname = TRIM(proname)//'.fch'
 r_fch = TRIM(proname)//'_UNO.fch'
 open(newunit=fid,file=TRIM(pyname),status='replace')

 write(fid,'(A)') 'from pyscf import gto, scf, lib'
 write(fid,'(A)') 'from mokit.lib.py2fch_direct import fchk'
 if(mult > 1) then
  write(fid,'(A)') 'from mokit.lib.rwwfn import update_density_using_mo_in_fch'
  write(fid,'(A)') 'from mokit.lib.gaussian import uno'
  write(fid,'(A)') 'from os import rename'
 end if
 write(fid,'(/,A,I0,A)') 'lib.num_threads(', nproc, ')'

 write(fid,'(A)') 'mol = gto.M()'
 write(fid,'(A,I0,A)') '# ', natom, ' atom(s)'
 write(fid,'(A)') "mol.atom = '''"
 do i = 1, natom, 1
  write(fid,'(A2,1X,3(1X,F17.8))') elem(i), coor(:,i)
 end do ! for i
 write(fid,'(A)') "'''"
 write(fid,'(A)') "mol.basis = 'STO-6G'"
 write(fid,'(A)') '# Remember to check the charge and spin'
 write(fid,'(A,I0)') 'mol.charge = ', charge
 write(fid,'(A,I0)') 'mol.spin = ', mult-1
 write(fid,'(A,I0)') 'mol.cart = False'
 write(fid,'(A)') 'mol.verbose = 4'
 write(fid,'(A)') 'mol.build(parse_arg=False)'

 write(fid,'(/,A)',advance='no') 'mf = scf.'
 if(mult == 1) then
  write(fid,'(A)',advance='no') 'RHF'
 else
  write(fid,'(A)',advance='no') 'UHF'
 end if

 if(dkh2_or_x2c) then
  write(fid,'(A)') '(mol).x2c1e()'
 else
  write(fid,'(A)') '(mol)'
 end if

 write(fid,'(A,I0,A)') 'mf.max_memory = ',mem*1000,' # MB'
 write(fid,'(A)') 'mf.max_cycle = 200'
 write(fid,'(A)') 'old_e = mf.kernel()'

 ! If SCF is not converged, use the Newton method to continue
 write(fid,'(/,A)') 'if mf.converged is False:'
 write(fid,'(A)')   '  mf = mf.newton()'
 write(fid,'(A,/)') '  old_e = mf.kernel()'

 if(mult > 1) then ! UHF
  ! loop to check wave function stability
  write(fid,'(A)') 'new_e = old_e + 2e-5'
  write(fid,'(A)') 'i = 0'
  write(fid,'(A)') 'while(i < 10):'
  write(fid,'(A)') '  mo1 = mf.stability()[0]'
  write(fid,'(A)') '  dm1 = mf.make_rdm1(mo1, mf.mo_occ)'
  write(fid,'(A)') '  mf = mf.newton()'
  write(fid,'(A)') '  new_e = mf.kernel(dm0=dm1)'
  write(fid,'(A)') '  if(abs(new_e-old_e) < 1e-5):'
  write(fid,'(A)') '    break # cannot find lower solution'
  write(fid,'(A)') '  old_e = new_e'
  write(fid,'(A)') '  i += 1'
  write(fid,'(A)') 'if i == 10:'
  write(fid,'(A)') "  raise OSError('PySCF stable=opt failed after 10 attempts.')"
  write(fid,'(/,A)') '# save UHF MOs into .fch file'
  write(fid,'(A)') "uhf_fch = '"//TRIM(fchname)//"'"
  write(fid,'(A)') "r_fch = '"//TRIM(r_fch)//"'"
  write(fid,'(A)') 'fchk(mf, uhf_fch)'
  write(fid,'(A)') 'update_density_using_mo_in_fch(uhf_fch)'
  write(fid,'(A)') 'uno(uhf_fch)'
  write(fid,'(A)') 'rename(r_fch, uhf_fch)'
 else         ! RHF
  write(fid,'(A)') 'if mf.converged is False:'
  write(fid,'(A)') "  raise OSError('PySCF RHF job failed.')"
  write(fid,'(/,A)') '# save RHF MOs into .fch file'
  write(fid,'(A)') "fchk(mf, '"//TRIM(fchname)//"', density=True)"
 end if

 close(fid)
 call submit_pyscf_job(pyname)
 call delete_files(2, [pyname, outname])
end subroutine do_min_bas_hf_pyscf

