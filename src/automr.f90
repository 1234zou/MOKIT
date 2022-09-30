! written by jxzou at 20200420: black-box multireference calculations
! updated by jxzou at 20200511: framework of the program
! updated by jxzou at 20201213: %casscf+noiter -> %mrci for correct CASCI NOONs in ORCA
! updated by jxzou at 20201222: import grad for CASSCF force in PySCF; read CASSCF force for Molpro
! updated by jxzou at 20210111: add subroutine do_mrpt3
! updated by jxzou at 20210119: add OpenMolcas-QCMaquis, OpenMolcas-CheMPS2 interface

! The input file is just like Gaussian .gjf format. MOKIT keywords should be
!  specified in the Title Card line like 'mokit{}'.

program main
 use print_id, only: iout
 implicit none
 integer :: i, j
 character(len=240) :: fname = ' '

 i = iargc()
 if(i /= 1) then
  write(iout,'(/,A)') ' ERROR in subroutine automr: wrong command line argument!'
  write(iout,'(A)')   " Example 1 (in bash): automr a.gjf >& a.out &"
  write(iout,'(A)')   " Example 2 (in dash): automr a.gjf >a.out 2>&1 &"
  write(iout,'(A,/)') ' See help: automr -h'
  stop
 end if

 call getarg(1, fname)

 select case(TRIM(fname))
 case('-v', '-V', '--version')
  write(iout,'(A)') 'AutoMR 1.2.4 :: MOKIT, release date: 2022-Sep-30'
  stop
 case('-h','-help','--help')
  write(6,'(/,A)') "Usage: automr [gjfname] >& [outname]"
  write(6,'(A)')   "  Example 1 (in bash): automr a.gjf >& a.out &"
  write(6,'(A)')   "  Example 2 (in dash): automr a.gjf >a.out 2>&1 &"
  write(6,'(/,A)') 'Options:'
  write(6,'(A)')   '  -h, -help, --help: Print this message and exit.'
  write(6,'(A)')   '  -v, -V, --version: Print the version number of automr and exit.'
  write(6,'(/,A)') 'Methods(#p ...):'
  write(6,'(A)')   '  GVB, CASCI, CASSCF, DMRGCI, DMRGSCF, NEVPT2, NEVPT3,&
                   & CASPT2, CASPT2-K,'
  write(6,'(A)')   '  CASPT3, MRMP2, OVBMP2, SDSPT2, MRCISD, MRCISDT, MCPDFT,&
                   & FICMRCCSD,'
  write(6,'(A)')   '  MkMRCCSD, MkMRCCSD(T), BWMRCCSD, BWMRCCSD(T), BCCC2b, BCCC3b'
  write(6,'(/,A)') 'Frequently used keywords in MOKIT{}:'
  write(6,'(A)')   '      HF_prog=Gaussian/PySCF/ORCA/PSI4'
  write(6,'(A)')   '     GVB_prog=GAMESS/Gaussian'
  write(6,'(A)')   '  CASSCF_prog=PySCF/OpenMolcas/ORCA/Molpro/GAMESS/Gaussian/BDF/PSI4/Dalton'
  write(6,'(A)')   '   CASCI_prog=PySCF/OpenMolcas/ORCA/Molpro/GAMESS/Gaussian/BDF/PSI4/Dalton'
  write(6,'(A)')   '  NEVPT2_prog=PySCF/OpenMolcas/ORCA/Molpro/BDF'
  write(6,'(A)')   '  CASPT2_prog=OpenMolcas/Molpro/ORCA'
  write(6,'(A)')   '  MCPDFT_prog=OpenMolcas/GAMESS'
  write(6,'(A)')   '  MRCISD_prog=OpenMolcas/Molpro/ORCA/Gaussian/GAMESS/PSI4/Dalton'
  write(6,'(A)')   '      CtrType=1/2/3 for uc-/ic-/FIC- MRCISD'
  write(6,'(A,/)') '    MRCC_prog=ORCA/NWChem'
  stop
 end select

 i = index(fname, '.gjf', back=.true.)
 j = index(fname, '.fch', back=.true.)
 if(i/=0 .and. j/=0) then
  write(iout,'(/,A)') "ERROR in subroutine automr: both '.gjf' and '.fch' keys&
                     & detected in filename "//TRIM(fname)//'.'
  write(iout,'(A)') "Better to use a filename only with suffix '.gjf'."
  stop
 else if(i == 0) then
  write(iout,'(/,A)') "ERROR in subroutine automr: '.gjf' key not found in&
                     & filename "//TRIM(fname)
  stop
 end if

 call require_file_exist(fname)
 call automr(fname)
 stop
end program main

! automatically do multireference calculations in a block-box way
subroutine automr(fname)
 use print_id, only: iout
 use mr_keyword, only: gjfname, read_program_path, parse_keyword, check_kywd_compatible
 implicit none
 character(len=24) :: data_string
 character(len=240), intent(in) :: fname

 gjfname = fname
 ! read paths of various programs from environment variables
 call read_program_path()
 call parse_keyword()
 call check_kywd_compatible()

 call do_hf()         ! RHF and/or UHF
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
 call do_PES_scan()   ! PES scan

 call fdate(data_string)
 write(iout,'(/,A)') 'Normal termination of AutoMR at '//TRIM(data_string)
end subroutine automr

! generate PySCF input file .py from Gaussian .fch(k) file, and get paired LMOs
subroutine get_paired_LMO()
 use print_id, only: iout
 use mr_keyword, only: eist, mo_rhf, ist, hf_fch, bgchg, chgname, dkh2_or_x2c,&
  nskip_uno
 use mol, only: nbf, nif, ndb, nacte, nacto, nacta, nactb, npair, npair0, nopen,&
  lin_dep, chem_core, ecp_core
 implicit none
 integer :: i, system, RENAME
 real(kind=8) :: unpaired_e
 character(len=24) :: data_string = ' '
 character(len=240) :: buf, proname, pyname, chkname, outname, fchname

 if(eist == 1) return ! excited state calculation
 if(ist == 5) return ! no need for this subroutine
 write(iout,'(//,A)') 'Enter subroutine get_paired_LMO...'

 if(ist == 0) then
  write(iout,'(A)') 'ERROR in subroutine get_paired_LMO: ist=0. It should be&
                   & non-zero before this subroutine.'
  stop
 end if

 call calc_ncore() ! calculate the number of core orbitals from array core_orb
 write(iout,'(3(A,I0))') 'chem_core=', chem_core, ', ecp_core=', ecp_core, &
                         ', Nskip_UNO=', nskip_uno

 i = index(hf_fch, '.fch', back=.true.)
 proname = hf_fch(1:i-1)

 if(mo_rhf) then
  if(ist == 6) then
   write(iout,'(A)') 'One set of MOs: GVB/STO-6G -> MO projection -> Rotate has&
                     & been invoked.'
  else
   write(iout,'(A)') 'One set of MOs: invoke RHF virtual MO projection ->&
                    & localization -> paring.'
   chkname = hf_fch(1:i-1)//'_proj.chk' ! this is PySCF chk file, not Gaussian
   i = system('bas_fch2py '//TRIM(hf_fch))
   pyname = TRIM(proname)//'_proj_loc_pair.py'
   i = RENAME(TRIM(proname)//'.py', TRIM(pyname))
   if(dkh2_or_x2c) call add_X2C_into_py(pyname)

   call prt_rhf_proj_script_into_py(pyname)
   call prt_auto_pair_script_into_py(pyname)
   write(buf,'(A)') 'python '//TRIM(pyname)//' >'//TRIM(proname)//&
                    '_proj_loc_pair.out 2>&1'
   write(iout,'(A)') '$'//TRIM(buf)
   i = system(TRIM(buf))
   call delete_file(chkname)
  end if

 else
  if(ist == 1) then
   write(iout,'(A)') 'Two sets of MOs, ist=1, invoke UNO associated rotation.'
  else if(ist == 2) then
   write(iout,'(A)') 'Two sets of MOs, ist=2, invoke UNO generation.'
  end if

  fchname = hf_fch(1:i-1)//'_uno.fch'
  i = system('bas_fch2py '//TRIM(hf_fch))
  if(ist == 1) then
   pyname = TRIM(proname)//'_uno_asrot.py'
   outname = TRIM(proname)//'_uno_asrot.out'
  else
   pyname = TRIM(proname)//'_uno.py'
   outname = TRIM(proname)//'_uno.out'
  end if

  i = RENAME(TRIM(proname)//'.py', TRIM(pyname))
  if(dkh2_or_x2c) call add_X2C_into_py(pyname)
  call prt_uno_script_into_py(pyname)

  if(ist == 1) call prt_assoc_rot_script_into_py(pyname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(pyname))
  write(buf,'(A)') 'python '//TRIM(pyname)//' >'//TRIM(outname)//" 2>&1"
  write(iout,'(A)') '$'//TRIM(buf)

  i = system(TRIM(buf))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine get_paired_LMO: PySCF job fails.'
   write(iout,'(A)') 'Please check file '//TRIM(outname)
   stop
  end if
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
 write(iout,'(A)') 'Leave subroutine get_paired_LMO at '//TRIM(data_string)
end subroutine get_paired_LMO

! print RHF virtual MOs projection scheme into a given .py file
subroutine prt_rhf_proj_script_into_py(pyname)
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, tencycle, hf_fch, dkh2_or_x2c
 use mol, only: natom, nuc, elem
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, allocatable :: ntimes(:)
 character(len=240) :: buf, pyname1, proj_fch, chkname
 character(len=240), intent(in) :: pyname

 pyname1 = TRIM(pyname)//'.tmp'
 buf = ' '
 i = index(hf_fch, '.fch')
 proj_fch = hf_fch(1:i-1)//'_proj.fch'
 chkname = hf_fch(1:i-1)//'_proj.chk' ! for PySCF STO-6G, not Gaussian .chk

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine prt_rhf_proj_script_into_py: end-of-file detected.'
  write(iout,'(A)') 'File may be incomplete: '//TRIM(pyname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A)') 'from pyscf import lib'
 write(fid2,'(A)') 'import numpy as np'
 write(fid2,'(A)') 'from mo_svd import mo_svd'
 write(fid2,'(A)') 'from py2fch import py2fch'
 write(fid2,'(A,/)') 'from shutil import copyfile'
 write(fid2,'(A,I0,A1,/)') 'lib.num_threads(',nproc,')'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:9) == 'mf.kernel') exit
  write(fid2,'(A)') TRIM(buf)
 end do
 write(fid2,'(A,I0,A)') 'mf.max_memory = ', mem*1000, ' # MB'
 write(fid2,'(A)') 'mf.kernel()'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:3) == '#dm') exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(tencycle) then
  write(fid2,'(A)') TRIM(buf(2:))
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf(2:))
  write(fid2,'(A)') "mf.chkfile = '"//TRIM(chkname)//"'"
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf(2:))
 else ! keep 10 cycle annotated
  write(fid2,'(A)') TRIM(buf)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
 end if

 close(fid1,status='delete')
 write(fid2,'(/,A)') '# copy this molecule at STO-6G'
 write(fid2,'(A)') 'mol2 = mol.copy()'

 if(ANY(nuc > 54)) then ! atoms > Xe
  allocate(ntimes(natom))
  call calc_ntimes(natom, elem, ntimes)
  write(fid2,'(A)') "mol2.basis = {'default':'STO-6G',"
  do i = 1, natom, 1
   if(nuc(i) > 54) write(fid2,'(A,I0,A)') "'"//TRIM(elem(i)),ntimes(i),"':'def2-SVP',"
  end do ! for i
  write(fid2,'(A)') '}'
  write(fid2,'(A)') 'mol2.ecp = mol.ecp.copy()'
  deallocate(ntimes)
 else
  write(fid2,'(A)') "mol2.basis = 'STO-6G'"
 end if

 write(fid2,'(A)') 'mol2.build()'
 if(dkh2_or_x2c) then
  write(fid2,'(A)') 'mf2 = scf.RHF(mol2).x2c1e()'
 else
  write(fid2,'(A)') 'mf2 = scf.RHF(mol2)'
 end if
 write(fid2,'(A)') 'mf2.max_cycle = 150'
 write(fid2,'(A)') "dm = mf2.from_chk('"//TRIM(chkname)//"')"
 write(fid2,'(A)') 'mf2.kernel(dm)'
 write(fid2,'(A)') 'nbf2 = mf2.mo_coeff.shape[0]'
 write(fid2,'(A)') 'nif2 = mf2.mo_coeff.shape[1]'
 write(fid2,'(/,A)') '# project virtual MOs onto those of RHF/STO-6G'
 write(fid2,'(A)') "cross_S = gto.intor_cross('int1e_ovlp', mol, mol2)"
 write(fid2,'(A)') 'idx = np.sum(mf.mo_occ > 0)'
 write(fid2,'(A)') 'npair = np.sum(mf2.mo_occ==0)'
 write(fid2,'(A)') 'nmo1 = nif - idx'
 write(fid2,'(A)') 'nmo2 = nif2 - idx'
 write(fid2,'(A)') 'coeff1 = mf.mo_coeff[:,range(idx,nif)]'
 write(fid2,'(A)') 'coeff2 = mf2.mo_coeff[:,range(idx,nif2)]'
 write(fid2,'(A)') 'mo_svd(nbf, nmo1, nbf2, nmo2, coeff1, coeff2, cross_S, False)'
 write(fid2,'(A)') 'mf.mo_coeff[:,range(idx,nif)] = coeff1.copy()'
 write(fid2,'(A)') '# project done'

 write(fid2,'(/,A)') '# save projected MOs into a new .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(hf_fch)//"', '"//TRIM(proj_fch)//"')"
 write(fid2,'(A)') "py2fch('"//TRIM(proj_fch)//"',nbf,nif,mf.mo_coeff,'a',mf.mo_occ,False)"
 write(fid2,'(A)') '# save done'

 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 return
end subroutine prt_rhf_proj_script_into_py

! print localization and automatically pairing information into a given .py file
subroutine prt_auto_pair_script_into_py(pyname)
 use print_id, only: iout
 use mr_keyword, only: localm, hf_fch
 use mol, only: chem_core, ecp_core
 implicit none
 integer :: i, ncore, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1, loc_fch
 character(len=240), intent(in) :: pyname

 ncore = chem_core - ecp_core
 pyname1 = TRIM(pyname)//'.tmp'
 buf = ' '
 i = index(hf_fch, '.fch')
 loc_fch = hf_fch(1:i-1)//'_proj_loc_pair.fch'

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine prt_auto_pair_script_into_py: end-of-file detected.'
  write(iout,'(A)') 'File may be incomplete: '//TRIM(pyname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A)') 'from rwwfn import get_1e_exp_and_sort_pair as sort_pair'
 if(localm == 'pm') then ! Pipek-Mezey localization
  write(fid2,'(A)') 'from lo import pm'
 else ! Boys localization
  write(fid2,'(A)') 'from lo import boys'
 end if
 write(fid2,'(A)') 'from pyscf.lo.boys import dipole_integral'
 write(fid2,'(A,/)') 'from auto_pair import pair_by_tdm'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')

 write(fid2,'(/,A)') '# localize the projected MOs of larger basis'
 write(fid2,'(A,I0)') 'ncore = ', ncore
 write(fid2,'(A)') 'idx2 = np.sum(mf.mo_occ==2)'
 write(fid2,'(A)') 'idx1 = min(npair, idx2-ncore)'
 write(fid2,'(A)') 'occ_idx = range(ncore,idx2)'
 write(fid2,'(A)') 'vir_idx = range(idx2,nif2)'
 write(fid2,'(A)') 'idx3 = npair # backup'
 write(fid2,'(A)') 'npair = idx1 # update'

 if(localm == 'pm') then ! Pipek-Mezey localization
  write(fid2,'(A)') "occ_loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],&
                   &mol._bas[:,3],mol.cart,nbf,idx2-ncore,mf.mo_coeff[:,occ_idx],&
                   &S,'mulliken')"
  write(fid2,'(A)') "vir_loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],&
                   &mol._bas[:,3],mol.cart,nbf,nif2-idx2,mf.mo_coeff[:,vir_idx],&
                   &S,'mulliken')"
 else ! Boys localization
  write(fid2,'(A)') 'mo_dipole = dipole_integral(mol, mf.mo_coeff[:,occ_idx])'
  write(fid2,'(A)') 'occ_loc_orb = boys(nbf, idx2-ncore, mf.mo_coeff[:,occ_idx], mo_dipole)'
  write(fid2,'(A)') 'mo_dipole = dipole_integral(mol, mf.mo_coeff[:,vir_idx])'
  write(fid2,'(A)') 'vir_loc_orb = boys(nbf, nif2-idx2, mf.mo_coeff[:,vir_idx], mo_dipole)'
 end if

 write(fid2,'(A)') 'mf.mo_coeff[:,occ_idx] = occ_loc_orb.copy()'
 write(fid2,'(A)') 'mf.mo_coeff[:,vir_idx] = vir_loc_orb.copy()'
 write(fid2,'(A)') '# localization done'

 write(fid2,'(/,A)') '# pair the active orbitals'
 write(fid2,'(A)') 'mo_dipole = dipole_integral(mol, mf.mo_coeff)'
 write(fid2,'(A)') 'nopen = np.sum(mf.mo_occ==1)'
 write(fid2,'(A)') 'nalpha = np.sum(mf.mo_occ > 0)'
 write(fid2,'(A)') 'nvir_lmo = idx3'
 write(fid2,'(A)') 'alpha_coeff = pair_by_tdm(ncore, npair, nopen, nalpha, nvir_lmo,&
                   &nbf, nif, mf.mo_coeff, mo_dipole)'
 write(fid2,'(A)') 'mf.mo_coeff = alpha_coeff.copy()'
 write(fid2,'(A)') '# pair done'

 write(fid2,'(/,A)') '# save the paired LMO into .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(hf_fch)//"', '"//TRIM(loc_fch)//"')"
 write(fid2,'(A)') "py2fch('"//TRIM(loc_fch)//"',nbf,nif,mf.mo_coeff,'a',mf.mo_occ,False)"
 write(fid2,'(A)') "sort_pair('"//TRIM(loc_fch)//"','"//TRIM(hf_fch)//"',npair)"
 write(fid2,'(A)') '# save done'

 write(fid2,'(/,A)') "f = open('uno.out', 'w+')"
 write(fid2,'(A)') "f.write('nbf=%i\n' %nbf)"
 write(fid2,'(A)') "f.write('nif=%i\n\n' %nif)"
 write(fid2,'(A)') 'idx1 = idx2 - npair'
 write(fid2,'(A)') "f.write('ndb=%i\n\n' %idx1)"
 write(fid2,'(A)') "f.write('idx=%i %i %i' %(idx1+1,nalpha+npair+1,nopen))"
 write(fid2,'(A)') 'f.close()'
 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 return
end subroutine prt_auto_pair_script_into_py

! print UNO script into a given .py file
subroutine prt_uno_script_into_py(pyname)
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, hf_fch, tencycle, uno_thres
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1, uno_fch
 character(len=240), intent(in) :: pyname

 pyname1 = TRIM(pyname)//'.tmp'
 buf = ' '

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine prt_uno_script_into_py: end-of-file detected.'
  write(iout,'(A)') 'File may be incomplete: '//TRIM(pyname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A)') 'from pyscf import lib'
 write(fid2,'(A)') 'import numpy as np'
 write(fid2,'(A)') 'import os'
 write(fid2,'(A)') 'from time import sleep'
 write(fid2,'(A)') 'from py2fch import py2fch'
 write(fid2,'(A)') 'from uno import uno'
 write(fid2,'(A,/)') 'from construct_vir import construct_vir'
 write(fid2,'(A,I0,A1,/)') 'lib.num_threads(',nproc,')'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:9) == 'mf.kernel') exit
  write(fid2,'(A)') TRIM(buf)
 end do
 write(fid2,'(A,I0,A)') 'mf.max_memory = ', mem*1000, ' # MB'
 write(fid2,'(A)') 'mf.kernel()'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:3) == '#dm') exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(tencycle) then
  write(fid2,'(A)') TRIM(buf(2:))
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf(2:))
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf(2:))
 else ! keep 10 cycle annotated
  write(fid2,'(A)') TRIM(buf)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
 end if
 close(fid1,status='delete')
 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 i = index(hf_fch, '.fch', back=.true.)
 uno_fch = hf_fch(1:i-1)//'_uno.fch'

 open(newunit=fid1,file=TRIM(pyname),status='old',position='append')
 write(fid1,'(/,A)') '# transform UHF canonical orbitals to UNO'
 write(fid1,'(A)') 'na = np.sum(mf.mo_occ[0]==1)'
 write(fid1,'(A)') 'nb = np.sum(mf.mo_occ[1]==1)'
 write(fid1,'(A,E12.5,A)') 'idx, noon, alpha_coeff = uno(nbf,nif,na,nb,mf.mo_coeff[0],&
                           &mf.mo_coeff[1],S,',uno_thres,')'
 write(fid1,'(A)') 'alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)'
 write(fid1,'(A)') 'mf.mo_coeff = (alpha_coeff, beta_coeff)'
 write(fid1,'(A)') '# done transform'

 write(fid1,'(/,A)') '# save the UNO into .fch file'
 write(fid1,'(A)') "os.system('fch_u2r "//TRIM(hf_fch)//' '//TRIM(uno_fch)//"')"
 write(fid1,'(A)') 'sleep(1) # in some node, py2fch begins when fch_u2r unfinished'
 write(fid1,'(A)') "py2fch('"//TRIM(uno_fch)//"',nbf,nif,mf.mo_coeff[0],'a',noon,True)"
 write(fid1,'(A)') '# save done'
 close(fid1)
 return
end subroutine prt_uno_script_into_py

! print associated rotation into a given .py file
subroutine prt_assoc_rot_script_into_py(pyname)
 use mol, only: chem_core, ecp_core
! use mol, only: natom, elem, nuc, chem_core, ecp_core
 use mr_keyword, only : localm, hf_fch, npair_wish, nskip_uno
 implicit none
 integer :: i, ncore, fid1, fid2, RENAME
! integer, allocatable :: ntimes(:)
 character(len=240) :: buf, pyname1, uno_fch, assoc_fch
 character(len=240), intent(in) :: pyname

 ncore = chem_core - ecp_core
 pyname1 = TRIM(pyname)//'.tmp'
 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine prt_assoc_rot_script_into_py: end-of-file detected.'
  write(6,'(A)') 'File may be incomplete: '//TRIM(pyname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 if(localm == 'pm') then
  write(fid2,'(A)') 'from lo import pm, boys'
 else
  write(fid2,'(A)') 'from lo import boys'
 end if
 write(fid2,'(A)') 'from pyscf.lo.boys import dipole_integral'
 write(fid2,'(A)') 'from auto_pair import pair_by_tdm'
 write(fid2,'(A)') 'from assoc_rot import assoc_rot'
 write(fid2,'(A)') 'from mo_svd import proj_occ_get_act_vir'
 write(fid2,'(A)') 'from rwwfn import get_1e_exp_and_sort_pair as sort_pair'
 write(fid2,'(A,/)') 'from shutil import copyfile'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')
 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 i = index(hf_fch, '.fch', back=.true.)
 uno_fch = hf_fch(1:i-1)//'_uno.fch'
 assoc_fch = hf_fch(1:i-1)//'_uno_asrot.fch'

 open(newunit=fid1,file=TRIM(pyname),status='old',position='append')
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
  write(fid1,'(A)') "  occ_loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],&
                    &mol._bas[:,3],mol.cart,nbf,"
  write(fid1,'(19X,A)') "npair-i,mf.mo_coeff[0][:,occ_idx],S,'mulliken')"
 else ! Boys localization
  write(fid1,'(A)') '  mo_dipole = dipole_integral(mol, mf.mo_coeff[0][:,occ_idx])'
  write(fid1,'(A)') '  occ_loc_orb = boys(nbf, npair-i, mf.mo_coeff[0][:,&
                     & occ_idx], mo_dipole)'
 end if

 write(fid1,'(A)') '  vir_loc_orb = assoc_rot(nbf, npair-i, mf.mo_coeff[0][:,occ_idx],&
                    & occ_loc_orb,'
 write(fid1,'(26X,A)') 'mf.mo_coeff[0][:,vir_idx])'
 write(fid1,'(A)') '  mf.mo_coeff[0][:,occ_idx] = occ_loc_orb.copy()'
 write(fid1,'(A)') '  mf.mo_coeff[0][:,vir_idx] = vir_loc_orb.copy()'
 write(fid1,'(A)') '# localization done'

! write(fid1,'(/,A)') '# using projection to supplement pairs (if needed)'
! write(fid1,'(A)') 'mol2 = mol.copy()'
! if(ANY(nuc > 54)) then ! atoms > Xe
!  allocate(ntimes(natom))
!  call calc_ntimes(natom, elem, ntimes)
!  write(fid1,'(A)') "mol2.basis = {'default':'STO-6G',"
!  do i = 1, natom, 1
!   if(nuc(i) > 54) write(fid2,'(A,I0,A)') "'"//TRIM(elem(i)),ntimes(i),"':'def2-SVP',"
!  end do ! for i
!  write(fid1,'(A)') '}'
!  write(fid1,'(A)') 'mol2.ecp = mol.ecp.copy()'
!  deallocate(ntimes)
! else
!  write(fid1,'(A)') "mol2.basis = 'STO-6G'"
! end if
! write(fid1,'(A)') 'mol2.build()'
! write(fid1,'(A)') 'nshl = mol2.nbas'
! write(fid1,'(A)') 'ang = mol2._bas[:,1].copy()'
! write(fid1,'(A)') 'if(mol2.cart == True):'
! write(fid1,'(A)') '  for i in range(0,nshl):'
! write(fid1,'(A)') '    ang[i] = (ang[i]+1)*(ang[i]+2)/2'
! write(fid1,'(A)') 'else:'
! write(fid1,'(A)') '  for i in range(0,nshl):'
! write(fid1,'(A)') '    ang[i] = 2*ang[i] + 1'
! write(fid1,'(A)') 'nbf2 = np.dot(ang, mol2._bas[:,3])'
! write(fid1,'(A)') 'na_np = idx[1] - 1'
! write(fid1,'(A)') 'nalpha = na_np - npair'
! write(fid1,'(A)') 'npair0 = nbf2 - nalpha # number of virtual MOs in STO-6G'
! write(fid1,'(A)') 'ncr = npair0 - npair'
! write(fid1,'(A)') 'if(ncr > 0):'
! write(fid1,'(A)') "  S2 = mol2.intor_symmetric('int1e_ovlp')"
! write(fid1,'(A)') "  cross_S = gto.intor_cross('int1e_ovlp', mol, mol2)"
! write(fid1,'(A)') '  alpha_coeff = proj_occ_get_act_vir(nbf, nif, nbf2, na_np,&
!                   & S2, cross_S, mf.mo_coeff[0])'
! write(fid1,'(A)') '  mf.mo_coeff = (alpha_coeff, beta_coeff)'
! write(fid1,'(A)') '  nopen = idx[2]'
! write(fid1,'(A)') '  ndb = nalpha - nopen - npair0'
! write(fid1,'(A,I0)') '  ncore = ', ncore
! write(fid1,'(A)') '  nloc = ndb + ncr - ncore'
! write(fid1,'(A)') "  f = open('uno.out', 'w+')"
! write(fid1,'(A)') "  f.write('nbf=%i\n' %nbf)"
! write(fid1,'(A)') "  f.write('nif=%i\n\n' %nif)"
! write(fid1,'(A)') "  f.write('ndb=%i\n\n' %ndb)"
! write(fid1,'(A)') "  f.write('idx=%i %i %i' %(ndb+1,nalpha+npair0+1,nopen))"
! write(fid1,'(A)') '  f.close()'
! write(fid1,'(A)') '  occ_idx = range(ncore, ndb+ncr)'
! write(fid1,'(A)') '  vir_idx = range(na_np, na_np+ncr)'
! write(fid1,'(A)') '  occ_idx1 = range(0, nloc)'
! write(fid1,'(A)') '  vir_idx1 = range(nloc, nloc+ncr)'
! write(fid1,'(A)') '  coeff = np.zeros((nbf, nloc+ncr))'
! ! These pairs are sigma bonds and lone electrons, just use boys localization
! write(fid1,'(A)') '  mo_dipole = dipole_integral(mol, mf.mo_coeff[0][:,occ_idx])'
! write(fid1,'(A)') '  coeff[:,occ_idx1] = boys(nbf, nloc, mf.mo_coeff[0][:,occ_idx], mo_dipole)'
! write(fid1,'(A)') '  mo_dipole = dipole_integral(mol, mf.mo_coeff[0][:,vir_idx])'
! write(fid1,'(A)') '  coeff[:,vir_idx1] = boys(nbf, ncr, mf.mo_coeff[0][:,vir_idx], mo_dipole)'
! write(fid1,'(A)') '  # pair the supplement active orbitals'
! write(fid1,'(A)') '  mo_dipole = dipole_integral(mol, coeff)'
! write(fid1,'(A)') '  alpha_coeff = pair_by_tdm(0,ncr,0,nloc,ncr,nbf,nloc+ncr,coeff,mo_dipole)'
! write(fid1,'(A)') '  mf.mo_coeff[0][:,occ_idx] = alpha_coeff[:,occ_idx1].copy()'
! write(fid1,'(A)') '  mf.mo_coeff[0][:,vir_idx] = alpha_coeff[:,vir_idx1].copy()'

 write(fid1,'(/,A)') '# save associated rotation MOs into .fch(k) file'
 write(fid1,'(A)') "copyfile('"//TRIM(uno_fch)//"', '"//TRIM(assoc_fch)//"')"
 write(fid1,'(A)') 'noon = np.zeros(nif)'
 write(fid1,'(A)') "py2fch('"//TRIM(assoc_fch)//"',nbf,nif,mf.mo_coeff[0],'a',noon,False)"
 write(fid1,'(A)') "sort_pair('"//TRIM(assoc_fch)//"','"//TRIM(uno_fch)//"',npair)"
 close(fid1)
end subroutine prt_assoc_rot_script_into_py

! calculate/determine the number of core orbitals
subroutine calc_ncore()
 use mr_keyword, only: hf_fch
 use mol, only: natom, nuc, chem_core, ecp_core
 use fch_content, only: core_orb, RNFroz
 implicit none
 integer :: i, fid
 character(len=240) :: buf

 chem_core = 0 ! initialization
 ecp_core = 0
 buf = ' '

 do i = 1, natom, 1
  chem_core = chem_core + core_orb(nuc(i))
 end do ! for i

 open(newunit=fid,file=TRIM(hf_fch),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:10) == 'ECP-RNFroz') exit
  if(buf(1:7) == 'Alpha O') then
   close(fid)
   return
  end if
 end do ! for while

 if(allocated(RNFroz)) deallocate(RNFroz)
 allocate(RNFroz(natom), source=0d0)
 read(fid,'(5(1X,ES15.8))') (RNFroz(i), i=1,natom)
 close(fid)

 ecp_core = INT(SUM(RNFroz))
end subroutine calc_ncore

! perform GVB/STO-6G, only valid for ist=6
subroutine do_minimal_basis_gvb()
 use mol, only: mult, nbf, nif, nopen, ndb, npair
 use mr_keyword, only: nproc, ist, npair_wish, gjfname, localm, hf_fch, mo_rhf,&
  nskip_uno, bgchg, inherit
 implicit none
 integer :: i, fid, system
 real(kind=8) :: e(3), uhf_s2 ! RHF/UHF/GVB energies and UHF spin mult
 real(kind=8) :: gvb_mult     ! GVB spin mult
 real(kind=8), allocatable :: noon(:)
 character(len=24) :: data_string
 character(len=240) :: buf, proname, mbgjf, gvb_nofch, mbout, pyname, outname
 ! mbgjf: minimal basis gjf

 if(ist /= 6) return
 i = (mult - 1)/2
 gvb_mult = DBLE(i*(i+1))

 i = index(gjfname, '.gjf', back=.true.)
 proname = gjfname(1:i-1)
 mbgjf = gjfname(1:i-1)//'_mb.gjf'
 mbout = gjfname(1:i-1)//'_mb.out'
 hf_fch = gjfname(1:i-1)//'_proj_rem.fch'
 pyname = gjfname(1:i-1)//'_proj_rem.py'
 outname = gjfname(1:i-1)//'_proj_rem.out'

 call prt_automr_mb_gvb_gjf(gjfname, mbgjf, npair_wish, nskip_uno, localm, &
                            bgchg, inherit)
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
 i = system(TRIM(buf))

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine do_minimal_basis_gvb: failed to compress&
                & minimal basis related files.'
  stop
 end if

 ! decompress the gvb_nofch file
 buf = 'tar -zxf '//TRIM(proname)//'_minb.tar.gz '//TRIM(gvb_nofch)
 write(6,'(A)') '$'//TRIM(buf)
 i = system(TRIM(buf))

 write(6,'(A)') 'GVB/STO-6G finished. Rotate MOs at target basis to resemble&
               & GVB/STO-6G orbitals...'

 call gen_fch_from_gjf(gjfname, hf_fch)
 call prt_orb_resemble_py_script(nproc, hf_fch, gvb_nofch, pyname)

 buf = 'python '//TRIM(pyname)//' >'//TRIM(outname)//" 2>&1"
 write(6,'(A)') '$'//TRIM(buf)
 i = system(TRIM(buf))
 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine do_minimal_basis_gvb: PySCF job failed.'
  write(6,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
  stop
 end if

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
                    &results shows that'
   write(6,'(A)') 'this is a closed-shell singlet molecule.'
   stop
  else ! nopen > 0
   write(6,'(/,A)') 'Warning from subroutine do_minimal_basis_gvb: GVB/STO-6G&
                   & shows that this molecule'
   write(6,'(A)') 'has npair=0. Calculation will be proceeded, but GVB and/or&
                 & CASSCF will be identical to ROHF.'
  end if
 end if

 call delete_file(gvb_nofch)
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
                                 bgchg, inherit)
 implicit none
 integer :: i, j, fid1, fid2
 integer, intent(in) :: npair, nskip_uno
 character(len=240) :: buf
 character(len=4), intent(in) :: localm
 character(len=240), intent(in) :: gjfname, mbgjf
 logical, intent(in) :: bgchg, inherit

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
  i = index(buf,'/')
  j = index(buf(i+1:),' ')
  buf = '/STO-6G'//buf(i+j:)
  write(fid2,'(A,//,A)',advance='no') TRIM(buf),'mokit{LocalM='//TRIM(localm)
  if(nskip_uno > 0) then
   write(fid2,'(A,I0)',advance='no') ',skip_uno=', nskip_uno
  end if
 else
  write(fid2,'(A,//,A)',advance='no') '/STO-6G','mokit{LocalM=PM'
 end if

 if(bgchg) write(fid2,'(A)',advance='no') ',charge'
 write(fid2,'(A,/)') ',GVB_conv=5D-4}'

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
 integer :: i, fid, system
 character(len=240), intent(in) :: proname
 character(len=240), intent(out) :: gvb_nofch

 i = system('ls '//TRIM(proname)//'_mb_*s.fch >mb.txt')
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
  write(6,'(A)') 'ERROR in subroutine read_hf_and_gvb_e_from_automr_out: failed&
                & to read HF energies'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 if(buf(1:6) == 'E(RHF)') then
  i = index(buf, '=')
  read(buf(i+1:),*) e(1)
  read(fid,'(A)') buf ! read E(UHF) line
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) e(2)
 i = index(buf, '=', back=.true.)
 read(buf(i+1:),*) uhf_s2

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:6) == 'E(GVB)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine read_hf_and_gvb_e_from_automr_out: failed&
                & to read GVB energies'
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
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
 read(fid,'(A)') buf
 close(fid)

 i = index(buf, '=')
 read(buf(i+1:),*) ndb
 buf(i:i) = ' '

 i = index(buf, '=')
 read(buf(i+1:),*) npair
 buf(i:i) = ' '

 i = index(buf, '=')
 read(buf(i+1:),*) nopen
end subroutine read_ndb_npair_nopen_from_automr_out

! call Gaussian to generate fch file from a given gjf file
subroutine gen_fch_from_gjf(gjfname, hf_fch)
 use util_wrapper, only: formchk
 use mr_keyword, only: gau_path
 implicit none
 integer :: i, j, mult, fid1, fid2, system
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

 i = index(buf,' ')
 j = index(buf,'/')
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

 i = system(TRIM(gau_path)//' '//TRIM(tmpgjf))
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
 integer :: i, charge, nblank, fid
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

 call bas_fch2py_wrap(fchname1)
 call bas_fch2py_wrap(fchname2)
 call read_nbf_and_nif_from_fch(fchname1, nbf1, nif1)
 call read_nbf_and_nif_from_fch(fchname2, nbf2, nif2)

 i = index(fchname1, '.fch', back=.true.)
 pyname = fchname1(1:i-1)//'.py0'
 pyname1 = fchname1(1:i-1)//'.py'

 i = index(fchname2, '.fch', back=.true.)
 pyname2 = fchname2(1:i-1)//'.py'

 open(newunit=fid1,file=TRIM(pyname1),status='old',position='rewind')
 open(newunit=fid3,file=TRIM(pyname),status='replace')
 write(fid3,'(A)') 'from mo_svd import orb_resemble'
 write(fid3,'(A)') 'from py2fch import py2fch'
 write(fid3,'(A)') 'from rwwfn import read_nbf_and_nif_from_fch'
 write(fid3,'(A)') 'import numpy as np'
 write(fid3,'(A)') 'from pyscf import lib'

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
 write(fid3,'(A)') 'mol2.build()'
 write(fid3,'(/,A)') "nbf1, nif1 = read_nbf_and_nif_from_fch('"//TRIM(fchname1)//"')"
 write(fid3,'(A)') "nbf2, nif2 = read_nbf_and_nif_from_fch('"//TRIM(fchname2)//"')"

 write(fid3,'(/,A,I0,A)') 'lib.num_threads(',nproc,')'
 write(fid3,'(A)') '# rotate MOs at target basis to resemble known orbitals'
 write(fid3,'(A)') "cross_S = gto.intor_cross('int1e_ovlp', mol, mol2)"
 write(fid3,'(A)') "mo1 = fch2py('"//TRIM(fchname1)//"', nbf1, nif1, 'a')"
 write(fid3,'(A)') "mo2 = fch2py('"//TRIM(fchname2)//"', nbf2, nif2, 'a')"
 write(fid3,'(A)') "mo3 = orb_resemble(nbf1, nif1, mo1, nbf2, nif2, mo2, cross_S)"
 write(fid3,'(A)') 'noon = np.zeros(nif1)'
 write(fid3,'(A)') "py2fch('"//TRIM(fchname1)//"', nbf1, nif1, mo3, 'a', noon,&
                  & False)"
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

 call open_file(fchname, .true., fid)
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

