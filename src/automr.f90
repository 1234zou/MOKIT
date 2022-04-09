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
  write(iout,'(A)') 'AutoMR 1.2.3 :: MOKIT, release date: 2022-Apr-9'
  stop
 case('-h','-help','--help')
  write(iout,'(/,A)')  "Usage: automr [gjfname] >& [outname]"
  write(iout,'(A)')    "  Example 1 (in bash): automr a.gjf >& a.out &"
  write(iout,'(A)')    "  Example 2 (in dash): automr a.gjf >a.out 2>&1 &"
  write(iout,'(/,A)')  'Options:'
  write(iout,'(A)')    '  -h, -help, --help: Print this message and exit.'
  write(iout,'(A)')    '  -v, -V, --version: Print the version number of automr and exit.'
  write(iout,'(/,A)')  'Methods(#p ...):'
  write(iout,'(A)')    '  GVB, CASCI, CASSCF, DMRGCI, DMRGSCF, NEVPT2, NEVPT3, CASPT2,&
                      & CASPT2K, CASPT3,'
  write(iout,'(A)')    '  MRMP2, OVBMP2, SDSPT2, MRCISD, MCPDFT, FICMRCCSD, MkMRCCSD,&
                      & BWMRCCSD,'
  write(iout,'(A)')    '  MkMRCCSD(T), BWMRCCSD(T)'
  write(iout,'(/,A)')  'Frequently used keywords in MOKIT{}:'
  write(iout,'(A)')    '      HF_prog=Gaussian/PySCF/ORCA/PSI4'
  write(iout,'(A)')    '     GVB_prog=GAMESS/Gaussian'
  write(iout,'(A)')    '  CASSCF_prog=PySCF/OpenMolcas/ORCA/Molpro/GAMESS/Gaussian/BDF/PSI4/Dalton'
  write(iout,'(A)')    '   CASCI_prog=PySCF/OpenMolcas/ORCA/Molpro/GAMESS/Gaussian/BDF/PSI4/Dalton'
  write(iout,'(A)')    '  NEVPT2_prog=PySCF/OpenMolcas/ORCA/Molpro/BDF'
  write(iout,'(A)')    '  CASPT2_prog=OpenMolcas/Molpro/ORCA'
  write(iout,'(A)')    '  MCPDFT_prog=OpenMolcas/GAMESS'
  write(iout,'(A)')    '  MRCISD_prog=OpenMolcas/Molpro/ORCA/Gaussian/GAMESS/PSI4/Dalton'
  write(iout,'(A)')    '      CtrType=1/2/3 for uc-/ic-/FIC- MRCISD'
  write(iout,'(A,/)')  '    MRCC_prog=ORCA/NWChem'
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
 call get_paired_LMO()
 call do_gvb()        ! GVB
 call do_cas(.false.) ! CASCI/DMRG-CASCI
 call do_cas(.true.)  ! CASSCF/DMRG-CASSCF
 call do_mrpt2()      ! CASPT2/NEVPT2/SDSPT2/MRMP2
 call do_mrpt3()      ! CASPT3/NEVPT3
 call do_mrcisd()     ! uncontracted/ic-/FIC- MRCISD
 call do_mcpdft()     ! MC-PDFT
 call do_mrcc()       ! MRCC
 call do_sa_cas()     ! SA-CASSCF
 call do_PES_scan()   ! PES scan

 call fdate(data_string)
 write(iout,'(/,A)') 'Normal termination of AutoMR at '//TRIM(data_string)
 return
end subroutine automr

! generate PySCF input file .py from Gaussian .fch(k) file, and get paired LMOs
subroutine get_paired_LMO()
 use print_id, only: iout
 use mr_keyword, only: mo_rhf, ist, hf_fch, bgchg, chgname, dkh2_or_x2c
 use mol, only: nbf, nif, ndb, nacte, nacto, nacta, nactb, npair, npair0, nopen,&
  lin_dep, chem_core, ecp_core
 implicit none
 integer :: i, system, RENAME
 real(kind=8) :: unpaired_e
 character(len=24) :: data_string = ' '
 character(len=240) :: buf, proname, pyname, chkname, outname, fchname

 if(ist > 4) return ! no need for this subroutine
 write(iout,'(//,A)') 'Enter subroutine get_paired_LMO...'

 if(ist == 0) then
  write(iout,'(A)') 'ERROR in subroutine get_paired_LMO: ist=0. It should be&
                   & non-zero before this subroutine.'
  stop
 end if

 call calc_ncore() ! calculate the number of core orbitals from array core_orb
 write(iout,'(2(A,I0))') 'chem_core=', chem_core, ', ecp_core=', ecp_core
 i = index(hf_fch, '.fch', back=.true.)
 proname = hf_fch(1:i-1)

 if(mo_rhf) then
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
 return
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
 write(fid2,'(A)') "sort_pair('"//TRIM(loc_fch)//"','"//TRIM(hf_fch)//"',idx2,nopen,npair)"
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
 use mr_keyword, only: mem, nproc, hf_fch, tencycle, ON_thres
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
                           &mf.mo_coeff[1],S,',ON_thres,')'
 write(fid1,'(A)') 'alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)'
 write(fid1,'(A)') 'mf.mo_coeff = (alpha_coeff, beta_coeff)'
 write(fid1,'(A)') '# done transform'

 write(fid1,'(/,A)') '# save the UNO into .fch file'
 write(fid1,'(A)') "os.system('fch_u2r "//TRIM(hf_fch)//"')"
 write(fid1,'(A)') "os.rename('"//hf_fch(1:i-1)//"_r.fch', '"//TRIM(uno_fch)//"')"
 write(fid1,'(A)') "py2fch('"//TRIM(uno_fch)//"',nbf,nif,mf.mo_coeff[0],'a',noon,True)"
 write(fid1,'(A)') '# save done'
 close(fid1)
 return
end subroutine prt_uno_script_into_py

! print associated rotation into a given .py file
subroutine prt_assoc_rot_script_into_py(pyname)
 use print_id, only: iout
 use mol, only: chem_core, ecp_core
! use mol, only: natom, elem, nuc, chem_core, ecp_core
 use mr_keyword, only : localm, hf_fch, npair_wish
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
  write(iout,'(A)') 'ERROR in subroutine prt_assoc_rot_script_into_py: end-of-file detected.'
  write(iout,'(A)') 'File may be incomplete: '//TRIM(pyname)
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
 if(npair_wish > 0) write(fid1,'(A,I0,A1)') '  npair = min(npair,', npair_wish, ')'
 write(fid1,'(A)') '  idx1 = idx2 - npair'
 write(fid1,'(A)') '  idx4 = idx3 + npair'
 write(fid1,'(A)') '  occ_idx = range(idx1,idx2)'
 write(fid1,'(A)') '  vir_idx = range(idx3,idx4)'
! write(fid1,'(A)') '  print(idx1, idx2, idx3, idx4)'

 if(localm == 'pm') then ! Pipek-Mezey localization
  write(fid1,'(A)') "  occ_loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],&
                    &mol._bas[:,3],mol.cart,nbf,npair,mf.mo_coeff[0][:,occ_idx],&
                    &S,'mulliken')"
 else ! Boys localization
  write(fid1,'(A)') '  mo_dipole = dipole_integral(mol, mf.mo_coeff[0][:,occ_idx])'
  write(fid1,'(A)') '  occ_loc_orb = boys(nbf, npair, mf.mo_coeff[0][:,occ_idx], mo_dipole)'
 end if

 write(fid1,'(A)')  '  vir_loc_orb = assoc_rot(nbf, npair, mf.mo_coeff[0][:,occ_idx],&
                   & occ_loc_orb, mf.mo_coeff[0][:,vir_idx])'
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
 close(fid1)
 return
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
 return
end subroutine calc_ncore

