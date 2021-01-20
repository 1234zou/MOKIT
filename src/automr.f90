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
  write(iout,'(A,/)') " Example: automr a.gjf >& a.out &"
  stop
 end if

 call getarg(1, fname)
 select case(TRIM(fname))
 case('-v', '-V', '--version')
  write(iout,'(A)') 'AutoMR 1.2.2 :: MOKIT'
  stop
 end select

 i = index(fname, '.gjf', back=.true.)
 j = index(fname, '.fch', back=.true.)
 if(i/=0 .and. j/=0) then
  write(iout,'(A)') "ERROR in subroutine automr: both '.gjf' and '.fch' keys detected&
                   & in filename "//TRIM(fname)//'.'
  write(iout,'(A)') "Better to use a filename only with suffix '.gjf'."
  stop
 else if(i == 0) then
  write(iout,'(A)') "ERROR in subroutine automr: '.gjf' key not found in filename "//TRIM(fname)
  stop
 end if

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
 logical :: alive

 inquire(file=TRIM(fname), exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine automr: input file does not exist!'
  write(iout,'(A)') 'Filename = '//TRIM(fname)
  stop
 end if

 gjfname = fname
 call read_program_path() ! read in paths in $MOKIT_ROOT/program.info
 call parse_keyword()
 call check_kywd_compatible()

 call do_hf()
 call get_paired_LMO()
 call do_gvb()
 call do_cas(.false.)! CASCI/DMRG-CASCI
 call do_cas(.true.) ! CASSCF/DMRG-CASSCF
 call do_mrpt2()     ! CASPT2/NEVPT2/SDSPT2/MRMP2
 call do_mrpt3()     ! CASPT3/NEVPT3
 call do_mrcisd()    ! uncontracted/ic-/FIC- MRCISD
 call do_mcpdft()    ! MC-PDFT
 call do_mrcc()      ! ic-MRCC

 call fdate(data_string)
 write(iout,'(/,A)') 'Normal termination of AutoMR at '//TRIM(data_string)
 return
end subroutine automr

! perform RHF/UHF computation
subroutine do_hf()
 use print_id, only: iout
 use mol, only: natom, coor, elem, nuc, charge, mult, rhf_e, uhf_e
 use mr_keyword, only: skiphf, mem, nproc, basis, cart, gau_path, hf_fch, ist,&
  mo_rhf, bgchg, read_bgchg_from_gjf, gjfname, chgname, uno, dkh2_or_x2c, &
  vir_proj, prt_strategy, gau_path
 use util_wrapper, only: formchk
 implicit none
 integer :: i, system
 real(kind=8) :: ssquare = 0.0d0
 character(len=24) :: data_string = ' '
 character(len=240) :: rhf_gjfname, uhf_gjfname, chkname
 logical :: alive

 write(iout,'(//,A)') 'Enter subroutine do_hf...'
 call check_exe_exist(gau_path)

 if(skiphf) then
  write(iout,'(A)') 'Skip the RHF/UHF step...'
  call read_natom_from_fch(hf_fch, natom)
  allocate(coor(3,natom), elem(natom), nuc(natom))
  call read_elem_and_coor_from_fch(hf_fch, natom, elem, nuc, coor, charge, mult)
  if(bgchg) call read_bgchg_from_gjf(.true.)
  call fdate(data_string)
  write(iout,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)
  return
 end if

 call read_natom_from_gjf(gjfname, natom)
 allocate(coor(3,natom), elem(natom), nuc(natom))
 call read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
 if(bgchg) call read_bgchg_from_gjf(.false.)

 i = index(gjfname, '.gjf', back=.true.)
 rhf_gjfname = gjfname(1:i-1)//'_rhf.gjf'
 uhf_gjfname = gjfname(1:i-1)//'_uhf.gjf'

 if(mult == 1) then ! singlet, perform RHF and UHF
  call generate_hf_gjf(rhf_gjfname, natom, elem, coor, charge, mult, basis,&
                       .false., cart, dkh2_or_x2c, mem, nproc)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(rhf_gjfname))
  call perform_scf_and_read_e(gau_path, rhf_gjfname, rhf_e, ssquare)

  call generate_hf_gjf(uhf_gjfname, natom, elem, coor, charge, mult, basis,&
                       .true., cart, dkh2_or_x2c, mem, nproc)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_gjfname))
  call perform_scf_and_read_e(gau_path, uhf_gjfname, uhf_e, ssquare)

 else               ! not singlet, only perform UHF

  call generate_hf_gjf(uhf_gjfname, natom, elem, coor, charge, mult, basis,&
                       .true., cart, dkh2_or_x2c, mem, nproc)
  call perform_scf_and_read_e(gau_path, uhf_gjfname, uhf_e, ssquare)
 end if

 if(mult == 1) then
  write(iout,'(/,A,F18.8,1X,A,F7.3)') 'E(RHF) = ',rhf_e,'a.u., <S**2>=',0.0
  write(iout,'(A,F18.8,1X,A,F7.3)') 'E(UHF) = ',uhf_e,'a.u., <S**2>=',ssquare
 else
  write(iout,'(/,A,F18.8,1X,A,F7.3)') 'E(UHF) = ',uhf_e,'a.u., <S**2>=',ssquare
 end if

 if(rhf_e - uhf_e > 1.0d-4) then
  write(iout,'(A)') 'UHF energy is lower, choose UHF wave function.'
  ist = 1
  mo_rhf = .false.
  i = index(gjfname, '.gjf', back=.true.)
  chkname = gjfname(1:i-1)//'_uhf.chk'
  hf_fch = gjfname(1:i-1)//'_uhf.fch'
 else
  write(iout,'(A)') 'RHF/UHF is equal, or has little difference, choose RHF.'
  ist = 3
  vir_proj = .true.; mo_rhf = .true.; uno = .false.
  i = index(gjfname, '.gjf', back=.true.)
  chkname = gjfname(1:i-1)//'_rhf.chk'
  hf_fch = gjfname(1:i-1)//'_rhf.fch'
 end if

 call formchk(chkname, hf_fch)

 write(iout,'(A)') 'Strategy updated:'
 call prt_strategy()

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)
 return
end subroutine do_hf

! generate PySCF input file .py from Gaussian .fch(k) file, and get paired LMOs
subroutine get_paired_LMO()
 use print_id, only: iout
 use mr_keyword, only: mo_rhf, ist, hf_fch, bgchg, chgname, dkh2_or_x2c
 use mol, only: nbf, nif, ndb, nacte, nacto, nacta, nactb, npair, npair0, nopen,&
                lin_dep
 implicit none
 integer :: i, j, fid, system, RENAME
 character(len=24) :: data_string = ' '
 character(len=240) :: proname, pyname, chkname, outname, fchname

 if(ist > 4) return ! no need for this subroutine
 write(iout,'(//,A)') 'Enter subroutine get_paired_LMO...'

 if(ist == 0) then
  write(iout,'(A)') 'ERROR in subroutine get_paired_LMO: ist=0. It should be&
                   & non-zero before this subroutine.'
  stop
 end if

 i = index(hf_fch,'.fch',back=.true.)
 proname = hf_fch(1:i-1)
 chkname = hf_fch(1:i-1)//'_proj.chk' ! this is PySCF chk file, not Gaussian
 fchname = hf_fch(1:i-1)//'_uno.fch'

 if(mo_rhf) then
  write(iout,'(A)') 'One set of MOs: invoke RHF virtual MO projection ->&
                   & localization -> paring.'

  i = system('bas_fch2py '//TRIM(hf_fch))
  pyname = TRIM(proname)//'_proj_loc_pair.py'
  call delete_file(pyname) ! if already exists, delete it
  i = RENAME(TRIM(proname)//'.py', TRIM(pyname))
  if(dkh2_or_x2c) call add_X2C_into_py(pyname)

  call prt_rhf_proj_script_into_py(pyname)
  call prt_auto_pair_script_into_py(pyname)
  i = system('python '//TRIM(pyname)//' >& '//TRIM(proname)//'_proj_loc_pair.out')
  call delete_file(chkname)

 else
  if(ist == 1) then
   write(iout,'(A)') 'Two sets of MOs, ist=1, invoke UNO associated rotation.'
  else if(ist == 2) then
   write(iout,'(A)') 'Two sets of MOs, ist=2, invoke UNO generation.'
  end if

  i = system('bas_fch2py '//TRIM(hf_fch)//' -uhf')
  if(ist == 1) then
   pyname = TRIM(proname)//'_uno_asrot.py'
   outname = TRIM(proname)//'_uno_asrot.out'
  else
   pyname = TRIM(proname)//'_uno.py'
   outname = TRIM(proname)//'_uno.out'
  end if
  call delete_file(pyname) ! if already exists, delete it
  i = RENAME(TRIM(proname)//'.py', TRIM(pyname))
  if(dkh2_or_x2c) call add_X2C_into_py(pyname)
  call prt_uno_script_into_py(pyname)

  if(ist == 1) call prt_assoc_rot_script_into_py(pyname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(pyname))
  i = system('python '//TRIM(pyname)//' >& '//TRIM(outname))

  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine get_paired_LMO: PySCF job fails.'
   write(iout,'(A)') 'Please check file '//TRIM(outname)
   stop
  end if

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
 use mol, only: ndb, npair
 implicit none
 integer :: i, fid1, fid2, RENAME
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
 write(fid2,'(A)') "mol2.basis = 'STO-6G'"
 write(fid2,'(A)') 'mol2.build()'
 if(dkh2_or_x2c) then
  write(fid2,'(A)') 'mf2 = scf.RHF(mol2).x2c()'
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
 write(fid2,'(A)') "py2fch('"//TRIM(proj_fch)//"', nbf, nif, mf.mo_coeff, Sdiag, 'a', mf.mo_occ)"
 write(fid2,'(A)') '# save done'

 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 return
end subroutine prt_rhf_proj_script_into_py

! print localization and automatically pairing information into a given .py file
subroutine prt_auto_pair_script_into_py(pyname)
 use print_id, only: iout
 use mr_keyword, only: localm, hf_fch
 use mol, only: ndb, npair
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1, loc_fch
 character(len=240), intent(in) :: pyname

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
 write(fid2,'(A)') 'idx2 = np.sum(mf.mo_occ==2)'
 write(fid2,'(A)') 'idx1 = idx2 - npair'
 write(fid2,'(A)') 'idx3 = 2*idx2 - idx1'
 write(fid2,'(A)') 'occ_idx = range(idx1,idx2)'
 write(fid2,'(A)') 'vir_idx = range(idx2,idx3)'

 if(localm == 'pm') then ! Pipek-Mezey localization
  write(fid2,'(A)') "occ_loc_orb = pm(mol.nbas, mol._bas[:,0], mol._bas[:,1], &
                   & mol._bas[:,3], mol.cart, nbf, npair, mf.mo_coeff[:,occ_idx],&
                   & S, 'mulliken')"
  write(fid2,'(A)') "vir_loc_orb = pm(mol.nbas, mol._bas[:,0], mol._bas[:,1], &
                   & mol._bas[:,3], mol.cart, nbf, npair, mf.mo_coeff[:,vir_idx],&
                   & S, 'mulliken')"
 else ! Boys localization
  write(fid2,'(A)') 'mo_dipole = dipole_integral(mol, mf.mo_coeff[:,occ_idx])'
  write(fid2,'(A)') 'occ_loc_orb = boys(nbf, npair, mf.mo_coeff[:,occ_idx], mo_dipole)'
  write(fid2,'(A)') 'mo_dipole = dipole_integral(mol, mf.mo_coeff[:,vir_idx])'
  write(fid2,'(A)') 'vir_loc_orb = boys(nbf, npair, mf.mo_coeff[:,vir_idx], mo_dipole)'
 end if

 write(fid2,'(A)') 'mf.mo_coeff[:,occ_idx] = occ_loc_orb.copy()'
 write(fid2,'(A)') 'mf.mo_coeff[:,vir_idx] = vir_loc_orb.copy()'
 write(fid2,'(A)') '# localization done'

 write(fid2,'(/,A)') '# pair the active orbitals'
 write(fid2,'(A)') 'mo_dipole = dipole_integral(mol, mf.mo_coeff)'
 write(fid2,'(A)') 'ncore = idx1'
 write(fid2,'(A)') 'nopen = np.sum(mf.mo_occ==1)'
 write(fid2,'(A)') 'nalpha = np.sum(mf.mo_occ > 0)'
 write(fid2,'(A)') 'nvir_lmo = npair'
 write(fid2,'(A)') 'alpha_coeff = pair_by_tdm(ncore, npair, nopen, nalpha, nvir_lmo, nbf, nif, mf.mo_coeff, mo_dipole)'
 write(fid2,'(A)') 'mf.mo_coeff = alpha_coeff.copy()'
 write(fid2,'(A)') '# pair done'

 write(fid2,'(/,A)') '# save the paired LMO into .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(hf_fch)//"', '"//TRIM(loc_fch)//"')"
 write(fid2,'(A)') "py2fch('"//TRIM(loc_fch)//"', nbf, nif, mf.mo_coeff, Sdiag, 'a', mf.mo_occ)"
 write(fid2,'(A)') '# save done'

 write(fid2,'(/,A)') "f = open('uno.out', 'w+')"
 write(fid2,'(A)') "f.write('nbf=%i\n' %nbf)"
 write(fid2,'(A)') "f.write('nif=%i\n\n' %nif)"
 write(fid2,'(A)') "f.write('ndb=%i\n\n' %ncore)"
 write(fid2,'(A)') "f.write('idx=%i %i %i' %(ncore+1,nalpha+npair+1,nopen))"
 write(fid2,'(A)') 'f.close()'
 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 return
end subroutine prt_auto_pair_script_into_py

! print UNO script into a given .py file
subroutine prt_uno_script_into_py(pyname)
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, hf_fch, tencycle
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
 write(fid1,'(A)') 'idx, noon, alpha_coeff = uno(nbf, nif, na, nb, mf.mo_coeff[0], mf.mo_coeff[1], S)'
 write(fid1,'(A)') 'alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)'
 write(fid1,'(A)') 'mf.mo_coeff = (alpha_coeff, beta_coeff)'
 write(fid1,'(A)') '# done transform'

 write(fid1,'(/,A)') '# save the UNO into .fch file'
 write(fid1,'(A)') "os.system('fch_u2r "//TRIM(hf_fch)//"')"
 write(fid1,'(A)') "os.rename('"//hf_fch(1:i-1)//"_r.fch', '"//TRIM(uno_fch)//"')"
 write(fid1,'(A)') "py2fch('"//TRIM(uno_fch)//"', nbf, nif, mf.mo_coeff[0], Sdiag, 'a', noon)"
 write(fid1,'(A)') '# save done'

 close(fid1)
 return
end subroutine prt_uno_script_into_py

! print associated rotation into a given .py file
subroutine prt_assoc_rot_script_into_py(pyname)
 use print_id, only: iout
 use mr_keyword, only : localm, hf_fch, npair_wish
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1, uno_fch, assoc_fch
 character(len=240), intent(in) :: pyname

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
  write(fid2,'(A)') 'from lo import pm'
 else
  write(fid2,'(A)') 'from lo import boys'
  write(fid2,'(A)') 'from pyscf.lo.boys import dipole_integral'
 end if
 write(fid2,'(A)') 'from assoc_rot import assoc_rot'
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
 write(fid1,'(A)') 'idx2 = idx[0] + npair - 1'
 write(fid1,'(A)') 'idx3 = idx2 + idx[2]'
 if(npair_wish > 0) write(fid1,'(A,I0,A1)') 'npair = min(npair,', npair_wish, ')'
 write(fid1,'(A)') 'idx1 = idx2 - npair'
 write(fid1,'(A)') 'idx4 = idx3 + npair'
 write(fid1,'(A)') 'occ_idx = range(idx1,idx2)'
 write(fid1,'(A)') 'vir_idx = range(idx3,idx4)'
 write(fid1,'(A)') 'print(idx1, idx2, idx3, idx4)'

 if(localm == 'pm') then ! Pipek-Mezey localization
  write(fid1,'(A)') "occ_loc_orb = pm(mol.nbas, mol._bas[:,0], mol._bas[:,1], &
                    & mol._bas[:,3], mol.cart, nbf, npair, mf.mo_coeff[0][:,occ_idx],&
                    & S, 'mulliken')"
 else ! Boys localization
  write(fid1,'(A)') 'mo_dipole = dipole_integral(mol, mf.mo_coeff[0][:,occ_idx])'
  write(fid1,'(A)') 'occ_loc_orb = boys(nbf, npair, mf.mo_coeff[0][:,occ_idx], mo_dipole)'
 end if

 write(fid1,'(A)')  'vir_loc_orb = assoc_rot(nbf, npair, mf.mo_coeff[0][:,occ_idx],&
                   & occ_loc_orb, mf.mo_coeff[0][:,vir_idx])'
 write(fid1,'(A)') 'mf.mo_coeff[0][:,occ_idx] = occ_loc_orb.copy()'
 write(fid1,'(A)') 'mf.mo_coeff[0][:,vir_idx] = vir_loc_orb.copy()'
 write(fid1,'(A)') '# localization done'

 write(fid1,'(/,A)') '# save associated rotation MOs into .fch(k) file'
 write(fid1,'(A)') "copyfile('"//TRIM(uno_fch)//"', '"//TRIM(assoc_fch)//"')"
 write(fid1,'(A)') 'noon = np.zeros(nif)'
 write(fid1,'(A)') "py2fch('"//TRIM(assoc_fch)//"', nbf, nif, mf.mo_coeff[0], Sdiag, 'a', noon)"
 close(fid1)
 return
end subroutine prt_assoc_rot_script_into_py

! print CASCI/DMRG-CASCI or CASSCF/DMRG-CASSCF script into a given .py file
subroutine prt_cas_script_into_py(pyname, gvb_fch, scf)
 use mol, only: npair, nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, casci, dmrgci, casscf, dmrgscf, maxM,&
                       hardwfn, crazywfn, hf_fch, casnofch, CIonly, &
                       casscf_force, dkh2_or_x2c
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1, cmofch
 character(len=240), intent(in) :: pyname, gvb_fch
 logical, intent(in) :: scf

 pyname1 = TRIM(pyname)//'.tmp'
 i = index(pyname, '.py', back=.true.)

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(dmrgci .or. dmrgscf) then
  write(fid2,'(A)') 'from pyscf import mcscf, dmrgscf, lib'
 else
  write(fid2,'(A)') 'from pyscf import mcscf, lib'
 end if
 write(fid2,'(A)') 'from py2fch import py2fch'
 write(fid2,'(A)') 'from shutil import copyfile'
 write(fid2,'(A,/)') 'import numpy as np'
 if(dmrgci .or. dmrgscf) then
  write(fid2,'(A,I0,A)') "dmrgscf.settings.MPIPREFIX ='mpirun -n ",nproc,"'"
 end if
 write(fid2,'(A,I0,A1,/)') 'lib.num_threads(',nproc,')'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:9) == 'mf.kernel') exit
  write(fid2,'(A)') TRIM(buf)
 end do
 write(fid2,'(A,I0,A)') 'mf.max_memory = ', mem*1000, ' # MB'
 write(fid2,'(A)') 'mf.kernel()'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')

 ! mem*500 is in fact mem*1000/2. The mc.max_memory and fcisolver.max_memory seem
 ! not to share common memory, they used two memory, so I have to make them half
 if(scf) then ! CASSCF/DMRG-CASSCF
  if(casscf) then
   if(dkh2_or_x2c) then
    write(fid2,'(A,3(I0,A))') 'mc = mcscf.CASSCF(mf,',nacto,',(',nacta,',',nactb,')).x2c()'
   else
    write(fid2,'(A,3(I0,A))') 'mc = mcscf.CASSCF(mf,',nacto,',(',nacta,',',nactb,'))'
   end if
   write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ',mem*500,' # MB'
  else ! DMRG-CASSCF
   if(dkh2_or_x2c) then
    write(fid2,'(A,3(I0,A))') 'mc = dmrgscf.DMRGSCF(mf,',nacto,',(',nacta,',',nactb,')).x2c()'
   else
    write(fid2,'(A,3(I0,A))') 'mc = dmrgscf.DMRGSCF(mf,',nacto,',(',nacta,',',nactb,'))'
   end if
   write(fid2,'(A,I0)') 'mc.fcisolver.maxM = ', maxM
   write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ',CEILING(DBLE(mem)/DBLE((2*nproc))),' # GB'
  end if
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*500, ' # MB'
  write(fid2,'(A)') 'mc.max_cycle = 200'
 else ! CASCI/DMRG-CASCI
  if(dkh2_or_x2c) then
   write(fid2,'(A,3(I0,A))') 'mc = mcscf.CASCI(mf,',nacto,',(',nacta,',',nactb,')).x2c()'
  else
   write(fid2,'(A,3(I0,A))') 'mc = mcscf.CASCI(mf,',nacto,',(',nacta,',',nactb,'))'
  end if
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*500, ' # MB'
  if(casci) then
   write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*500, ' # MB'
  else
   write(fid2,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=', maxM, ')'
   write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ',CEILING(DBLE(mem)/DBLE((2*nproc))),' # GB'
  end if
 end if

 write(fid2,'(A,I0)') 'mc.fcisolver.spin = ', nacta-nactb

 if(hardwfn) then
  write(fid2,'(A,I0,A)') 'mc.fix_spin_(ss=',nacta-nactb,')'
 else if(crazywfn) then
  write(fid2,'(A,I0,A)') 'mc.fix_spin_(ss=',nacta-nactb,')'
  write(fid2,'(A)') 'mc.fcisolver.level_shift = 0.2'
  write(fid2,'(A)') 'mc.fcisolver.pspace_size = 1200'
  write(fid2,'(A)') 'mc.fcisolver.max_space = 100'
 end if

 ! For CASCI/CASSCF, NOs will be generated; but for DMRGCI/DMRGSCF, both
 ! the (canonical) MOs and NOs will be generated.
 ! Since DMRG is not invariant to unitary rotations of orbitals, I hope
 ! the (canonical) MOs to be used in DMRG-NEVPT2/CASPT2 computations.
 if(.not.(dmrgci .or. dmrgscf)) write(fid2,'(A)') 'mc.natorb = True'
 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.kernel()'

 if(dmrgci .or. dmrgscf) then
  i = index(casnofch, '_NO', back=.true.)
  cmofch = casnofch(1:i)//'CMO.fch'
  write(fid2,'(/,A)') '# save CMOs into .fch file'
  write(fid2,'(A)') "copyfile('"//TRIM(gvb_fch)//"', '"//TRIM(cmofch)//"')"
  write(fid2,'(A)') 'noon = np.zeros(nif)'
  write(fid2,'(A)') "py2fch('"//TRIM(cmofch)//"',nbf,nif,mc.mo_coeff,Sdiag,'a',noon)"

  write(fid2,'(/,A)') 'mc.natorb = True'
  write(fid2,'(A)') 'mc.kernel()'
 end if

 write(fid2,'(/,A)') '# save NOs into .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(gvb_fch)//"', '"//TRIM(casnofch)//"')"
 write(fid2,'(A)') 'noon = np.zeros(nif)'
 write(fid2,'(A)') "py2fch('"//TRIM(casnofch)//"',nbf,nif,mc.mo_coeff,Sdiag,'a',noon)"
 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 return
end subroutine prt_cas_script_into_py

! print a CASCI or CASSCF .gjf file
subroutine prt_cas_gjf(gjfname, nacto, nacte, scf, force)
 use mr_keyword, only: mem, nproc, dkh2_or_x2c
 implicit none
 integer :: i, fid
 integer, intent(in) :: nacto, nacte
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: scf, force

 i = index(gjfname, '.gjf', back=.true.)
 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(A,I0,A,I0,A)',advance='no') '#p CAS(',nacte,',',nacto,')'
 if(dkh2_or_x2c) then
  write(fid,'(A)',advance='no') ' chkbasis nosymm guess=read geom=allcheck&
                                  & int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A)',advance='no') ' chkbasis nosymm guess=read geom=allcheck int=nobasistransform'
 end if
 if(force) write(fid,'(A)',advance='no') ' force'

 if(scf) then ! CASSCF
  write(fid,'(A,/)') ' scf(maxcycle=128)'
 else         ! CASCI
  write(fid,'(A,/)') ' scf(maxcycle=-2)'
  ! to obtain CASCI NOs, we need to use -2. Using -1 only calculates the CASCI energy.
 end if

 write(fid,'(A)') '--Link1--'
 write(fid,'(A)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(A)',advance='no') '#p chkbasis nosymm guess(read,only,save,NaturalOrbitals) geom=allcheck'
 if(dkh2_or_x2c) then
  write(fid,'(A)') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A)') ' int=nobasistransform'
 end if

 close(fid)
 return
end subroutine prt_cas_gjf

subroutine prt_cas_gms_inp(inpname, ncore, scf)
 use mol, only: nacto, nacte, charge, mult
 use mr_keyword, only: mem, nproc, hardwfn, crazywfn, dkh2_or_x2c, cart
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: ncore
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 if(scf) then   ! CASSCF
  write(fid2,'(A)',advance='no') ' $CONTRL SCFTYP=MCSCF RUNTYP=ENERGY ICHARG='
 else           ! CASCI
  write(fid2,'(A)',advance='no') ' $CONTRL SCFTYP=NONE CITYP=ALDET RUNTYP=ENERGY ICHARG='
 end if
 write(fid2,'(2(I0,A))') charge, ' MULT=', mult, ' NOSYM=1'

 write(fid2,'(A)',advance='no') '  ICUT=11'
 if(dkh2_or_x2c) write(fid2,'(A)',advance='no') ' RELWFN=DK'

 if(.not. cart) then
  write(fid2,'(A)') ' ISPHER=1 $END'
 else
  write(fid2,'(A)') ' $END'
 end if

 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',CEILING(DBLE(mem*125)/DBLE(nproc)), ' $END'

 if(scf) then   ! CASSCF
  write(fid2,'(A)',advance='no') ' $DET'
 else
  write(fid2,'(A)',advance='no') ' $CIDET'
 end if

 write(fid2,'(3(A,I0),A)',advance='no') ' NCORE=',ncore,' NELS=',nacte,' NACT=',nacto,' ITERMX=500'
 if(hardwfn) then
  write(fid2,'(A)') ' NSTATE=5 $END'
 else if(crazywfn) then
  write(fid2,'(A)') ' NSTATE=10 $END'
 else
  write(fid2,'(A)') ' $END'
 end if

 if(scf) write(fid2,'(A)') ' $MCSCF MAXIT=200 $END'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:7) == '$GUESS') exit
 end do

 BACKSPACE(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine prt_cas_gms_inp

! print CASCI/CASSCF keywords in to a given (Open)Molcas input file
subroutine prt_cas_molcas_inp(inpname, scf)
 use print_id, only: iout
 use mol, only: charge, mult, nacte, nacto
 use mr_keyword, only: maxM, dmrgci, dmrgscf
 implicit none
 integer :: i, fid1, fid2, RENAME
 logical, intent(in) :: scf
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(2:7) == 'SEWARD') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_cas_molcas_inp: no 'SEWARD'&
                   & found in file "//TRIM(inpname)
  stop
 end if

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:4) == "SCF") exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_cas_molcas_inp: no 'SCF'&
                   & found in file "//TRIM(inpname)
  stop
 end if

 write(fid2,'(A)') "&RASSCF"
 write(fid2,'(A,I0)') 'Spin = ', mult
 write(fid2,'(A,I0)') 'Charge = ', charge
 write(fid2,'(A,I0)') 'nActEl = ', nacte
 write(fid2,'(A,I0)') 'RAS2 = ', nacto
 if(.not. scf) write(fid2,'(A)') 'CIonly'
 i = index(inpname, '.input', back=.true.)
 write(fid2,'(A,/)') 'FILEORB = '//inpname(1:i-1)//'.INPORB'

 if(dmrgci .or. dmrgscf) then
  write(fid2,'(A)') 'DMRG'
  write(fid2,'(A)') 'RGinput'
  write(fid2,'(A)') ' conv_thresh = 1E-7'
  write(fid2,'(A)') ' nsweeps = 5'
  write(fid2,'(A,I0)') ' max_bond_dimension = ', MaxM
  write(fid2,'(A)') 'endRG'
 end if

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine prt_cas_molcas_inp

! print CASCI/CASSCF keywords in to a given ORCA .inp file
subroutine prt_cas_orca_inp(inpname, scf)
 use print_id, only: iout
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc, dkh2_or_x2c
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 read(fid1,'(A)') buf   ! skip nproc
 read(fid1,'(A)') buf   ! skip memory
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid2,'(A,I0,A)') '%maxcore ', CEILING(1000.0d0*DBLE(mem)/DBLE(nproc))

 read(fid1,'(A)') buf   ! skip '!' line
 if(scf) then
  write(fid2,'(A)') '! TightSCF'
 else
  write(fid2,'(A)') '! TightSCF noiter'
 end if

 if(dkh2_or_x2c) then
  write(fid2,'(A)') '%rel'
  write(fid2,'(A)') ' method DKH'
  write(fid2,'(A)') ' order 2'
  write(fid2,'(A)') 'end'
 end if

 if(scf) then ! CASSCF
  write(fid2,'(A)') '%casscf'
  write(fid2,'(A,I0)') ' nel ', nacte
  write(fid2,'(A,I0)') ' norb ', nacto
  write(fid2,'(A)') ' maxiter 200'
  write(fid2,'(A)') ' ActOrbs NatOrbs'
  write(fid2,'(A)') 'end'
 else ! CASCI
  write(fid2,'(A)') '%mrci'
  write(fid2,'(A)') ' tsel 0.0'
  write(fid2,'(A)') ' tpre 0.0'
  write(fid2,'(A)') ' Etol 1e-7'
  write(fid2,'(A)') ' Rtol 1e-7'
  write(fid2,'(A)') ' MaxIter 100'
  write(fid2,'(A)',advance='no') ' NewBlock 1 * nroots 1 excitations none refs cas('
  write(fid2,'(2(I0,A))') nacte, ',', nacto, ') end end'
  write(fid2,'(A)') ' doNatOrbs 2'
  write(fid2,'(A)') 'end'
 end if
 ! Q: Why not use %casscf plus NoIter for CASCI?
 ! A: The CASCI NOONs cannot be obtained in that way, so I have to use %mrci

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine prt_cas_orca_inp

! print CASCI/CASSCF keywords into a given Molpro input file
subroutine prt_cas_molpro_inp(inpname, scf, force)
 use print_id, only: iout
 use mol, only: ndb, npair, npair0, nacto
 use mr_keyword, only: mem, nproc, dkh2_or_x2c
 implicit none
 integer :: i, fid, nclosed, nocc
 character(len=240) :: buf, orbfile, put
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf, force

 orbfile = inpname
 call convert2molpro_fname(orbfile, '.a')

 put = ' '
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 BACKSPACE(fid)
 read(fid,'(A)') put

 if(put(1:4) /= '{put') then
  close(fid)
  write(iout,'(A)') 'ERROR in subroutine prt_cas_molpro_inp: wrong content found&
                   & in the final line of file '//TRIM(inpname)
  stop
 end if

 rewind(fid)
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:2) == 'wf') exit
 end do ! for while

 write(fid,'(A)') '{matrop;'
 write(fid,'(A)') 'read,mo,ORB,file='//TRIM(orbfile)//';'
 write(fid,'(A)') 'save,mo,2140.2,ORBITALS}'

 nclosed = ndb + npair - npair0
 nocc = nclosed + nacto

 ! Note: we need 'NoExtra' to completely close symmetry.
 ! Otherwise the CASCI energy is slightly different to that of other programs
 if(scf) then
  write(fid,'(2(A,I0),A)') '{CASSCF;closed,',nclosed,';occ,',nocc,';NoExtra}'
 else
  write(fid,'(2(A,I0),A)') '{CASSCF;closed,',nclosed,';occ,',nocc,&
                           ';DONT,ORBITAL;NoExtra}'
 end if

 if(force) write(fid,'(A)') 'Forces'
 write(fid,'(A)') TRIM(put)
 close(fid)
 return
end subroutine prt_cas_molpro_inp

! print CASCI/CASSCF keywords into a given BDF input file
subroutine prt_cas_bdf_inp(inpname, scf, force)
 use print_id, only: iout
 use mol, only: charge, mult, ndb, npair, npair0, nacto, nacte
 implicit none
 integer :: i, nclosed, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf, force

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == '$SCF') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_cas_bdf_inp: '$SCF' not found in&
                  & file "//TRIM(inpname)
  close(fid)
  stop
 end if

 nclosed = ndb + npair - npair0
 BACKSPACE(fid)
 write(fid,'(A)') '$MCSCF'
 write(fid,'(A)') 'CheckLin'
 write(fid,'(A)') 'TolLin'
 write(fid,'(A)') '1.D-6'
 write(fid,'(A,/,I0)') 'Charge', charge
 write(fid,'(A,/,I0)') 'Spin', mult
 write(fid,'(A,/,I0)') 'Close', nclosed
 write(fid,'(A,/,I0)') 'Active', nacto
 write(fid,'(A,/,I0)') 'Actel', nacte
 write(fid,'(A)') 'Guess'
 write(fid,'(A)') 'read'
 if(.not. scf) write(fid,'(A)') 'CASCI'
 write(fid,'(A)') '$END'

 if(force) then
  write(fid,'(/,A)') '$GRAD'
  write(fid,'(A,/,A1)') 'nrootgrad','1'
  write(fid,'(A)') '$END'
 end if

 close(fid)
 return
end subroutine prt_cas_bdf_inp

! print NEVPT2 script into a given .py file
subroutine prt_nevpt2_script_into_py(pyname)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, casnofch, casci, casscf, maxM, X2C
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1
 character(len=240), intent(in) :: pyname

 pyname1 = TRIM(pyname)//'.tmp'
 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 write(fid2,'(A)') 'from pyscf import mcscf, mrpt, lib'

 if(.not. (casci .or. casscf)) then
  write(fid2,'(A,/)') 'from pyscf import dmrgscf'
  write(fid2,'(A,I0,A)') "dmrgscf.settings.MPIPREFIX ='mpirun -n ",nproc,"'"
 end if
 write(fid2,'(A,I0,A1)') 'lib.num_threads(',nproc,')'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')

 write(fid2,'(A)') '# generate CASCI wfn'
 if(X2C) then
  write(fid2,'(A,3(I0,A))') 'mc = mcscf.CASCI(mf,',nacto,',(',nacta,',',nactb,')).x2c()'
 else
  write(fid2,'(A,3(I0,A))') 'mc = mcscf.CASCI(mf,',nacto,',(',nacta,',',nactb,'))'
 end if

 write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*500, ' # MB'

 if(casci .or. casscf) then
  write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*500, ' # MB'
 else
  write(fid2,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=', maxM, ')'
  write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ', CEILING(DBLE(mem)/2.0d0), ' # GB'
 end if

 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.kernel()'
 if(casci .or. casscf) then
  write(fid2,'(A)') 'mrpt.NEVPT(mc).kernel()'
 else
  write(fid2,'(A,I0,A)') 'mrpt.NEVPT(mc).compress_approx(maxM=',maxM,').kernel()'
 end if
 close(fid2)

 i = RENAME(pyname1, pyname)
 return
end subroutine prt_nevpt2_script_into_py

! print NEVPT2 keywords in to a given ORCA .inp file
subroutine prt_nevpt2_orca_inp(inpname)
 use print_id, only: iout
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc, DKH2, X2C, CIonly
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 read(fid1,'(A)') buf   ! skip nproc
 read(fid1,'(A)') buf   ! skip memory
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid2,'(A,I0,A)') '%maxcore ', CEILING(1000.0d0*DBLE(mem)/DBLE(nproc))

 read(fid1,'(A)') buf   ! skip '!' line
 if(CIonly) then
  write(fid2,'(A)') '! TightSCF noiter'
 else
  write(fid2,'(A)') '! TightSCF'
 end if

 if(DKH2) then
  write(fid2,'(A)') '%rel'
  write(fid2,'(A)') ' method DKH'
  write(fid2,'(A)') ' order 2'
  write(fid2,'(A)') 'end'
 else if(X2C) then
  write(iout,'(A)') 'ERROR in subroutine prt_nevpt2_orca_inp: NEVPT2 with X2C&
                   & is not supported in ORCA.'
  write(iout,'(A)') 'You can specify NEVPT2_prog=Molpro or OpenMolcas.'
  stop
 end if

 write(fid2,'(A)') '%casscf'
 write(fid2,'(A,I0)') ' nel ', nacte
 write(fid2,'(A,I0)') ' norb ', nacto
 write(fid2,'(A)') ' PTMethod SC_NEVPT2'
 write(fid2,'(A)') 'end'
 write(fid2,'(A)') '%method'
 write(fid2,'(A)') ' FrozenCore FC_NONE'
 write(fid2,'(A)') 'end'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine prt_nevpt2_orca_inp

! print NEVPT2 keywords into OpenMolcas .input file
! It seems that OpenMolcas does not support CASSCF-NEVPT2. So I have to use
! DMRG-NEVPT2.
subroutine prt_nevpt2_molcas_inp(inpname)
 use print_id, only: iout
 use mr_keyword, only: nproc, CIonly, maxM
 use mol, only: nacte, nacto, charge, mult
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == "&SCF") exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_nevpt2_molcas_inp: no '&SCF'&
                   & found in file "//TRIM(inpname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if
 close(fid1,status='delete')

 write(fid2,'(/,A)') "&DMRGSCF"
 write(fid2,'(A)') 'ActiveSpaceOptimizer=QCMaquis'
 write(fid2,'(A)') 'DMRGSettings'
 write(fid2,'(A,I0)') ' max_bond_dimension = ', maxM
 write(fid2,'(A)') ' nsweeps = 5'
 write(fid2,'(A)') 'EndDMRGSettings'
 write(fid2,'(A)') 'OOptimizationSettings'
 write(fid2,'(A,I0)') 'Charge = ', charge
 write(fid2,'(A,I0)') 'Spin = ', mult
 write(fid2,'(A,I0,A)') 'nActEl = ', nacte, ',0,0'
 write(fid2,'(A,I0)') 'RAS2 = ', nacto
 i = index(inpname,'.input',back=.true.)
 write(fid2,'(A)') 'FILEORB = '//inpname(1:i-1)//'.INPORB'
 if(CIonly) write(fid2,'(A)') 'CIonly'
 write(fid2,'(A)') 'NEVPT2Prep'
 write(fid2,'(A)') 'EvRDM'
 write(fid2,'(A)') 'EndOOptimizationSettings'

 write(fid2,'(/,A)') "&MOTRA"
 write(fid2,'(A)') 'Frozen = 0'
 write(fid2,'(A)') 'HDF5'

 write(fid2,'(/,A)') "&NEVPT2"
 write(fid2,'(A)') 'Frozen = 0'
 close(fid2)
 i = RENAME(inpname1, inpname)
 return
end subroutine prt_nevpt2_molcas_inp

! print CASTP2 keywords into OpenMolcas .input file
subroutine prt_caspt2_molcas_inp(inputname)
 use print_id, only: iout
 use mr_keyword, only: CIonly, maxM, dmrgci, dmrgscf
 use mol, only: nacte, nacto, charge, mult
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inputname1
 character(len=240), intent(in) :: inputname

 inputname1 = TRIM(inputname)//'.tmp'
 open(newunit=fid1,file=TRIM(inputname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inputname1),status='replace')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(2:7) == 'SEWARD') exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_caspt2_molcas_inp: no 'SEWARD'&
                   & found in file "//TRIM(inputname)
  stop
 end if

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  if(buf(1:4) == "&SCF") exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid1,status='delete')

 write(fid2,'(/,A)') "&RASSCF"
 write(fid2,'(A,I0)') 'Charge = ', charge
 write(fid2,'(A,I0)') 'Spin = ', mult
 write(fid2,'(A,I0,A)') 'nActEl= ', nacte, ' 0 0'
 write(fid2,'(A,I0)') 'RAS2 = ', nacto
 i = index(inputname,'.input',back=.true.)
 write(fid2,'(A)') 'FILEORB = '//inputname(1:i-1)//'.INPORB'

 if(CIonly) then
  if(dmrgci) then
   write(iout,'(A)') 'ERROR in subroutine prt_caspt2_molcas_inp: CIonly is not&
                    & allowed in DMRG-CASPT2.'
   write(iout,'(A)') 'It must be based on converged (DMRG-)CASSCF orbitals.'   
   stop
  else
   write(fid2,'(A)') 'CIonly'
  end if
 end if

 if(dmrgscf) then
  write(fid2,'(A,I0)') 'DMRG = ', maxM
  write(fid2,'(A)') '3RDM'
 end if

 write(fid2,'(/,A)') "&CASPT2"
 write(fid2,'(A)') 'MultiState= 1 1'
 if(dmrgscf) write(fid2,'(A)') 'CheMPS2'
 write(fid2,'(A,/)') 'Frozen= 0'

 close(fid2)
 i = RENAME(inputname1, inputname)
 return
end subroutine prt_caspt2_molcas_inp

! print MRMP2 keywords into GAMESS .inp file
subroutine prt_mrmp2_gms_inp(inpname)
 use mol, only: charge, mult, ndb, npair, npair0, nacte, nacto
 use mr_keyword, only: mem, nproc, hardwfn, crazywfn, DKH2, cart
 implicit none
 integer :: i, ncore, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 ncore = ndb + npair - npair0

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')
 write(fid2,'(A)',advance='no') ' $CONTRL SCFTYP=MCSCF RUNTYP=ENERGY ICHARG='
 write(fid2,'(2(I0,A))') charge, ' MULT=', mult, ' NOSYM=1'

 write(fid2,'(A)',advance='no') '  ICUT=11 MPLEVL=2'
 if(DKH2) write(fid2,'(A)',advance='no') ' RELWFN=DK'

 if(.not. cart) then
  write(fid2,'(A)') ' ISPHER=1 $END'
 else
  write(fid2,'(A)') ' $END'
 end if

 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',CEILING(DBLE(mem*125)/DBLE(nproc)), ' $END'

 write(fid2,'(A)',advance='no') ' $DET'
 write(fid2,'(3(A,I0))',advance='no') ' NCORE=',ncore,' NELS=',nacte,' NACT=',nacto

 if(hardwfn) then
  write(fid2,'(A)',advance='no') ' NSTATE=5'
 else if(crazywfn) then
  write(fid2,'(A)',advance='no') ' NSTATE=10'
 end if
 write(fid2,'(A)') ' ITERMX=500 $END'
 write(fid2,'(A,I0,A)') ' $DETPT NVAL=', ncore,' $END'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:7) == '$GUESS') exit
 end do

 BACKSPACE(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(inpname1, inpname)
 return
end subroutine prt_mrmp2_gms_inp

! print MRCISD keywords into OpenMolcas .input file
subroutine prt_mrcisd_molcas_inp(inpname)
 use print_id, only: iout
 use mol, only: nif, ndb, nopen, nacta, nactb, npair, npair0, charge, mult
 use mr_keyword, only: CtrType
 implicit none
 integer :: i, idx, nvir, ne, fid
 character(len=240), intent(in) :: inpname
 character(len=240) :: buf, inporb

 i = index(inpname, '.input', back=.true.)
 inporb = inpname(1:i-1)//'.INPORB'

 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(1:4)=="&SCF" .or. buf(1:7)=="&RASSCF") exit
 end do ! for while

 BACKSPACE(fid)
 write(fid,'(/,A)') "&RASSCF"
 write(fid,'(A,I0)') 'Charge = ', charge
 write(fid,'(A,I0)') 'Spin = ', mult
 write(fid,'(A)') 'FILEORB = '//TRIM(inporb)

 idx = ndb + npair - npair0 ! doubly occupied orbitals in CAS
 nvir = nif - idx - 2*npair0 - nopen ! virtual orbitals in CAS
 ne = 2*idx + nacta + nactb ! all electrons

 if(CtrType == 1) then ! uncontracted MRCISD
  write(fid,'(A,I0,A)') 'nActEl = ', ne, ' 2 2'
  write(fid,'(A,I0)') 'RAS1 = ', idx
  write(fid,'(A,I0)') 'RAS2 = ', 2*npair0+nopen
  write(fid,'(A,I0)') 'RAS3 = ', nvir
  write(fid,'(A,/)') 'CIonly'
 else if(CtrType == 2) then ! ic-MRCISD
  write(fid,'(A,I0,A)') 'nActEl = ', nacta+nactb, ' 0 0'
  write(fid,'(A,I0)') 'RAS2 = ', 2*npair0+nopen
  write(fid,'(A)') 'CIonly'

  write(fid,'(/,A)') "&MOTRA"
  write(fid,'(A)') 'LUMORB'
  write(fid,'(A)') 'Frozen = 0'

  write(fid,'(/,A)') "&GUGA"
  write(fid,'(A,I0)') 'Spin = ', mult
  write(fid,'(A,I0)') 'Inactive = ', idx
  write(fid,'(A,I0)') 'Active = ', 2*npair0+nopen
  write(fid,'(A,I0)') 'NACTel = ', nacta+nactb
  write(fid,'(A,I0)') 'CIAll = 1'

  write(fid,'(/,A)') "&MRCI"
  write(fid,'(A,/)') 'MAXIterations = 49'
 end if

 close(fid)
 return
end subroutine prt_mrcisd_molcas_inp

! print MRCISD keywords into ORCA .inp file
subroutine prt_mrcisd_orca_inp(inpname1)
 use print_id, only: iout
 use mol, only: nopen, nacta, nactb, npair0, mult
 use mr_keyword, only: mem, nproc, CtrType, DKH2
 implicit none
 integer :: i, fid1, fid2
 integer :: RENAME
 character(len=240), intent(in) :: inpname1
 character(len=240) :: buf, inpname2

 inpname2 = TRIM(inpname1)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname2),status='replace')
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid2,'(A,I0,A)') '%maxcore ', CEILING(1000.0d0*DBLE(mem)/DBLE(nproc))
 write(fid2,'(A)') '! NoIter'

 if(DKH2) then
  write(fid2,'(A)') '%rel'
  write(fid2,'(A)') ' method DKH'
  write(fid2,'(A)') ' order 2'
  write(fid2,'(A)') 'end'
 end if

 if(CtrType == 1) then ! uncontracted MRCISD
  write(fid2,'(A)') '%mrci'
  write(fid2,'(2(A,I0),A)') ' NewBlock 1 * nroots 1 refs cas(',nacta+nactb,',',&
                             2*npair0+nopen,') end end'
  write(fid2,'(A)') ' tsel 0.0'
  write(fid2,'(A)') ' tpre 0.0'
  write(fid2,'(A)') ' Etol 1e-7'
  write(fid2,'(A)') ' Rtol 1e-7'
 else if(CtrType == 3) then ! FIC-MRCISD
  write(fid2,'(A)') '%autoci'
  write(fid2,'(A)') ' CItype FICMRCI'
  write(fid2,'(A,I0)') ' nel ',  nacta+nactb
  write(fid2,'(A,I0)') ' norb ', 2*npair0+nopen
  write(fid2,'(A,I0)') ' mult ', mult
  write(fid2,'(A)') ' nroots 1'
  write(fid2,'(A)') ' DavidsonOpt 1'
 end if

 write(fid2,'(A)') ' MaxIter 100'
 write(fid2,'(A)') 'end'
 write(fid2,'(A)') '%method'
 write(fid2,'(A)') ' FrozenCore FC_NONE'
 write(fid2,'(A)') 'end'
 write(fid2,'(A)') '%coords'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:6) == '%coord') exit
 end do ! for while

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname2), TRIM(inpname1))
 return
end subroutine prt_mrcisd_orca_inp

! print MRCISD keywords into Gaussian .gjf file
subroutine prt_mrcisd_gau_inp(gjfname)
 use mol, only: nif, ndb, nopen, nacta, nactb, npair, npair0
 use mr_keyword, only: mem, nproc, DKH2
 implicit none
 integer :: i, ne, nvir, idx, fid
 character(len=240), intent(in) :: gjfname

 idx = ndb + npair - npair0 ! doubly occupied orbitals in CAS
 nvir = nif - idx - 2*npair0 - nopen ! virtual orbitals in CAS
 ne = 2*idx + nacta + nactb ! all electrons

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 i = index(gjfname, '.gjf', back=.true.)
 write(fid,'(A,I0)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A5,I0,A2)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(4(A,I0),A,/)') '#p CAS(',ne,',',nif,',ras(2,',idx,',2,',nvir,'))/chkbasis&
                           & nosymm guess=read geom=allcheck scf(maxcycle=-1)'
 if(DKH2) then
  write(fid,'(A)') 'int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A)') 'int=nobasistransform'
 end if

 close(fid)
 return
end subroutine prt_mrcisd_gau_inp

! print MRCISD keywords into Molpro input file
subroutine prt_mrcisd_molpro_inp(inpname)
 implicit none
 integer :: i, fid
 character(len=240), intent(in) :: inpname

 call prt_cas_molpro_inp(inpname, .false., .false.)
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 BACKSPACE(fid)
 write(fid,'(A)') '{MRCIC;CORE}'
 close(fid)

 return
end subroutine prt_mrcisd_molpro_inp

! print MC-PDFT or DMRG-PDFT keywords into OpenMolcas .input file
subroutine prt_mcpdft_molcas_inp(inpname)
 use print_id, only: iout
 use mol, only: charge, mult, nacte, nacto
 use mr_keyword, only: CIonly, dmrgci, dmrgscf, maxM, otpdf, DKH2
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.tmp'

 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(2:7) == 'SEWARD') exit
 end do
 close(fid1,status='delete')

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_mcpdft_molcas_inp: no 'SEWARD'&
                   & found in file "//TRIM(inpname)
  stop
 end if

 write(fid2,'(A)') ' Grid input'
 write(fid2,'(A)') '  grid=ultrafine'
 write(fid2,'(A)') ' End of grid input'
 if(DKH2) write(fid2,'(A)') ' Relativistic = R02O'

 write(fid2,'(/,A)') "&RASSCF"
 write(fid2,'(A,I0)') 'Spin = ', mult
 write(fid2,'(A,I0)') 'Charge = ', charge
 write(fid2,'(A,I0)') 'nActEl= ', nacte
 write(fid2,'(A,I0)') 'RAS2 = ', nacto
 i = index(inpname, '.input', back=.true.)
 write(fid2,'(A)') 'FILEORB = '//inpname(1:i-1)//'.INPORB'

 if(dmrgci .or. dmrgscf) then
  write(fid2,'(A)') 'DMRG'
  write(fid2,'(A)') 'RGinput'
  write(fid2,'(A)') ' conv_thresh = 1E-7'
  write(fid2,'(A)') ' nsweeps = 5'
  write(fid2,'(A,I0)') ' max_bond_dimension = ', MaxM
  write(fid2,'(A)') 'endRG'
 end if

 if(CIonly) write(fid2,'(A)') 'CIonly'

 write(fid2,'(/,A)') "&MCPDFT"
 write(fid2,'(A)') 'KSDFT='//TRIM(otpdf)
 close(fid2)

 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine prt_mcpdft_molcas_inp

! print MC-PDFT keywords into GAMESS .inp file
subroutine prt_mcpdft_gms_inp(inpname)
 use print_id, only: iout
 use mol, only: charge, mult, ndb, nacte, nacto, npair, npair0
 use mr_keyword, only: mem, nproc, cart, otpdf, DKH2, hardwfn, crazywfn, CIonly
 implicit none
 integer :: i, ncore, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 ncore = ndb + npair - npair0
 inpname1 = TRIM(inpname)//'.tmp'

 open(newunit=fid2,file=TRIM(inpname1),status='replace')
 if(CIonly) then
  write(fid2,'(A)',advance='no') ' $CONTRL SCFTYP=NONE CITYP=ALDET'
 else
  write(fid2,'(A)',advance='no') ' $CONTRL SCFTYP=MCSCF'
 end if

 write(fid2,'(2(A,I0),A)') ' RUNTYP=ENERGY ICHARG=',charge,' MULT=',mult,' NOSYM=1 ICUT=11'
 write(fid2,'(A)',advance='no') '  PDFTYP='//TRIM(otpdf)

 if(DKH2) write(fid2,'(A)',advance='no') ' RELWFN=DK'
 if(.not. cart) then
  write(fid2,'(A)') ' ISPHER=1 $END'
 else
  write(fid2,'(A)') ' $END'
 end if

 ! MC-PDFT in GAMESS cannot run in parallel currently. All memory given to 1 core.
 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',mem*125,' $END'

 if(CIonly) then
  write(fid2,'(A)',advance='no') ' $CIDET'
 else
  write(fid2,'(A)',advance='no') ' $DET'
 end if
 write(fid2,'(3(A,I0))',advance='no') ' NCORE=',ncore,' NELS=',nacte,' NACT=',nacto

 if(hardwfn) then
  write(fid2,'(A)',advance='no') ' NSTATE=5'
 else if(crazywfn) then
  write(fid2,'(A)',advance='no') ' NSTATE=10'
 end if
 write(fid2,'(A)') ' ITERMX=500 $END'
 write(fid2,'(A)') ' $DFT NRAD=99 $END'

 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:6) == '$GUES') exit
 end do

 write(fid2,'(A)') TRIM(buf)
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')
 close(fid2)

 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine prt_mcpdft_gms_inp

! perform GVB computation (only in Strategy 1,3)
subroutine do_gvb()
 use print_id, only: iout
 use mr_keyword, only: nproc, gms_path, gms_scr_path, mo_rhf, ist, hf_fch, gvb,&
  datname, npair_wish, bgchg, chgname, cart, check_gms_path
 use mol, only: nbf, nif, ndb, nopen, npair, lin_dep, gvb_e, nacto, nacta, &
                nactb, nacte, npair0
 implicit none
 integer :: i, j, system, RENAME
 character(len=24) :: data_string = ' '
 character(len=240) :: buf, proname, inpname, gmsname, pair_fch
 character(len=300) :: longbuf = ' '

 if(.not. gvb) return
 write(iout,'(//,A)') 'Enter subroutine do_gvb...'
 call check_gms_path()

 if(.not. (ist==1 .or. ist==3)) then
  write(iout,'(A)') 'ERROR in subroutine do_gvb: current ist /= 1,3.'
  stop
 end if

 i = index(hf_fch, '.fch', back=.true.)
 proname = hf_fch(1:i-1)
 ! In RHF virtual MO projection, it will generate a file uno.out additionally
 call read_npair_from_uno_out(nbf, nif, ndb, npair, nopen, lin_dep)

 if(mo_rhf) then ! paired LMOs obtained from RHF virtual projection
  if(npair_wish>0 .and. npair_wish/=npair) then
   write(iout,'(A)') 'ERROR in subroutine do_gvb: npair_wish cannot be assigned in this strategy.'
   stop
  end if
  write(iout,'(A,I0,A1)') 'GVB(', npair, ')'

  pair_fch = TRIM(proname)//'_proj_loc_pair.fch'
  write(buf,'(2(A,I0))') 'fch2inp '//TRIM(pair_fch)//' -gvb ',npair,' -open ',nopen
  i = system(TRIM(buf))
  write(inpname,'(A,I0,A)') TRIM(proname)//'_proj_loc_pair2gvb',npair,'.inp'
  write(gmsname,'(A,I0,A)') TRIM(proname)//'_proj_loc_pair2gvb',npair,'.gms'
  write(datname,'(A,I0,A)') TRIM(proname)//'_proj_loc_pair2gvb',npair,'.dat'
  i = RENAME(TRIM(proname)//'_proj_loc_pair.inp', inpname)

 else ! paired LMOs obtained from associated rotation of UNOs

  if(npair_wish>0 .and. npair_wish/=npair) then
   write(iout,'(2(A,I0),A)') 'Warning: AutoMR recommends GVB(',npair,'), but&
    & user specifies GVB(',npair_wish,'). Try to fulfill...'
   if(npair_wish < npair) then
    ndb = ndb + npair - npair_wish
    npair = npair_wish
    write(iout,'(A)') 'OK, fulfilled.'
   else if(npair_wish > npair) then
    write(iout,'(A)') 'ERROR in subroutine do_gvb: too large space specified.&
                     & Cannot fulfilled.'
    stop
   end if
  end if

  write(iout,'(A,I0,A1)') 'GVB(', npair, ')'
  pair_fch = TRIM(proname)//'_uno_asrot.fch'
  write(buf,'(2(A,I0))') 'fch2inp '//TRIM(pair_fch)//' -gvb ', npair,' -open ',nopen
  i = system(TRIM(buf))
  write(inpname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.inp'
  write(gmsname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.gms'
  write(datname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.dat'
  i = RENAME(TRIM(proname)//'_uno_asrot.inp', inpname)
 end if

 call modify_memory_in_inp(inpname)

 ! call GAMESS to do GVB computations (delete .dat file first, if any)
 buf = TRIM(gms_scr_path)//'/'//TRIM(datname)
 call delete_file(buf)
 if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

 write(longbuf,'(A,I0,A)') TRIM(gms_path)//' '//TRIM(inpname)//' 01 ',nproc,' >& '//TRIM(gmsname)
 i = system(TRIM(longbuf))
 call read_gvb_energy_from_gms(gmsname, gvb_e)

 ! move the .dat file into current directory
 i = system('mv '//TRIM(gms_scr_path)//'/'//TRIM(datname)//' .')
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine do_gvb: fail to move file. Possibly wrong gms_scr_path.'
  stop
 end if

 ! sort the GVB pairs by CI coefficients of the 1st NOs
 if(cart) then ! Cartesian functions
  write(longbuf,'(A,5(1X,I0))') 'gvb_sort_pairs '//TRIM(datname), nbf, nif, ndb, nopen, npair
 else          ! spherical harmonic functions
  call read_nbf_from_dat(datname, i)
  write(longbuf,'(A,5(1X,I0))') 'gvb_sort_pairs '//TRIM(datname), i, nif, ndb, nopen, npair
 end if
 i = system(TRIM(longbuf))

 ! generate corresponding .fch file from _s.dat file
 i = index(datname, '.dat')
 inpname = datname(1:i-1)//'_s.fch'
 datname = datname(1:i-1)//'_s.dat'
 call copy_file(pair_fch, inpname, .false.)
 write(longbuf,'(2(A,I0))') 'dat2fch '//TRIM(datname)//' '//TRIM(inpname)//' -gvb ',&
                             npair, ' -open ', nopen
 i = system(TRIM(longbuf))
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine do_gvb: failed to call utility dat2fch.'
  write(iout,'(A)') 'Did you delete it or forget to compile it?'
  stop
 end if

 ! extract NOONs from the above .dat file and print them into .fch file
 write(longbuf,'(A,3(1X,I0),A5)') 'extract_noon2fch '//TRIM(datname)//' '//&
                     TRIM(inpname), ndb+1, ndb+nopen+2*npair, nopen, ' -gau'
 i = system(TRIM(longbuf))

 ! find npair0: the number of active pairs (|C2| > 0.1)
 call find_npair0_from_dat(datname, npair, npair0)

 ! determine the number of orbitals/electrons in following CAS/DMRG computations
 nacta = npair0 + nopen
 nactb = npair0
 nacte = nacta + nactb
 nacto = nacte

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_gvb at '//TRIM(data_string)
 return
end subroutine do_gvb

! do CASCI(npair<=7) or DMRG-CASCI(npair>7) when scf=.False.
! do CASSCF(npair<=7) or DMRG-CASSCF(npair>7) when scf=.True.
subroutine do_cas(scf)
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, casci, dmrgci, casscf, dmrgscf, ist, hf_fch,&
  datname, nacte_wish, nacto_wish, gvb, casnofch, casci_prog, casscf_prog, &
  dmrgci_prog, dmrgscf_prog, gau_path, gms_path, molcas_path, orca_path, &
  gms_scr_path, molpro_path, bdf_path, bgchg, chgname, casscf_force, dkh2_or_x2c,&
  check_gms_path, prt_strategy
 use mol, only: nbf, nif, npair, nopen, npair0, ndb, casci_e, casscf_e, nacta, &
                nactb, nacto, nacte, gvb_e, mult, ptchg_e, nuc_pt_e, natom, grad
 use util_wrapper, only: formchk, unfchk, gbw2mkl, mkl2gbw, fch2inp_wrap
 implicit none
 integer :: i, j, idx1, idx2, nvir, system, RENAME
 integer :: indb(2), inacta(2), inactb(2), inacto(2), inacte(2)
 real(kind=8) :: e(2)   ! e(1) is CASCI enery, e(2) is CASSCF energy
 character(len=10) :: cas_prog = ' '
 character(len=24) :: data_string = ' '
 character(len=240) :: buf, fchname, pyname, inpname, outname, proname, mklname
 character(len=240) :: orbname, xmlname
 logical, intent(in) :: scf
 logical :: rel

 if(scf) then
  if((.not. casscf) .and. (.not.dmrgscf)) return
 else
  if((.not. casci) .and. (.not.dmrgci)) return
 end if

 write(iout,'(//,A)') 'Enter subroutine do_cas...'
 if(ist == 5) then
  call read_no_info_from_fch(hf_fch,nbf,nif,indb,nopen,inacta,inactb,inacto,inacte)
  ! if ist=5, read nbf, nif, nopen, nacto, ... variables from NO .fch(k) file
  i = inacte(1); j = inacto(1)
 else
  i = 2*npair0 + nopen; j = i
 end if

 if(nacte_wish==0 .and. ist==5) then
  ndb = indb(1); nacta = inacta(1); nactb = inactb(1)
  nacto = inacto(1); nacte = inacte(2)
  npair0 = nactb; npair = npair0
 end if

 if(nacte_wish>0 .and. i/=nacte_wish) then
  write(iout,'(4(A,I0),A)') 'Warning: AutoMR recommends CAS(',i,'e,',j,'o), but&
   & user specifies CAS(',nacte_wish,'e,',nacto_wish, 'o). Trying to fulfill...'

  if(ist == 5) then
   if(nacte_wish > inacte(2)) then
    write(iout,'(A)') 'ERROR in subroutine do_cas: too large active space required.'
    write(iout,'(2(A,I0),A)') 'Maximum allowed: CAS(',inacte(2),',',inacte(2),')'
    stop
   else
    i = nacte_wish; j = nacto_wish
    npair0 = (nacte_wish-nopen)/2; npair = npair0
    ndb = indb(1) + inactb(1) - npair
    nacta = npair0 + nopen; nactb = npair0
    nacto = nacto_wish; nacte = nacto_wish
   end if
   write(iout,'(A)') 'OK, fulfilled.'

  else ! ist /= 5
   if(2*npair+nopen < nacte_wish) then
    write(iout,'(A)') 'ERROR in subroutine do_cas: too large space specified. Cannot fulfilled.'
    write(iout,'(2(A,I0))') '2*npair+nopen=', 2*npair+nopen, ', nacte_wish=', nacte_wish
    stop
   else ! 2*npair+nopen >= nacte_wish
    if(MOD(nacte_wish-nopen,2) /= 0) then
     write(iout,'(A)') 'ERROR in subroutine do_cas: wrong space specified. Cannot fulfilled.'
     write(iout,'(A)') 'nacte_wish-nopen is not an even integer.'
     write(iout,'(2(A,I0))') 'nopen=', nopen, ', nacte_wish=', nacte_wish
     stop
    end if

    write(iout,'(A)') 'OK, fulfilled.'
    npair0 = (nacte_wish-nopen)/2
    i = 2*npair0 + nopen
    nacta = npair0 + nopen; nactb = npair0
    nacto = nacto_wish; nacte = nacto_wish
   end if
  end if
 end if

 if(scf) then
  if(casscf) then
   data_string = 'CASSCF'
   cas_prog = casscf_prog
  else
   data_string = 'DMRG-CASSCF'
   cas_prog = dmrgscf_prog
  end if
 else ! scf = .False.
  if(casci) then
   data_string = 'CASCI'
   cas_prog = casci_prog
  else
   data_string = 'DMRG-CASCI'
   cas_prog = dmrgci_prog
  end if
 end if
 write(iout,'(A)',advance='no') TRIM(data_string)

 idx1 = ndb + npair - npair0 + 1
 idx2 = idx1 + 2*npair0 + nopen - 1
 nvir = nif - (idx1-1) - 2*npair0 - nopen
 write(iout,'(A,2(I0,A))') '(',i,',',i,') using program '//TRIM(cas_prog)

 if(i == 0) then
  write(iout,'(/,A)') 'There is no active orbital/electron. AutoMR terminated.'
  write(iout,'(A)') 'The reason is this molecule has no multi-configurational&
                   & or multi-reference character.'
  write(iout,'(A)') 'This molecule can be well described by single reference&
                   & methods, e.g. MP2, CCSD(T).'
  write(iout,'(A)') 'Thus no need for multi-reference computation. But if you&
                   & have to do it, you can manually'
  write(iout,'(A)') 'specify the size of acitve space in .gjf file. For example&
                   & CASSCF(8,8) for methane(CH4).'
  write(iout,'(A)') 'The maximum of active orbitals is 2*npair. You can find&
                   & npair in GVB computations above.'
  stop
 end if
 write(iout,'(A,I4,4X,A,I4)') 'doubly_occ=', idx1-1, 'nvir=', nvir
 write(iout,'(2(A,I0))') 'No. of active alpha/beta e = ', nacta,'/',nactb

 if(nopen+2*npair0 > 15) then
  if(scf) then
   casscf = .false.
   dmrgscf = .true.
   cas_prog = dmrgscf_prog
   write(iout,'(A)') 'Warning: CASSCF is switched to DMRG-CASSCF due to active&
                    & space larger than (15,15).'
   if(casscf_prog /= 'pyscf') then
    write(iout,'(A)') 'ERROR in subroutine do_cas: DMRGSCF required. But&
                     & CASSCF_prog='//TRIM(casscf_prog)//'.'
    stop
   end if

  else ! not scf
   casci = .false.
   dmrgci = .true.
   cas_prog = dmrgci_prog
   write(iout,'(A)') 'Warning: CASCI is switched to DMRG-CASCI due to active&
                    & space larger than (15,15).'
   if(casci_prog /= 'pyscf') then
    write(iout,'(A)') 'ERROR in subroutine do_cas: DMRGCI required. But&
                     & CASCI_prog='//TRIM(casci_prog)//'.'
    stop
   end if
  end if

  write(iout,'(A)') 'Strategy updated:'
  call prt_strategy()
 end if

 if(ist<1 .or. ist>5) then
  write(iout,'(A)') 'ERROR in subroutine do_cas: ist out of range.'
  write(iout,'(A,I0)') 'ist=', ist
  stop
 end if

 if((dmrgci .or. dmrgscf) .and. cas_prog/='pyscf') then
  write(iout,'(A)') 'ERROR in subroutine do_cas: DMRG-CASCI/CASSCF calculation&
                   & is only supported by PySCF.'
  write(iout,'(A)') 'Wrong casci_prog or casscf_prog: '//TRIM(cas_prog)
  stop
 end if

 if(ist==1 .or. ist==3) then
  i = index(datname, '.dat', back=.true.)
  fchname = datname(1:i-1)//'.fch'
  pyname = datname(1:i-1)//'.py'
 else if(ist == 2) then ! UHF -> UNO -> CASCI/CASSCF
  i = index(hf_fch, '.fch', back=.true.)
  fchname = hf_fch(1:i-1)//'_uno.fch'
  pyname = hf_fch(1:i-1)//'_uno.py'
  inpname = hf_fch(1:i-1)//'_uno.py2'
  i = RENAME(TRIM(pyname), TRIM(inpname))
  ! bas_fch2py will generate file '_uno.py', so we need to rename it to another filename
 else if(ist == 5) then
  i = index(hf_fch, '.fch', back=.true.)
  fchname = hf_fch
  pyname = hf_fch(1:i-1)//'.py'
 end if

 proname = ' '
 i = index(hf_fch, '.fch', back=.true.)
 if(ist==1 .or. ist==3) then
  if(scf) then
   write(proname,'(A,I0,A)') hf_fch(1:i-1)//'_gvb', npair, '_2CASSCF'
  else
   write(proname,'(A,I0,A)') hf_fch(1:i-1)//'_gvb', npair, '_2CASCI'
  end if
 else if(ist == 2) then
  if(scf) then
   write(proname,'(A)') hf_fch(1:i-1)//'_uno2CASSCF'
  else
   write(proname,'(A)') hf_fch(1:i-1)//'_uno2CASCI'
  end if
 end if
 casnofch = TRIM(proname)//'_NO.fch'

 select case(TRIM(cas_prog))
 case('pyscf')
  i = system('bas_fch2py '//TRIM(fchname))
  inpname = TRIM(proname)//'.py'
  i = RENAME(TRIM(pyname), TRIM(inpname))
  call prt_cas_script_into_py(inpname, fchname, scf)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  if(casscf_force) then
   i = system("echo 'from pyscf import grad' >> "//TRIM(inpname))
   i = system("echo 'mc.Gradients().kernel()' >> "//TRIM(inpname))
  end if
  j = index(inpname, '.py', back=.true.)
  outname = inpname(1:j-1)//'.out'
  i = system('python '//TRIM(inpname)//' >& '//TRIM(outname))

  write(buf,'(A,2(1X,I0))') 'extract_noon2fch '//TRIM(outname)//' '//&
                             TRIM(casnofch), idx1, idx2
  i = system(TRIM(buf))

 case('gaussian')
  call check_exe_exist(gau_path)

  inpname = TRIM(proname)//'.gjf'
  outname = TRIM(proname)//'.log'
  mklname = TRIM(proname)//'.chk'
  call prt_cas_gjf(inpname, nacto, nacte, scf, casscf_force)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call unfchk(fchname, mklname)

  i = system(TRIM(gau_path)//' '//TRIM(inpname))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_cas: Gaussian CAS job failed!'
   write(iout,'(A)') 'Filename='//TRIM(inpname)
   stop
  end if
  call formchk(mklname, casnofch)

 case('gamess')
  call check_gms_path()
  inpname = TRIM(proname)//'.dat'
  buf = TRIM(gms_scr_path)//'/'//TRIM(inpname) ! delete the possible .dat file
  call delete_file(buf)
  ! do not use datname in the above three lines! because datname may be that of a GVB job

  call fch2inp_wrap(fchname, .false., .false., 0, 0)

  i = index(fchname, '.fch', back=.true.)
  outname = fchname(1:i-1)//'.inp'
  inpname = TRIM(proname)//'.inp'
  i = RENAME(TRIM(outname), TRIM(inpname))
  outname = TRIM(proname)//'.gms'
  call prt_cas_gms_inp(inpname, idx1-1, scf)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  if(casscf_force) i = system("sed -i '1,1s/ENERGY/GRADIENT/' "//TRIM(inpname))
  write(buf,'(A,I0,A)') TRIM(gms_path)//' '//TRIM(inpname)//' 01 ', nproc, " >&"//TRIM(outname)
  i = system(TRIM(buf))

  ! make a copy of the .fch file to save NOs
  if(ist /= 2) then ! datname is a GVB job .dat file
   i = index(datname,'.dat')
   inpname = datname(1:i-1)//'.fch'
  else              ! no GVB job
   i = index(hf_fch, '.fch')
   inpname = hf_fch(1:i-1)//'_uno.fch'
  end if
  call copy_file(inpname, casnofch, .false.)

  ! move the .dat file into current directory
  datname = TRIM(proname)//'.dat'   ! update datname
  i = system('mv '//TRIM(gms_scr_path)//'/'//TRIM(datname)//' .')

  ! transfer CASSCF pseudo-canonical MOs from .dat to .fch
  if(scf) then
   i = system('dat2fch '//TRIM(datname)//' '//TRIM(casnofch))
   write(buf,'(A,2(1X,I0))') 'extract_noon2fch '//TRIM(outname)//' '//&
                              TRIM(casnofch), idx1, idx2
   i = system(TRIM(buf))
  end if

  ! transfer NOs from .dat to .fch
  write(buf,'(A,I0,1X,I0)') 'dat2fch '//TRIM(datname)//' '//TRIM(casnofch)//' -no ',1,idx2
  i = system(TRIM(buf))

 case('openmolcas')
  call check_exe_exist(molcas_path)

  i = system('fch2inporb '//TRIM(fchname))
  i = index(fchname, '.fch', back=.true.)
  outname = fchname(1:i-1)//'.INPORB'
  inpname = TRIM(proname)//'.INPORB'
  i = RENAME(TRIM(outname),TRIM(inpname))
  i = index(fchname, '.fch', back=.true.)
  outname = fchname(1:i-1)//'.input'
  inpname = TRIM(proname)//'.input'
  i = RENAME(TRIM(outname),TRIM(inpname))
  outname = TRIM(proname)//'.out'
  orbname = TRIM(proname)//'.RasOrb.1'
  call prt_cas_molcas_inp(inpname, scf)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  if(casscf_force) i = system("echo '&ALASKA' >> "//TRIM(inpname))
  i = system(TRIM(molcas_path)//' '//TRIM(inpname)//" >& "//TRIM(outname))

  ! make a copy of the .fch file to save NOs
  call copy_file(fchname, casnofch, .false.)

  ! transfer NOs from .dat to .fch
  i = system('orb2fch '//TRIM(orbname)//' '//TRIM(casnofch)//' -no')

 case('orca')
  call check_exe_exist(orca_path)
  i = system('fch2mkl '//TRIM(fchname))

  i = index(fchname, '.fch', back=.true.)
  pyname = fchname(1:i-1)//'.inp'
  orbname = fchname(1:i-1)//'.mkl'
  inpname = TRIM(proname)//'.inp'
  mklname = TRIM(proname)//'.mkl'
  outname = TRIM(proname)//'.out'
  i = RENAME(TRIM(orbname),TRIM(mklname))
  i = RENAME(TRIM(pyname),TRIM(inpname))

  call prt_cas_orca_inp(inpname, scf)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(mklname)
  if(casscf_force) i = system("sed -i '3,3s/TightSCF/TightSCF EnGrad/' "//TRIM(inpname))
  i = system(TRIM(orca_path)//' '//TRIM(inpname)//" >& "//TRIM(outname))

  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  if(scf) then ! CASSCF
   orbname = TRIM(proname)//'.gbw'
  else         ! CASCI
   pyname = TRIM(proname)//'.b0_s0.nat'
   orbname = TRIM(proname)//'.gbw'
   i = RENAME(TRIM(pyname),TRIM(orbname))
  end if

  call gbw2mkl(orbname)
  mklname = TRIM(proname)//'.mkl'
  i = system('mkl2fch '//TRIM(mklname)//' '//TRIM(casnofch)//' -no')

 case('molpro')
  call check_exe_exist(molpro_path)

  i = system('fch2com '//TRIM(fchname)) ! generate .com and .txt
  i = index(fchname, '.fch', back=.true.)
  mklname = fchname(1:i-1)//'.com'
  pyname = fchname
  call convert2molpro_fname(pyname, '.a')
  inpname = TRIM(proname)//'.com'
  orbname = TRIM(proname)//'.fch'
  call convert2molpro_fname(orbname, '.a')
  outname = TRIM(proname)//'.out'
  xmlname = TRIM(proname)//'.xml'
  i = RENAME(TRIM(mklname), TRIM(inpname))
  i = RENAME(TRIM(pyname), TRIM(orbname))
  call prt_cas_molpro_inp(inpname, scf, casscf_force)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  i = CEILING(DBLE(mem*125)/DBLE(nproc))
  write(buf,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -m ', i,&
                           'm '//TRIM(inpname)
  i = system(TRIM(buf))
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine do_cas: Molpro CASCI/CASSCF job failed.'
   stop
  end if
  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  i = system('xml2fch '//TRIM(xmlname)//' '//TRIM(casnofch)//' -no')

 case('bdf')
  call check_exe_exist(bdf_path)

  i = system('fch2bdf '//TRIM(fchname)//' -no') ! generate _bdf.inp, .BAS, .inporb
  i = index(fchname, '.fch', back=.true.)
  mklname = fchname(1:i-1)//'_bdf.inp'
  pyname  = fchname(1:i-1)//'_bdf.inporb'
  inpname = TRIM(proname)//'.inp'
  xmlname = TRIM(proname)//'.inporb'
  orbname = TRIM(proname)//'.casorb'
  outname = TRIM(proname)//'.out'
  i = RENAME(TRIM(mklname), TRIM(inpname))
  i = RENAME(TRIM(pyname), TRIM(xmlname))
  call prt_cas_bdf_inp(inpname, scf, casscf_force)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  i = system(TRIM(bdf_path)//' '//TRIM(proname))
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine do_cas: BDF CASCI/CASSCF job failed.'
   stop
  end if
  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  i = system('bdf2fch '//TRIM(orbname)//' '//TRIM(casnofch)//' -no')

 case default
  write(iout,'(A)') 'ERROR in subroutine do_cas: invalid CAS_prog. Allowed&
                   & programs are Gaussian, GAMESS, PySCF, ORCA, Molpro.'
  write(iout,'(A)') 'CAS_prog='//TRIM(cas_prog)
  stop
 end select

 ! i is 'extract NOONs from the output file and print them into .fch file'
 if(i /= 0) then
  write(iout,'(A)') 'Warning in subroutine do_cas: possibly extract_noon2fch failed.&
                   & Filename = '//TRIM(outname)
  write(iout,'(A)') 'This does not affect the CASCI/CASSCF energy. So continue.'
 end if

 if(ist == 2) then
  i = index(hf_fch, '.fch', back=.true.)
  pyname = hf_fch(1:i-1)//'_uno.py'
  inpname = hf_fch(1:i-1)//'_uno.py2'
  i = RENAME(TRIM(inpname), TRIM(pyname)) ! rename it back
 end if

 ! read energy, check convergence and check spin
 call read_cas_energy_from_output(cas_prog, outname, e, scf, nacta-nactb, &
                                  (dmrgci.or.dmrgscf), ptchg_e, nuc_pt_e)

 if(gvb .and. 2*npair+nopen==nacto .and. e(1)-gvb_e>2D-6) then
  write(iout,'(A)') 'ERROR in subroutine do_cas: active space of GVB and CAS&
                   & are equal, but CASCI/CASSCF energy'
  write(iout,'(A)') 'is higher than that of GVB. This is probably due to: (1)&
                   & CASCI stucks in a higher energy'
  write(iout,'(A)') 'local minimum or not pure spin state; (2) GVB MOs are di&
                   &sordered (GAMESS bug).'
  stop
 end if

 casci_e = e(1)
 write(iout,'(/,A,F18.8,1X,A4)') 'E(CASCI)  = ', e(1), 'a.u.'
 if(scf) then
  casscf_e = e(2)
  write(iout,'(A,F18.8,1X,A4)') 'E(CASSCF) = ', e(2), 'a.u.'
 end if

 if(casscf_force) then
  allocate(grad(3*natom))

  select case(cas_prog)
  case('pyscf')
   call read_grad_from_pyscf_out(outname, natom, grad)
  case('gaussian')
   call read_grad_from_gau_log(outname, natom, grad)
  case('gamess')
   call read_grad_from_gms_dat(datname, natom, grad)
  case('openmolcas')
   call read_grad_from_molcas_out(outname, natom, grad)
  case('orca')
   call read_grad_from_orca_out(outname, natom, grad)
  case('molpro')
   call read_grad_from_molpro_out(outname, natom, grad)
  case('bdf')
   call read_grad_from_bdf_out(outname, natom, grad)
  case default
   write(iout,'(A)') 'ERROR in subroutine do_cas: program cannot be identified.'
   write(iout,'(A)') 'cas_prog='//TRIM(cas_prog)
   stop
  end select

  write(iout,'(A)') 'Cartesian gradient (HARTREE/BOHR):'
  write(iout,'(5(1X,ES15.8))') (grad(i),i=1,3*natom)
 end if

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_cas at '//TRIM(data_string)
 return
end subroutine do_cas

! do uncontracted/ic-/FIC- MRCISD(+Q) for npair<=7, or <=CAS(14,14)
subroutine do_mrcisd()
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, casci, casscf, CIonly, ist, hf_fch, mrcisd,&
  mrcisd_prog, CtrType, casnofch, molcas_path, orca_path, gau_path, molpro_path,&
  bgchg, casci_prog, casscf_prog, chgname
 use mol, only: nbf, nif, npair, nopen, npair0, ndb, casci_e, casscf_e, davidson_e,&
  mrcisd_e, ptchg_e
 use util_wrapper, only: unfchk, mkl2gbw
 integer :: i, system
 real(kind=8) :: e
 character(len=24) :: data_string
 character(len=240) :: string, chkname, inpname, outname, mklname
 character(len=47) :: error_warn='ERROR in subroutine do_mrcisd: invalid CtrType='

 if(.not. mrcisd) return
 write(iout,'(//,A)') 'Enter subroutine do_mrcisd...'

 if(npair0 > 7) then
  write(iout,'(A)') 'ERROR in subroutine do_mrcisd: reference wfn larger than (14,14).'
  write(iout,'(A)') 'DMRG-MRCISD is not supported currently.'
  stop
 end if

 if(.not. CIonly) then
  if(TRIM(casscf_prog) == 'orca') then
   write(iout,'(A)') 'Warning: ORCA is used as the CASSCF solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(iout,'(A)') 'This will affect the CI energy up to 10^-5 a.u.&
                    & Such small error is usually not important.'
   write(iout,'(A)') 'If you care about the accuracy, please use another CASSCF solver.'
  end if

  string = 'MRCISD based on optimized CASSCF orbitals.'
 else ! CIonly = .True.

  if(TRIM(casci_prog) == 'orca') then
   write(iout,'(A)') 'Warning: ORCA is used as the CASCI solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(iout,'(A)') 'This will affect the CI energy up to 10^-5 a.u.&
                    & Such small error is usually not important.'
   write(iout,'(A)') 'If you care about the accuracy, please use another CASCI solver.'
  end if

  string = 'MRCISD based on CASCI orbitals.'
  write(iout,'(A)') 'Warning: the CASSCF orbital optimization is strongly recommended'
  write(iout,'(A)') 'to be performed before MRCI, unless it is too time-consuming.'
 end if
 write(iout,'(A)') TRIM(string)

 write(iout,'(A)') 'Frozen_Core = F, MRCISD computation using program '//TRIM(mrcisd_prog)

 select case(TRIM(mrcisd_prog))
 case('openmolcas')
  call check_exe_exist(molcas_path)
  i = system('fch2inporb '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  chkname = casnofch(1:i-1)//'.INPORB'
  string  = casnofch(1:i-1)//'.input'
  i = index(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'MRCISD.INPORB'
  inpname = casnofch(1:i)//'MRCISD.input'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(chkname), TRIM(mklname))
  i = RENAME(TRIM(string), TRIM(inpname))
  chkname = ' '
  call prt_mrcisd_molcas_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  i = system(TRIM(molcas_path)//' '//TRIM(inpname)//" >& "//TRIM(outname))

 case('orca')
  call check_exe_exist(orca_path)
  i = system('fch2mkl '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  chkname = casnofch(1:i-1)//'.mkl'
  string  = casnofch(1:i-1)//'.inp'
  i = index(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'MRCISD.mkl'
  inpname = casnofch(1:i)//'MRCISD.inp'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(chkname), TRIM(mklname))
  i = RENAME(TRIM(string), TRIM(inpname))
  chkname = ' '
  call prt_mrcisd_orca_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(mklname)
  i = system(TRIM(orca_path)//' '//TRIM(inpname)//" >& "//TRIM(outname))

 case('gaussian')
  call check_exe_exist(gau_path)
  i = index(casnofch, '_NO', back=.true.)
  chkname = casnofch(1:i)//'MRCISD.chk'
  inpname = casnofch(1:i)//'MRCISD.gjf'
  outname = casnofch(1:i)//'MRCISD.log'
  call unfchk(casnofch, chkname)
  call prt_mrcisd_gau_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  i = system(TRIM(gau_path)//' '//TRIM(inpname))

 case('molpro')
  call check_exe_exist(molpro_path)

  i = system('fch2com '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'.com'
  chkname = casnofch
  call convert2molpro_fname(chkname, '.a')

  i = index(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.com'
  outname = casnofch(1:i)//'MRCISD.out'
  mklname = inpname
  call convert2molpro_fname(mklname, '.a')
  i = RENAME(TRIM(string), TRIM(inpname))
  i = RENAME(TRIM(chkname), TRIM(mklname))

  call prt_mrcisd_molpro_inp(inpname)
  i = CEILING(DBLE(mem*125)/DBLE(nproc))
  write(string,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -m ', i,'m '//TRIM(inpname)
  i = system(TRIM(string))
 end select

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine do_mrcisd: error termination.'
  write(iout,'(A)') 'Filename='//TRIM(outname)
  stop
 end if

 ! read Davidson correction and MRCISD energy from OpenMolcas/ORCA/Gaussian output file
 call read_mrcisd_energy_from_output(CtrType, mrcisd_prog, outname, ptchg_e, davidson_e, e)
 mrcisd_e = e + davidson_e ! E(MRCISD+Q)
 if(CIonly) then   ! E(MRCISD) - (E(CASCI) or E(CASSCF))
  e = e - casci_e
 else
  e = e - casscf_e
 end if

 select case(TRIM(mrcisd_prog))
 case('openmolcas')
  select case(CtrType)
  case(1) ! uncontracted MRCISD
   write(iout,'(/,A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E(MRCISD)      =', mrcisd_e, 'a.u.'
  case(2) ! ic-MRCISD
   write(iout,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E_corr(icMRCISD) =', e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E(icMRCISD+Q)    =', mrcisd_e, 'a.u.'
  case default
   write(iout,'(A,I0)') error_warn, CtrType
   stop
  end select
 case('orca')
  write(iout,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
  select case(CtrType)
  case(1) ! uncontracted MRCISD
   write(iout,'(A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E(MRCISD+Q)    =', mrcisd_e, 'a.u.'
  case(3) ! FIC-MRCISD
   write(iout,'(A,F18.8,1X,A4)') 'E_corr(FIC-MRCISD) =', e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E(FIC-MRCISD+Q)    =', mrcisd_e, 'a.u.'
  case default
   write(iout,'(A,I0)') error_warn, CtrType
   stop
  end select
 case('gaussian') ! only uncontracted MRCISD
  if(CtrType /= 1) then
   write(iout,'(A,I0)') error_warn, CtrType
   stop
  end if
  write(iout,'(/,A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(MRCISD)      =', mrcisd_e, 'a.u.'
 case('molpro')
  write(iout,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(MRCISD)      =', mrcisd_e, 'a.u.'
 end select

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_mrcisd at '//TRIM(data_string)
 return
end subroutine do_mrcisd

! do ic-MRCC based on CASSCF, npair<=5
subroutine do_mrcc()
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, casci, casscf, CIonly
 implicit none
! integer :: i
! character(len=24) :: data_string

! if(.not. mrcc) return
! write(iout,'(//,A)') 'Enter subroutine do_mrcc...'

! write(iout,'(/,A)') 'Enter subroutine do_mrcc...'
! call fdate(data_string)
! write(iout,'(A)') 'ERROR in subroutine do_mrcc: not supported currently.'
 return
end subroutine do_mrcc

! modify the memory in a given .inp file
subroutine modify_memory_in_inp(inpname)
 use print_id, only: iout
 use mr_keyword, only: mem, nproc
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf,'MWORDS') /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine modify_memory_in_inp: no 'MWORDS' found&
                   & in file "//TRIM(inpname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',FLOOR(DBLE(mem)*1000.0d0/(8.0d0*DBLE(nproc))),' $END'

 ! copy the remaining content
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine modify_memory_in_inp

