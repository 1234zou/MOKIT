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
  write(iout,'(A,/)') " Example 2 (in dash): automr a.gjf >a.out 2>&1 &"
  stop
 end if

 call getarg(1, fname)

 select case(TRIM(fname))
 case('-v', '-V', '--version')
  write(iout,'(A)') 'AutoMR 1.2.2 :: MOKIT'
  stop
 case('-h','-help','--help')
  write(iout,'(/,A)')  "Usage: automr [gjfname] >& [outname]"
  write(iout,'(A)')    "  Example 1 (in bash): automr a.gjf >& a.out &"
  write(iout,'(A)')    "  Example 2 (in dash): automr a.gjf >a.out 2>&1 &"
  write(iout,'(/,A)')  'Options:'
  write(iout,'(A)')    '  -h, -help, --help: Print this message and exit.'
  write(iout,'(A)')    '  -v, -V, --version: Print the version number of automr and exit.'
  write(iout,'(/,A)')  'Methods(#p ...):'
  write(iout,'(A)')    '  GVB, CASCI, CASSCF, NEVPT2, NEVPT3, CASPT2, CASPT3, MRMP2, MRCISD,'
  write(iout,'(A)')    '  MCPDFT, SDSPT2, DMRGCI, DMRGSCF'
  write(iout,'(/,A)')  'Keywords in MOKIT{}:'
  write(iout,'(A)')    '  CASSCF_prog=PySCF, OpenMolcas, ORCA, Molpro, GAMESS, Gaussian, BDF, PSI4'
  write(iout,'(A)')    '   CASCI_prog=PySCF, OpenMolcas, ORCA, Molpro, GAMESS, Gaussian, BDF, PSI4'
  write(iout,'(A)')    '  NEVPT2_prog=PySCF, OpenMolcas, ORCA, Molpro, BDF'
  write(iout,'(A)')    '  CASPT2_prog=OpenMolcas, Molpro'
  write(iout,'(A)')    '  MCPDFT_prog=OpenMolcas, GAMESS'
  write(iout,'(A)')    '  MRCISD_prog=OpenMolcas, Molpro, ORCA, Gaussian, PSI4'
  write(iout,'(A,/)')  '  CtrType=1/2/3 for uc-/ic-/FIC- MRCISD'
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
 use mol, only: natom, atom2frag, nfrag, frag_char_mult, coor, elem, nuc, charge,&
  mult, rhf_e, uhf_e
 use mr_keyword, only: readuhf, readrhf, skiphf, mem, nproc, basis, cart, &
  gau_path, hf_fch, ist, mo_rhf, bgchg, read_bgchg_from_gjf, gjfname, chgname,&
  uno, dkh2_or_x2c, vir_proj, prt_strategy, gau_path, frag_guess
 use util_wrapper, only: formchk
 implicit none
 integer :: i, system
 real(kind=8) :: ssquare = 0d0
 character(len=24) :: data_string = ' '
 character(len=240) :: rhf_gjfname, uhf_gjfname, chkname
 logical :: eq

 write(iout,'(//,A)') 'Enter subroutine do_hf...'

 if(skiphf) then
  write(iout,'(A)') 'Provided .fch(k) file. Skip the RHF/UHF step...'
  call read_natom_from_fch(hf_fch, natom)
  allocate(coor(3,natom), elem(natom), nuc(natom))
  call read_elem_and_coor_from_fch(hf_fch, natom, elem, nuc, coor, charge, mult)

  if(readuhf .and. mult==1) then
   write(iout,'(A)') 'Check whether provided UHF is equivalent to RHF...'
   call check_if_uhf_equal_rhf(hf_fch, eq)
   if(eq) then
    write(iout,'(A)') 'This is actually a RHF wave function. Alpha=Beta.&
                     & Switching to ist=3.'
    i = system('fch_u2r '//TRIM(hf_fch))
    i = index(hf_fch, '.fch', back=.true.)
    hf_fch = hf_fch(1:i-1)//'_r.fch'
    readuhf = .false.; readrhf = .true.; ist = 3
    vir_proj = .true.; mo_rhf = .true. ; uno = .false.
    write(iout,'(A)') 'Strategy updated:'
    call prt_strategy()
   else
    write(iout,'(A)') 'This seems a truly UHF wave function.'
   end if
  end if

  if(bgchg) call read_bgchg_from_gjf(.true.)
  call fdate(data_string)
  write(iout,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)
  return
 end if

 call check_exe_exist(gau_path)
 call read_natom_from_gjf(gjfname, natom)
 allocate(coor(3,natom), elem(natom), nuc(natom))
 call read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
 if(frag_guess) then
  allocate(atom2frag(natom), frag_char_mult(2,nfrag))
  call read_frag_guess_from_gjf(gjfname, natom, atom2frag, nfrag, frag_char_mult)
 end if
 if(bgchg) call read_bgchg_from_gjf(.false.)

 i = index(gjfname, '.gjf', back=.true.)
 rhf_gjfname = gjfname(1:i-1)//'_rhf.gjf'
 uhf_gjfname = gjfname(1:i-1)//'_uhf.gjf'

 if(mult == 1) then ! singlet, perform RHF and UHF
  call generate_hf_gjf(rhf_gjfname, .false.)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(rhf_gjfname))
  call perform_scf_and_read_e(gau_path, rhf_gjfname, rhf_e, ssquare)

  call generate_hf_gjf(uhf_gjfname, .true.)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_gjfname))
  call perform_scf_and_read_e(gau_path, uhf_gjfname, uhf_e, ssquare)

 else               ! not singlet, only perform UHF

  call generate_hf_gjf(uhf_gjfname, .true.)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_gjfname))
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
  write(iout,'(A)') 'RHF/UHF energy is equal, or has little difference, choose RHF.'
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
 character(len=240) :: buf, proname, pyname, chkname, outname, fchname

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
  write(buf,'(A)') 'python '//TRIM(pyname)//' >'//TRIM(outname)//" 2>&1"
  write(iout,'(A)') '$'//TRIM(buf)
  i = system(TRIM(buf))

  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine get_paired_LMO: PySCF job fails.'
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
 write(fid2,'(A)') "py2fch('"//TRIM(proj_fch)//"',nbf,nif,mf.mo_coeff,Sdiag,'a',mf.mo_occ,False)"
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
  write(fid2,'(A)') "occ_loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],&
                   &mol._bas[:,3],mol.cart,nbf,npair,mf.mo_coeff[:,occ_idx],&
                   &S,'mulliken')"
  write(fid2,'(A)') "vir_loc_orb = pm(mol.nbas,mol._bas[:,0],mol._bas[:,1],&
                   &mol._bas[:,3],mol.cart,nbf,npair,mf.mo_coeff[:,vir_idx],&
                   &S,'mulliken')"
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
 write(fid2,'(A)') 'alpha_coeff = pair_by_tdm(ncore, npair, nopen, nalpha, nvir_lmo,&
                   &nbf, nif, mf.mo_coeff, mo_dipole)'
 write(fid2,'(A)') 'mf.mo_coeff = alpha_coeff.copy()'
 write(fid2,'(A)') '# pair done'

 write(fid2,'(/,A)') '# save the paired LMO into .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(hf_fch)//"', '"//TRIM(loc_fch)//"')"
 write(fid2,'(A)') "py2fch('"//TRIM(loc_fch)//"',nbf,nif,mf.mo_coeff,Sdiag,'a',mf.mo_occ,False)"
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
 write(fid1,'(A)') 'idx, noon, alpha_coeff = uno(nbf,nif,na,nb,mf.mo_coeff[0],mf.mo_coeff[1],S)'
 write(fid1,'(A)') 'alpha_coeff = construct_vir(nbf, nif, idx[1], alpha_coeff, S)'
 write(fid1,'(A)') 'mf.mo_coeff = (alpha_coeff, beta_coeff)'
 write(fid1,'(A)') '# done transform'

 write(fid1,'(/,A)') '# save the UNO into .fch file'
 write(fid1,'(A)') "os.system('fch_u2r "//TRIM(hf_fch)//"')"
 write(fid1,'(A)') "os.rename('"//hf_fch(1:i-1)//"_r.fch', '"//TRIM(uno_fch)//"')"
 write(fid1,'(A)') "py2fch('"//TRIM(uno_fch)//"',nbf,nif,mf.mo_coeff[0],Sdiag,'a',noon,True)"
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
 write(fid1,'(A)') "py2fch('"//TRIM(assoc_fch)//"',nbf,nif,mf.mo_coeff[0],Sdiag,'a',noon,False)"
 close(fid1)
 return
end subroutine prt_assoc_rot_script_into_py

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

! print MRCISD keywords into a PSI4 input file
subroutine prt_mrcisd_psi4_inp(inpname)
 use mol, only: nif, ndb, nopen, nacto, npair, npair0
 use mr_keyword, only: iout, mem
 implicit none
 integer :: i, fid, idx, nvir
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 idx = ndb + npair - npair0 ! doubly occupied orbitals in CAS
 nvir = nif - idx - 2*npair0 - nopen ! virtual orbitals in CAS

 call modify_memory_in_psi4_inp(inpname, mem)
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(index(buf,'d_convergence') > 0) exit

  if(buf(1:6) == 'memory') then
   write(iout,'(A)') 'ERROR in subroutine prt_mrcisd_psi4_inp: incomplete&
                    & file '//TRIM(inpname)
   close(fid)
   stop
  end if
 end do ! for while

 write(fid,'(A)') ' ci_maxiter 100'
 write(fid,'(A,I0,A)') ' ras1 [',idx  ,']'
 write(fid,'(A,I0,A)') ' ras2 [',nacto,']'
 write(fid,'(A,I0,A)') ' ras3 [',nvir ,']'
 write(fid,'(A)') '}'

 write(fid,'(/,A)') "mrcisd_energy = energy('detci',ref_wfn=scf_wfn)"
 close(fid)
 return
end subroutine prt_mrcisd_psi4_inp

! perform GVB computation (only in Strategy 1,3)
subroutine do_gvb()
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, gms_path, gms_scr_path, mo_rhf, ist, hf_fch,&
  gvb, datname, npair_wish, bgchg, chgname, cart, check_gms_path
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

  pair_fch = TRIM(proname)//'_proj_loc_pair.fch'
  write(buf,'(2(A,I0))') 'fch2inp '//TRIM(pair_fch)//' -gvb ',npair,' -open ',nopen
  write(iout,'(A)') '$'//TRIM(buf)
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

  pair_fch = TRIM(proname)//'_uno_asrot.fch'
  write(buf,'(2(A,I0))') 'fch2inp '//TRIM(pair_fch)//' -gvb ',npair,' -open ',nopen
  write(iout,'(A)') '$'//TRIM(buf)
  i = system(TRIM(buf))
  write(inpname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.inp'
  write(gmsname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.gms'
  write(datname,'(A,I0,A)') TRIM(proname)//'_uno_asrot2gvb',npair,'.dat'
  i = RENAME(TRIM(proname)//'_uno_asrot.inp', inpname)
 end if

 call modify_memory_in_gms_inp(inpname, mem, nproc)

 ! call GAMESS to do GVB computations (delete .dat file first, if any)
 buf = TRIM(gms_scr_path)//'/'//TRIM(datname)
 call delete_file(buf)
 if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

 write(longbuf,'(A,I0,A)') TRIM(inpname)//' 01 ',nproc,' >'//TRIM(gmsname)//" 2>&1"
! write(iout,'(A)') '$$GMS '//TRIM(longbuf)
 i = system(TRIM(gms_path)//' '//TRIM(longbuf))
 call read_gvb_energy_from_gms(gmsname, gvb_e)

 ! move the .dat file into current directory
 i = system('mv '//TRIM(gms_scr_path)//'/'//TRIM(datname)//' .')
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine do_gvb: fail to move file. Possibly&
                   & wrong gms_scr_path.'
  write(iout,'(A)') 'gms_scr_path='//TRIM(gms_scr_path)
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
 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine do_gvb: failed to call utility extract_noon2fch.'
  write(iout,'(A)') 'Did you delete it or forget to compile it?'
  stop
 end if

 ! update Total SCF Density in .fch(k) file
 call update_density_using_no_and_on(inpname)

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

! do uncontracted/ic-/FIC- MRCISD(+Q) for npair<=7, or <=CAS(14,14)
subroutine do_mrcisd()
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, casci, casscf, CIonly, ist, hf_fch, mrcisd,&
  mrcisd_prog, CtrType, casnofch, molcas_path, orca_path, gau_path, molpro_path,&
  psi4_path, bgchg, casci_prog, casscf_prog, chgname, F12, RI
 use mol, only: nbf, nif, npair, nopen, npair0, ndb, casci_e, casscf_e, davidson_e,&
  mrcisd_e, ptchg_e, nuc_pt_e
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
  i = system(TRIM(molcas_path)//' '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1")

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
  i = system(TRIM(orca_path)//' '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1")

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
  write(iout,'(A)') '$'//TRIM(string)
  i = system(TRIM(string))

 case('psi4')
  call check_exe_exist(psi4_path)

  i = system('fch2psi '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'_psi.inp'

  i = index(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.inp'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(string), TRIM(inpname))

  call prt_mrcisd_psi4_inp(inpname)
  write(string,'(A,I0)') 'psi4 '//TRIM(inpname)//' '//TRIM(outname)&
                         //' -n ',nproc
  write(iout,'(A)') '$'//TRIM(string)
  i = system(TRIM(string))
 
 case default
  write(iout,'(A)') 'ERROR in subroutine do_mrcisd: invalid program='//TRIM(mrcisd_prog)
  stop
 end select

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine do_mrcisd: error termination.'
  write(iout,'(A)') 'Filename='//TRIM(outname)
  stop
 end if

 ! read Davidson correction and MRCISD energy from OpenMolcas/ORCA/Gaussian output file
 call read_mrcisd_energy_from_output(CtrType, mrcisd_prog, outname, ptchg_e,&
      nuc_pt_e, davidson_e, e)
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
 case('gaussian','psi4') ! only uncontracted MRCISD
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
subroutine modify_memory_in_gms_inp(inpname, mem, nproc)
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: mem, nproc
 integer, parameter :: iout = 6
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
  write(iout,'(A)') "ERROR in subroutine modify_memory_in_gms_inp: no 'MWORDS'&
                   & found in file "//TRIM(inpname)
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
end subroutine modify_memory_in_gms_inp

