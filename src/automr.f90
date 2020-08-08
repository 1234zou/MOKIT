! written by jxzou at 20200420: black-box multireference calculations
! updated by jxzou at 20200511: framework of the program

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
  write(iout,'(A,/)') ' Example: automr a.gjf'
  stop
 end if

 call getarg(1, fname)

 i = index(fname, '.gjf', back=.true.)
 j = index(fname, '.fch', back=.true.)
 if(i/=0 .and. j/=0) then
  write(iout,'(A)') "ERROR in subroutine automr: both '.gjf' and '.fch' key detected&
                   & in filename "//TRIM(fname)
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
 use mr_keyword, only: gjfname, read_program_path, parse_keyword, check_kywd_compatible, prt_strategy
 implicit none
 integer :: i
 character(len=24) :: data_string
 character(len=240), intent(in) :: fname

 gjfname = fname
 call read_program_path() ! read in paths in $MOKIT_ROOT/program.info
 call parse_keyword()
 call check_kywd_compatible()

 call do_hf()
 call get_paired_LMO()
 call do_gvb()
 call do_cas(.false.)! CASCI/DMRG-CASCI
 call do_cas(.true.) ! CASSCF/DMRG-CASSCF
 call do_mrpt2()     ! CASPT2/NEVPT2
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
  mo_rhf, bgchg, read_bgchg_from_gjf, gjfname, chgname, prt_strategy, formchk
 implicit none
 integer :: i, system
 real(kind=8) :: ssquare = 0.0d0
 character(len=24) :: data_string = ' '
 character(len=240) :: rhf_gjfname, uhf_gjfname, chkname

 write(iout,'(//,A)') 'Enter subroutine do_hf...'

 if(skiphf) then
  write(iout,'(A)') 'Skip the RHF/UHF step...'
  call read_natom_from_fch(hf_fch, natom)
  allocate(coor(3,natom), elem(natom), nuc(natom))
  call read_elem_and_coor_from_fch(hf_fch, natom, elem, nuc, coor, charge, mult)
  if(bgchg) call read_bgchg_from_gjf(gjfname, .true.)
  call fdate(data_string)
  write(iout,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)
  return
 end if

 call read_natom_from_gjf(gjfname, natom)
 allocate(coor(3,natom), elem(natom), nuc(natom))
 call read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
 if(bgchg) call read_bgchg_from_gjf(gjfname, .false.)

 i = index(gjfname, '.gjf', back=.true.)
 rhf_gjfname = gjfname(1:i-1)//'_rhf.gjf'
 uhf_gjfname = gjfname(1:i-1)//'_uhf.gjf'

 if(mult == 1) then ! singlet, perform RHF and UHF
  call generate_hf_gjf(rhf_gjfname, natom, elem, coor, charge, mult, basis,&
                       .false., cart, mem, nproc)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(rhf_gjfname))
  call perform_hf_and_read_e(gau_path, rhf_gjfname, rhf_e, ssquare)

  call generate_hf_gjf(uhf_gjfname, natom, elem, coor, charge, mult, basis,&
                       .true., cart, mem, nproc)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_gjfname))
  call perform_hf_and_read_e(gau_path, uhf_gjfname, uhf_e, ssquare)

 else               ! not singlet, only perform UHF

  call generate_hf_gjf(uhf_gjfname, natom, elem, coor, charge, mult, basis,&
                       .true., cart, mem, nproc)
  call perform_hf_and_read_e(gau_path, uhf_gjfname, uhf_e, ssquare)
 end if

 if(mult == 1) then
  write(iout,'(/,A,F18.8,1X,A,F7.3)') 'E(RHF) = ', rhf_e, 'a.u., <S**2>=', 0.0
  write(iout,'(A,F18.8,1X,A,F7.3)') 'E(UHF) = ', uhf_e, 'a.u., <S**2>=',ssquare
 else
  write(iout,'(/,A,F18.8,1X,A,F7.3)') 'E(UHF) = ', uhf_e, 'a.u., <S**2>=',ssquare
 end if

 if(rhf_e - uhf_e > 1.0d-4) then
  write(iout,'(/,A)') 'UHF energy is lower, choose UHF wave function.'
  ist = 1
  mo_rhf = .false.
  i = index(gjfname, '.gjf', back=.true.)
  chkname = gjfname(1:i-1)//'_uhf.chk'
  hf_fch = gjfname(1:i-1)//'_uhf.fch'
  call formchk(chkname, hf_fch)
 else
  write(iout,'(/,A)') 'RHF/UHF is equal, or has little difference, choose RHF.'
  ist = 3
  mo_rhf = .true.
  i = index(gjfname, '.gjf', back=.true.)
  chkname = gjfname(1:i-1)//'_rhf.chk'
  hf_fch = gjfname(1:i-1)//'_rhf.fch'
  call formchk(chkname, hf_fch)
 end if

 write(iout,'(A)') 'Strategy updated:'
 call prt_strategy()

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)
 return
end subroutine do_hf

! generate PySCF input file .py from Gaussian .fch(k) file, and get paired LMOs
subroutine get_paired_LMO()
 use print_id, only: iout
 use mr_keyword, only: mo_rhf, ist, hf_fch, bgchg, chgname
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

  call prt_rhf_proj_script_into_py(pyname)
  call prt_auto_pair_script_into_py(pyname)
  i = system('python '//TRIM(pyname)//' >& '//TRIM(proname)//'_proj_loc_pair.out')
  call delete_file(chkname)

 else
  if(ist == 1) then
   write(iout,'(A)') 'Two sets of MOs, invoke UNO associated rotation.'
  else if(ist == 2) then
   write(iout,'(A)') 'Two sets of MOs, invoke UNO generation.'
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
  call prt_uno_script_into_py(pyname)

  if(ist == 1) call prt_assoc_rot_script_into_py(pyname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(pyname))
  i = system('python '//TRIM(pyname)//' >& '//TRIM(outname))

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
 use mr_keyword, only: mem, nproc, tencycle, hf_fch
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
 write(fid2,'(A)') "mol2.basis = 'sto-6g'"
 write(fid2,'(A)') 'mol2.build()'
 write(fid2,'(A)') 'mf2 = scf.RHF(mol2)'
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
                       hardwfn, crazywfn, hf_fch, casnofch, CIonly, casscf_force
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, pyname1, no_fch
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
   write(fid2,'(A,3(I0,A))') 'mc = mcscf.CASSCF(mf,', nacto, ',(',nacta,',',nactb,'))'
   write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*500, ' # MB'
  else ! DMRG-CASSCF
   write(fid2,'(A,3(I0,A))') 'mc = dmrgscf.DMRGSCF(mf,', nacto, ',(',nacta,',',nactb,'))'
   write(fid2,'(A,I0)') 'mc.fcisolver.maxM = ', maxM
   write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ', CEILING(DBLE(mem)/2.0d0), ' # GB'
  end if
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*500, ' # MB'
  !write(fid2,'(A,I0)') 'mc.fcisolver.num_thrds = ', nproc
  write(fid2,'(A)') 'mc.max_cycle = 200'
 else ! CASCI/DMRG-CASCI
  write(fid2,'(A,3(I0,A))') 'mc = mcscf.CASCI(mf,', nacto, ',(',nacta,',',nactb,'))'
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*500, ' # MB'
  if(casci) then
   write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*500, ' # MB'
  else
   !write(fid2,'(A,I0)') 'mc.fcisolver.num_thrds = ', nproc
   write(fid2,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=', maxM, ')'
   write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ', CEILING(DBLE(mem)/2.0d0), ' # GB'
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

 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.natorb = True'
 write(fid2,'(A)') 'mc.kernel()'

 write(fid2,'(/,A)') '# save NOs into .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(gvb_fch)//"', '"//TRIM(casnofch)//"')"
 write(fid2,'(A)') 'noon = np.zeros(nif)'
 write(fid2,'(A)') "py2fch('"//TRIM(casnofch)//"', nbf, nif, mc.mo_coeff, Sdiag, 'a', noon)"
 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
 return
end subroutine prt_cas_script_into_py

! print a CASCI or CASSCF .gjf file
subroutine prt_cas_gjf(gjfname, mem, nproc, nacto, nacte, scf, force)
 implicit none
 integer :: i, fid
 integer, intent(in) :: mem, nproc, nacto, nacte
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: scf, force

 i = index(gjfname, '.gjf', back=.true.)
 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(A,I0,A,I0,A)',advance='no') '#p CAS(',nacte,',',nacto,')'
 write(fid,'(A)',advance='no') ' chkbasis nosymm guess=read geom=allcheck int=nobasistransform'
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
 write(fid,'(A,/)') '#p chkbasis nosymm guess(read,only,save,NaturalOrbitals) geom=allcheck int=nobasistransform'
 close(fid)
 return
end subroutine prt_cas_gjf

subroutine prt_cas_gms_inp(inpname, ncore, scf)
 use mol, only: nacto, nacte, charge, mult
 use mr_keyword, only: mem, nproc, hardwfn, crazywfn
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

 write(fid2,'(A)') '  ICUT=10 $END'
 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',CEILING(DBLE(mem*125)/DBLE(nproc)), ' $END'

 if(scf) then   ! CASSCF
  write(fid2,'(3(A,I0),A)',advance='no') ' $DET'
 else
  write(fid2,'(3(A,I0),A)',advance='no') ' $CIDET'
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

! print CASCI/CASSCF keywords in to a given Molcas/OpenMolcas input file
subroutine prt_cas_molcas_inp(inpname, scf)
 use print_id, only: iout
 use mol, only: charge, mult, nacte, nacto
 implicit none
 integer :: i, fid
 logical, intent(in) :: scf
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == "&SCF") exit
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_cas_molcas_inp: no ''&
                   & found in file "//TRIM(inpname)
  stop
 end if

 BACKSPACE(fid) ! overwrite the &SCF section
 write(fid,'(A)') "&RASSCF"
 write(fid,'(A,I0)') 'Spin = ', mult
 write(fid,'(A,I0)') 'Charge = ', charge
 write(fid,'(A,I0)') 'nActEl= ', nacte
 write(fid,'(A,I0)') 'RAS2 = ', nacto
 if(.not. scf) write(fid,'(A)') 'CIonly'
 i = index(inpname, '.input', back=.true.)
 write(fid,'(A,/)') 'FILEORB = '//inpname(1:i-1)//'.INPORB'

 close(fid)
 return
end subroutine prt_cas_molcas_inp

! print CASCI/CASSCF keywords in to a given Molcas/OpenMolcas input file
subroutine prt_cas_orca_inp(inpname, scf)
 use print_id, only: iout
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc
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
 write(fid2,'(A)') '%casscf'
 write(fid2,'(A,I0)') ' nel ', nacte
 write(fid2,'(A,I0)') ' norb ', nacto
 if(scf) write(fid2,'(A)') ' maxiter 200'
 write(fid2,'(A,I0)') ' ActOrbs NatOrbs'
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
end subroutine prt_cas_orca_inp

! print NEVPT2 script into a given .py file
subroutine prt_nevpt2_script_into_py(pyname)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, casnofch, casci, casscf, maxM
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
 end do
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
 write(fid2,'(A,3(I0,A))') 'mc = mcscf.CASCI(mf,', nacto, ',(',nacta,',',nactb,'))'
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

! print CASTP2 keywords into MOLCAS/OpenMolcas .input file
subroutine prt_caspt2_script_into_input(inputname)
 use print_id, only: iout
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
  if(buf(2:4) == 'SCF') exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine prt_caspt2_script_into_input: fileIO error.'
  write(iout,'(A)') 'inputname = '//TRIM(inputname)
  stop
 end if

 write(fid2,'(A)') "&RASSCF"
 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid1,status='delete')

 write(fid2,'(A,I0,A)') 'nActEl= ', nacte, ' 0 0'
 write(fid2,'(A,I0)') 'RAS2 = ', nacto
 write(fid2,'(A)') 'CIonly'
 write(fid2,'(/,A)') "&CASPT2"
 write(fid2,'(A)') 'MultiState= 1 1'
 write(fid2,'(A,/)') 'Frozen= 0'
 close(fid2)
 i = RENAME(inputname1, inputname)
 return
end subroutine prt_caspt2_script_into_input

! perform GVB computation (only in Strategy 1,3)
subroutine do_gvb()
 use print_id, only: iout
 use mr_keyword, only: nproc, gms_path, gms_scr_path, mo_rhf, ist, hf_fch, gvb,&
                       datname, npair_wish, bgchg, chgname
 use mol, only: nbf, nif, ndb, nopen, npair, lin_dep, gvb_e, nacto, nacta, &
                nactb, nacte, npair0
 implicit none
 integer :: i, j, system, RENAME
 character(len=24) :: data_string = ' '
 character(len=240) :: buf, proname, inpname, gmsname, pair_fch
 character(len=300) :: longbuf = ' '

 if(.not. gvb) return
 write(iout,'(//,A)') 'Enter subroutine do_gvb...'

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
   write(iout,'(2(A,I0),A)') 'Attention: AutoMR recommends GVB(',npair,'), but&
    & user specifies GVB(',npair_wish,'). Trying to fulfill...'
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

 ! sort the GVB pairs by CI coefficients of the 1st NOs
 write(longbuf,'(A,5(1X,I0))') 'gvb_sort_pairs '//TRIM(datname), nbf, nif, ndb, nopen, npair
 i = system(TRIM(longbuf))

 ! generate corresponding .fch file from _s.dat file
 i = index(datname, '.dat')
 inpname = datname(1:i-1)//'_s.fch'
 datname = datname(1:i-1)//'_s.dat'
 call copy_file(pair_fch, inpname, .false.)
 write(longbuf,'(2(A,I0))') 'dat2fch '//TRIM(datname)//' '//TRIM(inpname)//' -gvb ',&
                             npair, ' -open ', nopen
 i = system(TRIM(longbuf))

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
  gms_scr_path, bgchg, chgname, casscf_force, formchk, unfchk, gbw2mkl, mkl2gbw
 use mol, only: nbf, nif, npair, nopen, npair0, ndb, casci_e, casscf_e, nacta, &
                nactb, nacto, nacte, gvb_e, mult, ptchg_e, nuc_pt_e, natom, grad
 implicit none
 integer :: i, j, idx1, idx2, nvir, system, RENAME
 real(kind=8) :: e(2)   ! e(1) is CASCI enery, e(2) is CASSCF energy
 character(len=10) :: cas_prog = ' '
 character(len=24) :: data_string = ' '
 character(len=240) :: buf, fchname, pyname, inpname, outname, proname, mklname
 logical, intent(in) :: scf

 if(scf) then
  if((.not. casscf) .and. (.not.dmrgscf)) return
 else
  if((.not. casci) .and. (.not.dmrgci)) return
 end if

 write(iout,'(//,A)') 'Enter subroutine do_cas...'
 if(ist == 5) then
  call read_no_info_from_fch(hf_fch, nbf, nif, ndb, nopen, nacta, nactb, nacto, nacte)
  ! if ist=5, read nbf, nif, nopen, nacto, ... variables from NO .fch(k) file
  i = nacte; j = nacto
 else
  i = 2*npair0 + nopen; j = i
 end if

 if(nacte_wish>0 .and. i/=nacte_wish) then
  write(iout,'(4(A,I0),A)') 'Attention: AutoMR recommends CAS(',i,'e,',j,'o), but&
   & user specifies CAS(',nacte_wish,'e,',nacto_wish, 'o). Trying to fulfill...'

  if(ist == 5) then
   !if() then
   !end if

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
 write(iout,'(A,2(I0,A))') '(',i,',',i,')'
 write(iout,'(A,I5,4X,A,I5)') 'doubly_occ=', idx1-1, 'nvir=', nvir
 write(iout,'(2(A,I0))') 'No. of active alpha/beta e = ', nacta,'/',nactb

 if(npair0 > 7) then
  if(scf) then
   casscf = .false.
   dmrgscf = .true.
   write(iout,'(A)') 'Attention: CASSCF is switched to DMRG-CASSCF due to active&
                    & space larger than (14,14).'
   if(casscf_prog /= 'pyscf') then
    write(iout,'(A)') 'ERROR in subroutine do_cas: DMRGSCF required. But&
                     & CASSCF_prog='//TRIM(casscf_prog)//'.'
    stop
   end if

  else ! not scf
   casci = .false.
   dmrgci = .true.
   write(iout,'(A)') 'Attention: CASCI is switched to DMRG-CASCI due to active&
                    & space larger than (14,14).'
   if(casci_prog /= 'pyscf') then
    write(iout,'(A)') 'ERROR in subroutine do_cas: DMRGCI required. But&
                     & CASCI_prog='//TRIM(casci_prog)//'.'
    stop
   end if
  end if
 end if

 if(ist<1 .or. ist>3) then
  write(iout,'(A)') 'ERROR in subroutine do_cas: not supported currently.'
  write(iout,'(A,I0)') 'ist=', ist
  stop
 end if

 if(ist==1 .or. ist==3) then
  i = index(datname, '.dat', back=.true.)
  fchname = datname(1:i-1)//'.fch'
  pyname = datname(1:i-1)//'.py'
 else ! ist=2, UHF -> UNO -> CASCI/CASSCF
  i = index(hf_fch, '.fch', back=.true.)
  fchname = hf_fch(1:i-1)//'_uno.fch'
  pyname = hf_fch(1:i-1)//'_uno.py'
  inpname = hf_fch(1:i-1)//'_uno.py2'
  i = RENAME(TRIM(pyname), TRIM(inpname))
  ! bas_fch2py will generate file '_uno.py', so we need to rename it to another filename
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

 if(cas_prog == 'pyscf') then
  i = system('bas_fch2py '//TRIM(fchname))
  inpname = TRIM(proname)//'.py'
  i = RENAME(TRIM(pyname), TRIM(inpname))
  call prt_cas_script_into_py(inpname, fchname, scf)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  if(casscf_force) i = system("echo 'mc.Gradients().kernel()' >> "//TRIM(inpname))
  j = index(inpname, '.py', back=.true.)
  outname = inpname(1:j-1)//'.out'
  i = system('python '//TRIM(inpname)//' >& '//TRIM(outname))

 else if(cas_prog == 'gaussian') then
  inpname = TRIM(proname)//'.gjf'
  outname = TRIM(proname)//'.log'
  call prt_cas_gjf(inpname, mem, nproc, nacto, nacte, scf, casscf_force)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call unfchk(fchname, proname//'.chk')
  i = system(TRIM(gau_path)//' '//TRIM(inpname))
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine do_cas: Gaussian CAS job failed!'
   write(iout,'(A)') 'Filename='//TRIM(inpname)
   stop
  end if
  call formchk(TRIM(proname)//'.chk', casnofch)

  if(mult /= 1) i = system('fch_u2r '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  inpname = casnofch(1:i-1)//'_r.fch'
  i = RENAME(inpname, casnofch)
  inpname = TRIM(proname)//'.gjf'

 else if(cas_prog == 'gamess') then
  inpname = TRIM(proname)//'.dat'
  buf = TRIM(gms_scr_path)//'/'//TRIM(inpname) ! delete the possible .dat file
  call delete_file(buf)
  ! do not use datname in the above three lines! because datname may be that of a GVB job

  i = system('fch2inp '//TRIM(fchname))
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
  if(scf) i = system('dat2fch '//TRIM(datname)//' '//TRIM(casnofch))

  ! transfer NOs from .dat to .fch
  write(buf,'(A,I0,1X,I0)') 'dat2fch '//TRIM(datname)//' '//TRIM(casnofch)//' -no ',1,idx2
  i = system(TRIM(buf))

 else if(cas_prog == 'openmolcas') then
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
  call prt_cas_molcas_inp(inpname, scf)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  if(casscf_force) i = system("echo '&ALASKA' >> "//TRIM(inpname))
  i = system(TRIM(molcas_path)//' '//TRIM(inpname)//" >& "//TRIM(outname))

 else if(cas_prog == 'orca') then
  i = system('fch2mkl '//TRIM(fchname))

  i = index(fchname, '.fch', back=.true.)
  outname = fchname(1:i-1)//'.mkl'
  mklname = TRIM(proname)//'.mkl'
  i = RENAME(TRIM(outname),TRIM(mklname))

  i = index(fchname, '.fch', back=.true.)
  outname = fchname(1:i-1)//'.inp'
  inpname = TRIM(proname)//'.inp'
  i = RENAME(TRIM(outname),TRIM(inpname))

  call prt_cas_orca_inp(inpname, scf)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(TRIM(proname)//'.mkl')
  call delete_file(mklname)
  if(casscf_force) i = system("sed -i '3,3s/TightSCF/TightSCF EnGrad/' "//TRIM(inpname))

  outname = TRIM(proname)//'.out'
  i = system(TRIM(orca_path)//' '//TRIM(inpname)//" >& "//TRIM(outname))

  ! make a copy of the .fch file to save NOs
  call copy_file(fchname, casnofch, .false.)
  call gbw2mkl(TRIM(proname)//'.gbw')
  i = system('mkl2fch '//TRIM(proname)//'.mkl '//TRIM(casnofch))
 end if

 ! extract NOONs from the output file and print them into .fch file
 ! If the CASCI_prog = Gaussian, the NOs are already in .fch(k) file
 if(.not. (cas_prog=='gaussian' .or. cas_prog=='openmolcas')) then
  write(buf,'(A,2(1X,I0))') 'extract_noon2fch '//TRIM(outname)//' '//&
                             TRIM(casnofch), idx1, idx2
  i = system(TRIM(buf))
  if(i /= 0) then
   write(iout,'(A)') 'Warning in subroutine do_cas: extract_noon2fch failed.'
   write(iout,'(A)') 'Filename = '//TRIM(outname)
   write(iout,'(A)') 'This does not affect the energy. So continue.'
  end if
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

 if(gvb .and. 2*npair+nopen==nacto .and. gvb_e<e(1)) then
  write(iout,'(A)') 'ERROR in subroutine do_cas: active space of GVB and CAS&
   & are equal, but CASCI/CASSCF energy is higher than that of GVB.'
  write(iout,'(A)') 'This is probably due to: (1) CASCI stucks in a higher energy&
   & local minimum or not pure spin state; (2) GVB MOs are disordered (GAMESS bug).'
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

! do CASCI/CASPT2(npair<=7) or DMRG-CASCI/DMRG-NEVPT2 (npair>7) when scf=.False.
! do CASSCF/CASPT2(npair<=7) or DMRG-CASSCF/DMRG-NEVPT2 (npair>7) when scf=.True.
subroutine do_mrpt2()
 use print_id, only: iout
 use mr_keyword, only: casci, casscf, dmrgci, dmrgscf, CIonly, caspt2,&
                       nevpt2, casnofch, casscf_prog, casci_prog, bgchg,&
                       chgname
 use mol, only: casci_e, casscf_e, caspt2_e, nevpt2_e
 implicit none
 integer :: i, RENAME, system
 character(len=24) :: data_string
 character(len=240) :: data_string1, pyname, outname, inputname
 character(len=240) :: fchname
 real(kind=8) :: e = 0.0d0

 if((.not.caspt2) .and. (.not.nevpt2)) return
 write(iout,'(//,A)') 'Enter subroutine do_mrpt2...'

 if(.not. CIonly) then
  if(casscf_prog == 'orca') then
   write(iout,'(A)') 'Attention: ORCA is used as the CASSCF solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(iout,'(A)') 'This will affect the PT2 energy at about 10^-5 a.u.&
                    & such small error is usually not important.'
   write(iout,'(A)') 'If you care about the accuracy, please use another CASSCF solver.'
  end if

  if(casscf) then
   if(caspt2) then
    data_string1 = 'CASPT2 based on optimized CASSCF orbitals.'
   else
    data_string1 = 'CASSCF-NEVPT2 based on optimized CASSCF orbitals.'
   end if
  else ! DMRG-CASSCF
   if(caspt2) then
    data_string1 = 'DMRG-CASPT2 based on optimized DMRG-CASSCF orbitals.'
   else
    data_string1 = 'DMRG-NEVPT2 based on optimized DMRG-CASSCF orbitals.'
   end if
  end if

 else ! CIonly = .True.
  if(casci_prog == 'orca') then
   write(iout,'(A)') 'Attention: ORCA is used as the CASCI solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(iout,'(A)') 'This will affect the PT2 energy at about 10^-5 a.u.&
                    & such small error is usually not important.'
   write(iout,'(A)') 'If you care about the accuracy, please use another CASCI solver.'
  end if

  if(casci) then
   if(caspt2) then
    data_string1 = 'CASPT2 based on CASCI orbitals.'
   else
    data_string1 = 'NEVPT2 based on CASCI orbitals.'
   end if
  else ! DMRG-CASCI
   if(caspt2) then
    data_string1 = 'DMRG-CASPT2 based on DMRG-CASCI orbitals.'
   else
    data_string1 = 'DMRG-NEVPT2 based on DMRG-CASCI orbitals.'
   end if
  end if
  write(iout,'(A)') 'Attention: the CASSCF orbital optimization is strongly&
                     & recommended'
  write(iout,'(A)') 'to be performed before PT2, unless it is too time-consuming.'
 end if
 write(iout,'(A)') TRIM(data_string1)

 write(iout,'(A)') 'Frozen_core = F'

 if(nevpt2) then
  i = system('bas_fch2py '//TRIM(casnofch))
  i = index(casnofch, '.fch')
  inputname = casnofch(1:i-1)//'.py'
  pyname = casnofch(1:i-1)//'_NEVPT2.py'
  outname = casnofch(1:i-1)//'_NEVPT2.out'
  i = RENAME(inputname, pyname)
  call prt_nevpt2_script_into_py(pyname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(pyname))

  i = system('python '//TRIM(pyname)//" >& "//TRIM(outname))
 else ! CASPT2
  i = system('fch2inporb '//TRIM(casnofch)//' -no') ! generate .input and .INPORB
  i = index(casnofch, '.fch', back=.true.)
  pyname = casnofch(1:i-1)//'.input'
  i = index(casnofch, '_NO', back=.true.)
  inputname = casnofch(1:i-1)//'_CASPT2.input'
  outname = casnofch(1:i-1)//'_CASPT2.out'
  i = RENAME(pyname, inputname)
  call prt_caspt2_script_into_input(inputname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inputname))

  i = system('pymolcas '//TRIM(inputname)//" >& "//TRIM(outname))
 end if

 if(i /= 0) then
  if(nevpt2) then
   write(iout,'(A)') 'ERROR in subroutine do_mrpt2: NEVPT2 computation failed.'
   write(iout,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
  else
   write(iout,'(A)') 'ERROR in subroutine do_mrpt2: CASPT2 computation failed.'
   write(iout,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
  end if
  stop
 end if

 if(nevpt2) then ! read NEVPT2 correlation energy
  call read_mrpt2_energy_from_pyscf_out(outname, e)
 else            ! read CASPT2 total energy
  call read_mrpt2_energy_from_molcas_out(outname, e)
 end if

 if(caspt2) then
  write(iout,'(A)') 'IP-EA shift = 0.25 (default)'
  if(.not. CIonly) then
   write(iout,'(/,A,F18.8,1X,A4)') 'E_corr(CASPT2) = ', e-casscf_e, 'a.u.'
  else
   write(iout,'(/,A,F18.8,1X,A4)') 'E_corr(CASPT2) = ', e-casci_e, 'a.u.'
  end if
  caspt2_e = e
  write(iout,'(A,F18.8,1X,A4)') 'E(CASPT2)      = ', caspt2_e, 'a.u.'
 else ! NEVPT2
  if(.not. CIonly) then
   nevpt2_e = e + casscf_e
  else
   nevpt2_e = e + casci_e
  end if
  write(iout,'(A,F18.8,1X,A4)') 'E(SC-NEVPT2)      = ', nevpt2_e, 'a.u.'
 end if

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_mrpt2 at '//TRIM(data_string)
 return
end subroutine do_mrpt2

! do ic-MRCC based on CASSCF, npair<=5
subroutine do_mrcc()
 use print_id, only: iout
 use mr_keyword, only: casci, casscf, dmrgci, dmrgscf, CIonly
 implicit none
 integer :: i


! write(iout,'(/,A)') 'Enter subroutine do_mrcc...'
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

subroutine perform_hf_and_read_e(gau_path, gjfname, e, ssquare)
 use print_id, only: iout
 implicit none
 integer :: i, fid, system
 real(kind=8), intent(out) :: e, ssquare
 character(len=240) :: buf, logname
 character(len=240), intent(in) :: gau_path
 character(len=240), intent(in) :: gjfname

 e = 0.0d0
 ssquare = 1.0d0
 i = index(gjfname, '.gjf', back=.true.)
 logname = gjfname(1:i-1)//'.log'

 i = system(TRIM(gau_path)//' '//TRIM(gjfname))
 if(i /= 0) then
  write(fid,'(/,A)') 'ERROR in subroutine perform_hf_and_read_e: running Gaussian failed.'
  write(fid,'(A)') 'You can open file '//TRIM(logname)//' and check why.'
  stop
 end if

 open(newunit=fid,file=TRIM(logname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(index(buf,'SCF Done') /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine perform_hf_and_read_e: no 'SCF Done'&
                   & found in file "//TRIM(logname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) e

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 i = index(buf, 'S**2')
 if(i /= 0) read(buf(i+6:),*) ssquare

 close(fid)
 return
end subroutine perform_hf_and_read_e

