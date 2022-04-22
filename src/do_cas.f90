! written by jxzou at 20210215: moved here from file automr.f90
! updated by jxzou at 20210216: implement first correct version of ist=5

module icss_param
 implicit none
 integer :: ngrid(3)
 real(kind=8) :: origin(3)
 real(kind=8), parameter :: r_ext = 6.5d0 ! 6.5 Angstrom
 real(kind=8), parameter :: interval = 0.3d0 ! 0.3 Angstrom
 ! this can be viewed as coarse grids
end module icss_param

! do CASCI(npair<=7) or DMRG-CASCI(npair>7) when scf=.False.
! do CASSCF(npair<=7) or DMRG-CASSCF(npair>7) when scf=.True.
subroutine do_cas(scf)
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, casci, dmrgci, casscf, dmrgscf, ist, hf_fch,&
  datname, nacte_wish, nacto_wish, gvb, casnofch, casci_prog, casscf_prog, &
  dmrgci_prog, dmrgscf_prog, gau_path, gms_path, molcas_path, orca_path, &
  gms_scr_path, molpro_path, bdf_path, bgchg, chgname, casscf_force,&
  check_gms_path, prt_strategy, RI, nmr, ICSS, dryrun, nstate
 use mol, only: nbf, nif, npair, nopen, npair0, ndb, casci_e, casscf_e, nacta, &
  nactb, nacto, nacte, gvb_e, ptchg_e, nuc_pt_e, natom, grad
 use util_wrapper, only: formchk, unfchk, gbw2mkl, mkl2gbw, fch2inp_wrap, &
  mkl2fch_wrap
 implicit none
 integer :: i, j, idx1, idx2, nvir, nfile, system, RENAME
 real(kind=8) :: unpaired_e ! unpaired electrons
 real(kind=8) :: e(2)   ! e(1) is CASCI enery, e(2) is CASSCF energy
 character(len=10) :: cas_prog = ' '
 character(len=24) :: data_string = ' '
 character(len=240) :: buf, fchname, pyname, inpname, outname, proname, mklname
 character(len=240) :: orbname, xmlname
 logical, intent(in) :: scf

 if(scf) then
  if((.not. casscf) .and. (.not.dmrgscf)) return
 else
  if((.not. casci) .and. (.not.dmrgci)) return
 end if
 write(iout,'(//,A)') 'Enter subroutine do_cas...'

 if(ist == 5) then
  write(iout,'(A)') 'Radical index for input NOs:'
  call calc_unpaired_from_fch(hf_fch, 3, .false., unpaired_e)
  ! read nbf, nif, nopen, nacto, ... variables from NO .fch(k) file
  call read_no_info_from_fch(hf_fch,nbf,nif,ndb,nopen,nacta,nactb,nacto,nacte)
  i = nacte; j = nacto
  npair0 = nactb; npair = npair0
 else
  i = 2*npair0 + nopen; j = i
 end if

 if(nacte_wish>0 .and. i/=nacte_wish) then
  write(iout,'(4(A,I0),A)') 'Warning: AutoMR recommends CAS(',i,'e,',j,'o), but&
   & user specifies CAS(',nacte_wish,'e,',nacto_wish, 'o). Trying to fulfill...'

  if(ist == 5) then
   if(nacte_wish > nacte) then
    write(iout,'(A)') 'ERROR in subroutine do_cas: too large active space required.'
    write(iout,'(2(A,I0),A)') 'Maximum allowed: CAS(',nacte,',',nacte,')'
    stop
   else if(nacte_wish < nacte) then
    i = nacte_wish; j = nacto_wish
    npair0 = (i-nopen)/2;   npair = npair0
    ndb = ndb + nactb - npair
    nacta = npair0 + nopen; nactb = npair0
    nacto = nacto_wish;     nacte = nacto_wish
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
  write(iout,'(/,A)') 'The reason is this molecule has little multi-configurational&
                     & or multi-reference'
  write(iout,'(A)') 'character. This molecule can be well described by single&
                   & reference methods, e.g.'
  write(iout,'(A)') 'MP2, CCSD(T). Thus no need for multi-reference computation&
                    &. But if you still want'
  write(iout,'(A)') 'to do it, you can manually specify the size of acitve spac&
                    &e in .gjf file. For ex-'
  write(iout,'(A)') 'ample, CASSCF(4,4) for water(H2O), or CASSCF(8,8) for meth&
                    &ane(CH4). The maximum'
  write(iout,'(A)') 'number of active orbitals is 2*npair. You can find npair &
                    &in GVB computations above.'
  write(iout,'(A,I0)') 'For this molecule, npair=', npair
  stop
 end if
 write(iout,'(2(A,I4,3X),A,L1,3X,A,I0)') 'doubly_occ=', idx1-1, 'nvir=',nvir,&
                                         'RIJK=', RI, 'nstate=', nstate
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
  write(iout,'(A,I0)') 'Allowed values are 1~5. But got ist=', ist
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
  fchname = hf_fch
  i = index(hf_fch, '.fch', back=.true.)
  pyname = hf_fch(1:i-1)//'.py'
 end if

 if(dryrun) return ! do not perform electronic structure calculations

 proname = ' '
 i = index(hf_fch, '.fch', back=.true.)
 if(ist==1 .or. ist==3) then
  if(scf) then
   write(proname,'(A,I0,A)') hf_fch(1:i-1)//'_gvb', npair, '_CASSCF'
  else
   write(proname,'(A,I0,A)') hf_fch(1:i-1)//'_gvb', npair, '_CASCI'
  end if
 else if(ist == 2) then
  if(scf) then
   write(proname,'(A)') hf_fch(1:i-1)//'_uno2CASSCF'
  else
   write(proname,'(A)') hf_fch(1:i-1)//'_uno2CASCI'
  end if
 else if(ist == 5) then
  if(scf) then
   proname = hf_fch(1:i-1)//'_CASSCF'
  else
   proname = hf_fch(1:i-1)//'_CASCI'
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

  write(buf,'(A)') 'python '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1"
  write(iout,'(A)') '$'//TRIM(buf)
  i = system(TRIM(buf))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_cas: PySCF CASCI/CASSCF job failed.'
   write(iout,'(A)') 'Please open file '//TRIM(outname)//' and check.'
   stop
  end if

 case('gaussian')
  call check_exe_exist(gau_path)

  inpname = TRIM(proname)//'.gjf'
  outname = TRIM(proname)//'.log'
  mklname = TRIM(proname)//'.chk'
  call prt_cas_gjf(inpname, nacto, nacte, scf, casscf_force)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call unfchk(fchname, mklname)

  write(iout,'(A)') '$'//TRIM(gau_path)//' '//TRIM(inpname)
  i = system(TRIM(gau_path)//' '//TRIM(inpname))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_cas: Gaussian CASCI/CASSCF job failed.'
   write(iout,'(A)') 'Please open file '//TRIM(inpname)//' and check.'
   stop
  end if
  call formchk(mklname, casnofch)
  call modify_IROHF_in_fch(casnofch, 0)

 case('gamess')
  call check_gms_path()
  inpname = TRIM(proname)//'.dat'
  buf = TRIM(gms_scr_path)//'/'//TRIM(inpname) ! delete the possible .dat file
  call delete_file(buf)
  ! do not use datname in the above three lines! because datname may be that of a GVB job

  call fch2inp_wrap(fchname, .false., 0, 0)
  i = index(fchname, '.fch', back=.true.)
  outname = fchname(1:i-1)//'.inp'
  inpname = TRIM(proname)//'.inp'
  i = RENAME(TRIM(outname), TRIM(inpname))
  outname = TRIM(proname)//'.gms'
  call prt_cas_gms_inp(inpname, idx1-1, scf)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  if(casscf_force) i = system("sed -i '1,1s/ENERGY/GRADIENT/' "//TRIM(inpname))

  write(buf,'(A,I0,A)') TRIM(inpname)//' 01 ',nproc,' >'//TRIM(outname)//" 2>&1"
  write(iout,'(A)') '$$GMS '//TRIM(buf)
  i = system(TRIM(gms_path)//' '//TRIM(buf))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_cas: GAMESS CASCI/CASSCF job failed.'
   write(iout,'(A)') 'Please open file '//TRIM(outname)//' and check.'
   stop
  end if

  ! make a copy of the .fch file to save NOs
  if(ist /= 2) then ! datname is a GVB job .dat file
   i = index(datname,'.dat')
   inpname = datname(1:i-1)//'.fch'
  else              ! no GVB job
   i = index(hf_fch, '.fch')
   inpname = hf_fch(1:i-1)//'_uno.fch'
  end if
  call copy_file(inpname, casnofch, .false.)

  ! move .dat file into current directory
  datname = TRIM(proname)//'.dat'
  i = system('mv '//TRIM(gms_scr_path)//'/'//TRIM(datname)//' .')

  ! transfer CASSCF pseudo-canonical MOs from .dat to .fch
  if(scf) then
   i = system('dat2fch '//TRIM(datname)//' '//TRIM(casnofch))
   write(buf,'(A,2(1X,I0))') 'extract_noon2fch '//TRIM(outname)//' '//&
                              TRIM(casnofch), idx1, idx2
   i = system(TRIM(buf))
  end if

  ! transfer NOs from .dat to .fch
  write(buf,'(A,I0)') 'dat2fch '//TRIM(datname)//' '//TRIM(casnofch)//' -no 1 ',idx2
  i = system(TRIM(buf))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_cas: failed to call utility dat2fch.'
   write(iout,'(A)') 'Related files: '//TRIM(datname)//', '//TRIM(casnofch)//'.'
   stop
  end if

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

  write(buf,'(A)') 'pymolcas '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1"
  write(iout,'(A)') '$'//TRIM(buf)
  i = system(TRIM(buf))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_cas: OpenMolcas CASCI/CASSCF job failed.'
   write(iout,'(A)') 'Please open file '//TRIM(outname)//' and check.'
   stop
  end if

  ! make a copy of the .fch file to save NOs
  call copy_file(fchname, casnofch, .false.)

  ! transfer NOs from .dat to .fch
  i = system('orb2fch '//TRIM(orbname)//' '//TRIM(casnofch)//' -no')

 case('orca')
  call check_exe_exist(orca_path)
  i = system('fch2mkl '//TRIM(fchname))

  i = index(fchname, '.fch', back=.true.)
  pyname  = fchname(1:i-1)//'_o.inp'
  orbname = fchname(1:i-1)//'_o.mkl'
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

  write(buf,'(A)') TRIM(inpname)//' >'//TRIM(outname)//" 2>&1"
  write(iout,'(A)') '$$ORCA '//TRIM(buf)
  i = system(TRIM(orca_path)//' '//TRIM(buf))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_cas: ORCA CASCI/CASSCF job failed.'
   write(iout,'(A)') 'Please open file '//TRIM(outname)//' and check.'
   stop
  end if

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
  call mkl2fch_wrap(mklname, casnofch, .true.)

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
  write(buf,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -t 1 -m ', i,&
                           'm '//TRIM(inpname)
  write(iout,'(2(A,I0),A)') '$molpro -n ',nproc,' -t 1 -m ',i,'m '//TRIM(inpname)
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

  write(iout,'(A)') '$$BDF '//TRIM(proname)
  i = system(TRIM(bdf_path)//' '//TRIM(proname))
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine do_cas: BDF CASCI/CASSCF job failed.'
   stop
  end if
  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  i = system('bdf2fch '//TRIM(orbname)//' '//TRIM(casnofch)//' -no')

 case('psi4')
  i = system('fch2psi '//TRIM(fchname))
  i = index(fchname, '.fch', back=.true.)
  mklname = fchname(1:i-1)//'_psi.inp'
  pyname  = fchname(1:i-1)//'_psi.A'
  inpname = TRIM(proname)//'.inp'
  xmlname = TRIM(proname)//'.A'
  outname = TRIM(proname)//'.out'
  i = RENAME(TRIM(mklname), TRIM(inpname))
  i = RENAME(TRIM(pyname), TRIM(xmlname))
  call prt_cas_psi4_inp(inpname, scf, casscf_force)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  write(buf,'(A,I0)') 'psi4 '//TRIM(inpname)//' '//TRIM(outname)//' -n ', nproc
  write(iout,'(A)') '$'//TRIM(buf)
  i = system(TRIM(buf))
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine do_cas: PSI4 CASCI/CASSCF job failed.'
   write(iout,'(A)') 'Please open file '//TRIM(outname)//' and check.'
   stop
  end if

  write(buf,'(A,2(1X,I0))') 'extract_noon2fch '//TRIM(outname)//' '//&
                             TRIM(casnofch), idx1, idx2
  i = system(TRIM(buf))

 case('dalton')
  write(6,'(/,A)') 'Warning: you should use OpenMP version of Dalton. If you &
                   &are using MPI-'
  write(6,'(A,/)') 'parallelized Dalton, you should kill the job.'
  i = system('fch2dal '//TRIM(fchname))
  i = index(fchname, '.fch', back=.true.)
  mklname = fchname(1:i-1)//'.dal'
  pyname  = fchname(1:i-1)//'.mol'
  inpname = TRIM(proname)//'.dal'
  xmlname = TRIM(proname)//'.mol'
  orbname = TRIM(proname)//'.MOPUN'
  outname = TRIM(proname)//'.out'
  i = RENAME(TRIM(mklname), TRIM(inpname))
  i = RENAME(TRIM(pyname), TRIM(xmlname))
  call prt_cas_dalton_inp(inpname, scf, casscf_force)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  call submit_dalton_job(proname, mem, nproc, .false., .false., .false.)

  ! untar/unzip the compressed package
  i = system('tar -xpf '//TRIM(proname)//'.tar.gz DALTON.MOPUN SIRIUS.RST')
  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  i = system('dal2fch DALTON.MOPUN '//TRIM(casnofch)//' -no')
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_cas: failed to call utility dal2fch.'
   write(iout,'(A)') 'Please open files DALTON.MOPUN and '//TRIM(casnofch)//&
                    &' and check.'
   stop
  end if
  orbname = 'DALTON.MOPUN'
  call delete_file(orbname)

 case default
  write(iout,'(A)') 'ERROR in subroutine do_cas: Allowed programs are Gaussian&
                   & Gaussian, GAMESS, PySCF,'
  write(iout,'(A)') 'OpenMolcas, ORCA, Molpro, BDF, PSI4 and Dalton. But got &
                    &CAS_prog='//TRIM(cas_prog)
  stop
 end select

 ! generate CASCI/CASSCF density from NOs and NOONs
 ! Note: density in .fch of Gaussian CASSCF job is wrong (Gaussian bug), I have
 !  to re-generate density
 if(TRIM(cas_prog) /= 'pyscf') call update_density_using_no_and_on(casnofch)

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
  write(iout,'(/,A)') 'ERROR in subroutine do_cas: active space of GVB and CAS&
                     & are equal, but CASCI/CASSCF'
  write(iout,'(A)') 'energy is higher than that of GVB. This is probably due to:&
                   & (1) CASCI stucks in a higher'
  write(iout,'(A)') 'energy local minimum or not pure spin state; (2) GVB MOs are&
                   & disordered (GAMESS bug).'
  stop
 end if

 casci_e = e(1)
 write(iout,'(/,A,F18.8,1X,A4)') 'E(CASCI)  = ', e(1), 'a.u.'
 if(scf) then
  casscf_e = e(2)
  write(iout,'(A,F18.8,1X,A4)') 'E(CASSCF) = ', e(2), 'a.u.'
 end if

 call calc_unpaired_from_fch(casnofch, 3, .false., unpaired_e)

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

 if(nmr) then
  call prt_cas_dalton_nmr_inp(casnofch, scf, ICSS, nfile)
  inpname = TRIM(proname)//'_NMR'
  call submit_dalton_job(inpname, mem, nproc, .false., .true., .false.)
  inpname = TRIM(proname)//'_NMR.out'
  call read_shieldings_from_dalton_out(inpname)
 end if
 if(ICSS) then
  call submit_dalton_icss_job(proname, mem, nproc, nfile)
  call gen_icss_cub(proname, nfile)
 end if

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_cas at '//TRIM(data_string)
 return
end subroutine do_cas

! print CASCI/DMRG-CASCI or CASSCF/DMRG-CASSCF script into a given .py file
subroutine prt_cas_script_into_py(pyname, gvb_fch, scf)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, casci, dmrgci, casscf, dmrgscf, maxM,&
  hardwfn, crazywfn, casnofch, dkh2_or_x2c, RI, RIJK_bas
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, pyname1, cmofch
 character(len=240), intent(in) :: pyname, gvb_fch
 logical, intent(in) :: scf

 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
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
 if(scf) then
  if(casscf) then ! CASSCF
   write(fid2,'(3(A,I0),A)',advance='no') 'mc = mcscf.CASSCF(mf,', nacto,',(',nacta,',',nactb,')'
   if(dkh2_or_x2c) then
    write(fid2,'(A)') ').x2c1e()'
   else
    if(RI) then
     write(fid2,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
    else
     write(fid2,'(A)') ')'
    end if
   end if
   write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ',mem*200,' # MB'
  else ! DMRG-CASSCF
   write(fid2,'(A)',advance='no') 'mc = dmrgscf.DMRGSCF(mf,'
   if(dkh2_or_x2c) then
    write(fid2,'(3(I0,A))') nacto,',(',nacta,',',nactb,')).x2c1e()'
   else
    write(fid2,'(3(I0,A))') nacto,',(',nacta,',',nactb,'))'
   end if
   write(fid2,'(A,I0)') 'mc.fcisolver.maxM = ', maxM
   write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ',CEILING(DBLE(mem)/DBLE((5*nproc))),' # GB'
  end if
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*800, ' # MB'
  write(fid2,'(A)') 'mc.max_cycle = 200'
 else ! CASCI/DMRG-CASCI
  write(fid2,'(3(A,I0),A)',advance='no') 'mc = mcscf.CASCI(mf,', nacto,',(',nacta,',',nactb,')'
  if(dkh2_or_x2c) then
   write(fid2,'(A)') ').x2c1e()'
  else
   if(RI) then
    write(fid2,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
   else
    write(fid2,'(A)') ')'
   end if
  end if
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*800, ' # MB'
  if(casci) then ! CASCI
   write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*200, ' # MB'
  else           ! DMRG-CASCI
   write(fid2,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=', maxM, ')'
   write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ',CEILING(DBLE(mem)/DBLE((5*nproc))),' # GB'
  end if
 end if

 write(fid2,'(A,I0)') 'mc.fcisolver.spin = ', nacta-nactb

 if(hardwfn) then
  write(fid2,'(A,I0,A)') 'mc.fix_spin_(ss=',nacta-nactb,')'
  write(fid2,'(A)') 'mc.fcisolver.max_cycle = 200'
 else if(crazywfn) then
  write(fid2,'(A,I0,A)') 'mc.fix_spin_(ss=',nacta-nactb,')'
  write(fid2,'(A)') 'mc.fcisolver.level_shift = 0.2'
  write(fid2,'(A)') 'mc.fcisolver.pspace_size = 1200'
  write(fid2,'(A)') 'mc.fcisolver.max_space = 100'
  write(fid2,'(A)') 'mc.fcisolver.max_cycle = 300'
 else
  write(fid2,'(A)') 'mc.fcisolver.max_cycle = 100'
 end if

 ! For DMRG-CASCI/CASSCF, both the original MOs and NOs will be saved/punched.
 ! Since DMRG is not invariant to unitary rotations of orbitals, I hope the
 ! original MOs to be used in DMRG-NEVPT2/CASPT2 computations. NOs are often
 ! more delocalized than original MOs.
 if(.not. dmrgscf) write(fid2,'(A)') 'mc.natorb = True'
 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.kernel()'

 i = index(casnofch, '_NO', back=.true.)
 cmofch = casnofch(1:i)//'CMO.fch'
 if(dmrgci) then
  write(fid2,'(/,A)') '# copy original MOs into a new file'
  write(fid2,'(A)') "copyfile('"//TRIM(gvb_fch)//"', '"//TRIM(cmofch)//"')"
 else if(dmrgscf) then
  write(fid2,'(/,A)') '# save CMOs into .fch file'
  write(fid2,'(A)') "copyfile('"//TRIM(gvb_fch)//"', '"//TRIM(cmofch)//"')"
  write(fid2,'(A)') 'noon = np.zeros(nif)'
  write(fid2,'(A)') "py2fch('"//TRIM(cmofch)//"',nbf,nif,mc.mo_coeff,'a',noon,False)"
  write(fid2,'(/,A)') 'mc.natorb = True'
  write(fid2,'(A)') 'mc.kernel()'
 end if

 write(fid2,'(/,A)') '# save NOs into .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(gvb_fch)//"', '"//TRIM(casnofch)//"')"
 write(fid2,'(A)') "py2fch('"//TRIM(casnofch)//"',nbf,nif,mc.mo_coeff,'a',mc.mo_occ,True)"
 ! mc.mo_occ only exists for PySCF >= 1.7.4
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
 write(fid,'(2(A,I0),A)',advance='no') '#p CAS(',nacte,',',nacto,')'
 write(fid,'(A)',advance='no') ' chkbasis nosymm guess=read geom=allcheck'

 if(dkh2_or_x2c) then
  write(fid,'(A)',advance='no') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A)',advance='no') ' int=nobasistransform'
 end if
 if(force) write(fid,'(A)',advance='no') ' force'

 if(scf) then ! CASSCF
  write(fid,'(A,/)') ' scf(maxcycle=128)'
 else         ! CASCI
  write(fid,'(A,/)') ' scf(maxcycle=-2)'
  ! to obtain CASCI NOs, we need to use -2, since -1 only calculates CASCI energy
 end if

 write(fid,'(A)') '--Link1--'
 write(fid,'(A)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(A)',advance='no') '#p chkbasis nosymm guess(read,only,save,&
                               &NaturalOrbitals) geom=allcheck'
 if(dkh2_or_x2c) then
  write(fid,'(A,/)') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A,/)') ' int=nobasistransform'
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
 use mr_keyword, only: maxM, dmrgci, dmrgscf, RI, RIJK_bas, mokit_root,&
  hardwfn, crazywfn
 implicit none
 integer :: i, j, fid1, fid2, RENAME, system
 logical, intent(in) :: scf
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)

 inpname1 = TRIM(inpname)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:7) == 'SEWARD') exit
  j = index(buf, '...')
  if(j > 0) then
   if(RI) then
    j = index(buf, '.')
    write(fid2,'(A)') buf(1:j)//TRIM(RIJK_bas1)//'..'//TRIM(buf(j+3:))
   else
    write(fid2,'(A)') TRIM(buf)
   end if
  else
   write(fid2,'(A)') TRIM(buf)
  end if
 end do

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_cas_molcas_inp: no 'SEWARD'&
                   & found in file "//TRIM(inpname)
  stop
 end if

 if(RI) write(fid2,'(A)') 'RIJK'
 write(fid2,'(A)') "&SEWARD"

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == "&SCF") exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

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
 write(fid2,'(A)') 'CIMX= 200'
 write(fid2,'(A)') 'Tight= 5d-8 5d-6'
 if(crazywfn) then
  write(fid2,'(A)') 'SDav= 500'
 else if(hardwfn) then
  write(fid2,'(A)') 'SDav= 300'
 end if
 if(.not. scf) write(fid2,'(A)') 'CIonly'
 i = index(inpname, '.input', back=.true.)
 write(fid2,'(A,/)') 'FILEORB = '//inpname(1:i-1)//'.INPORB'

 if(dmrgci .or. dmrgscf) then
  write(fid2,'(A)') 'DMRG'
  write(fid2,'(A)') 'RGinput'
  write(fid2,'(A)') ' conv_thresh = 1E-8'
  write(fid2,'(A)') ' nsweeps = 20'
  write(fid2,'(A,I0)') ' max_bond_dimension = ', MaxM
  write(fid2,'(A)') 'endRG'
 end if

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))

 ! if RIJK is on, we need to generate the fitting basis set file for OpenMolcas
 if(RI) then
  i = system('cp '//TRIM(mokit_root)//'/basis/'//TRIM(RIJK_bas1)//' .')
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine prt_cas_molcas_inp: failed to copy&
                    & file from'
   write(iout,'(A)') TRIM(mokit_root)//'/basis/'//TRIM(RIJK_bas1)//' to '//&
                     ' current directory.'
   stop
  end if
  i = system('bas_gau2molcas '//TRIM(RIJK_bas1))
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine prt_cas_molcas_inp: failed to call&
                    & utility bas_gau2molcas.'
   write(iout,'(A)') 'Did you forget to compile it?'
   stop
  end if

  call delete_file(RIJK_bas1)
  call upper(RIJK_bas1)
  i = system('mv '//TRIM(RIJK_bas1)//' $MOLCAS/basis_library/jk_Basis/')
 end if

 return
end subroutine prt_cas_molcas_inp

! print CASCI/CASSCF keywords in to a given ORCA .inp file
subroutine prt_cas_orca_inp(inpname, scf)
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc, dkh2_or_x2c, RI, RIJK_bas
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
 write(fid2,'(A,I0)') '%maxcore ', CEILING(1000d0*DBLE(mem)/DBLE(nproc))

 read(fid1,'(A)') buf   ! skip '!' line
 write(fid2,'(A)',advance='no') '!'
 if(.not. scf) write(fid2,'(A)',advance='no') ' noiter'
 ! RIJK in CASSCF must be combined with CONVentional
 if(RI) write(fid2,'(A)',advance='no') ' RIJK conv '//TRIM(RIJK_bas)
 write(fid2,'(A)') ' TightSCF'

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
 end if
 ! Q: Why not use %casscf plus NoIter for CASCI?
 ! A: The CASCI NOONs cannot be obtained in that way, so I have to use %mrci

 if(RI) write(fid2,'(A)') ' TrafoStep RI'
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

! print CASCI/CASSCF keywords into a given Molpro input file
subroutine prt_cas_molpro_inp(inpname, scf, force)
 use print_id, only: iout
 use mol, only: ndb, npair, npair0, nacto
 use mr_keyword, only: RI, RIJK_bas
 implicit none
 integer :: fid, nclosed, nocc
 character(len=21) :: RIJK_bas1
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

 write(fid,'(A)',advance='no') '{'
 if(RI) write(fid,'(A)',advance='no') 'DF-'
 write(fid,'(A)',advance='no') 'CASSCF'
 if(RI) then
  call auxbas_convert(RIJK_bas, RIJK_bas1, 2)
  write(fid,'(A)',advance='no') ',df_basis='//TRIM(RIJK_bas1)//',df_basis_exch='&
                                //TRIM(RIJK_bas1)
 end if

 ! Note: we need 'NoExtra' to completely close symmetry.
 ! Otherwise the CASCI energy is slightly different to that of other programs
 if(scf) then
  write(fid,'(2(A,I0),A)') ';closed,',nclosed,';occ,',nocc,';NoExtra}'
 else
  write(fid,'(2(A,I0),A)') ';closed,',nclosed,';occ,',nocc,&
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
 use mol, only: nbf, nif, charge, mult, ndb, npair, npair0, nacto, nacte
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
 if(nbf > nif) then
  write(fid,'(A)') 'CheckLin'
  write(fid,'(A)') 'TolLin'
  write(fid,'(A)') '1.D-6'
 end if
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

! print CASCI/CASSCF keywords into a given PSI4 input file
subroutine prt_cas_psi4_inp(inpname, scf, force)
 use mol, only: ndb, npair, npair0, nacte
 use mr_keyword, only: mem, hardwfn, crazywfn, RI, RIJK_bas
 implicit none
 integer :: i, nclosed, fid
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, casnofch
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf, force

 call modify_memory_in_psi4_inp(inpname, mem)
 i = index(inpname, '.inp', back=.true.)
 casnofch = inpname(1:i-1)//'_NO.fch'
 nclosed = ndb + npair - npair0

 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(1:1) == '}') exit
 end do ! for while

 BACKSPACE(fid)
 if(RI) then
  call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
  write(fid,'(A)') ' scf_type df'
  write(fid,'(A)') ' df_basis_scf '//TRIM(RIJK_bas1)
  write(fid,'(A)') ' mcscf_type df'
  write(fid,'(A)') ' df_basis_mcscf '//TRIM(RIJK_bas1)
 end if

 if(crazywfn) then
  write(fid,'(A,I0)') ' ci_maxiter 500'
  write(fid,'(A)') ' H0_guess_size 1500'
  write(fid,'(A)') ' H0_blocksize 1500'
 else if(hardwfn) then
  write(fid,'(A,I0)') ' ci_maxiter 300'
  write(fid,'(A)') ' H0_guess_size 1200'
  write(fid,'(A)') ' H0_blocksize 1200'
 else
  write(fid,'(A,I0)') ' ci_maxiter 200'
 end if

 write(fid,'(A,I0,A)') ' restricted_docc [', nclosed, ']'
 write(fid,'(A,I0,A)') ' active [', nacte, ']'
 write(fid,'(A)') ' canonicalize_inactive_favg true'
 write(fid,'(A)') ' nat_orbs true'
 write(fid,'(A)') '}'

 if(scf) then
  write(fid,'(/,A)') "casscf_energy, cas_wfn = energy('casscf',ref_wfn=scf_wfn,&
                     &return_wfn=True)"
 else
  write(fid,'(/,A)') "casci_energy, cas_wfn = energy('fci',ref_wfn=scf_wfn,&
                     &return_wfn=True)"
 end if

 write(fid,'(A)') "fchk(cas_wfn,'"//TRIM(casnofch)//"')"
 return
end subroutine prt_cas_psi4_inp

! print CASCI/CASSCF keywords into a given Dalton input file
subroutine prt_cas_dalton_inp(inpname, scf, force)
 use print_id, only: iout
 use mol, only: mult, ndb, npair, npair0, nacto, nacte
 implicit none
 integer :: i, nclosed, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf, force

 inpname1 = TRIM(inpname)//'.t'
 nclosed = ndb + npair - npair0

 open(newunit=fid1,file=TRIM(inpname1),status='replace')
 write(fid1,'(A)') '**DALTON INPUT'
 write(fid1,'(A)') '.RUN WAVE FUNCTION'
 write(fid1,'(A,/,A)') '**WAVE FUNCTIONS','.HF'
 if(scf) then
  write(fid1,'(A)') '.MCSCF'
 else
  write(fid1,'(A)') '.CI'
 end if
 write(fid1,'(A,/,A)') '*SCF INPUT','.NODIIS'
 write(fid1,'(A,/,A)') '.NOQCSCF','.NONCANONICAL'
 write(fid1,'(A,/,A)') '*CONFIGURATION INPUT', '.SPIN MULTIPLICITY'
 write(fid1,'(I0)') mult
 write(fid1,'(A,/,I0)') '.INACTIVE ORBITALS', nclosed
 write(fid1,'(A,/,I0)') '.CAS SPACE', nacto
 write(fid1,'(A,/,I0)') '.ELECTRONS', nacte
 if(scf) then
  write(fid1,'(A)') '*OPTIMIZATION'
  write(fid1,'(A,/,A)') '.MAX CI','500'
  write(fid1,'(A,/,A)') '.MAX MACRO ITERATIONS','50'
  write(fid1,'(A,/,A)') '.MAX MICRO ITERATIONS','200'
 else
  write(fid1,'(A)') '*CI INPUT'
  write(fid1,'(A,/,A)') '.MAX ITERATIONS','500'
 end if
 write(fid1,'(A,/,A)') '*PRINT LEVELS', '.CANONI'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:12) == '*ORBITAL INP') exit
 end do ! for while

 BACKSPACE(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12) == '**END OF INP') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_cas_dalton_inp: no '**END OF INP'&
                   & found in"
  write(iout,'(A)') 'file '//TRIM(inpname)//'.'
  close(fid1,status='delete')
  close(fid)
  stop
 end if

 write(fid1,'(A)') '**END OF INPUT'
 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine prt_cas_dalton_inp

! print CASSCF NMR keywords into a given Dalton input file
subroutine prt_cas_dalton_nmr_inp(fchname, scf, ICSS, nfile)
 use mol, only: mult, ndb, npair, npair0, nacto, nacte
 use util_wrapper, only: fch2dal_wrap
 implicit none
 integer :: i, nclosed, fid, fid1, RENAME
 integer, intent(out) :: nfile
 character(len=240) :: buf, dalname, molname, olddal, oldmol
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: scf, ICSS

 i = index(fchname, '.fch', back=.true.)
 olddal = fchname(1:i-1)//'.dal'
 oldmol = fchname(1:i-1)//'.mol'

 i = index(fchname, '_NO.fch', back=.true.)
 if(i == 0) i = index(fchname, '.fch')
 dalname = fchname(1:i-1)//'_NMR.dal'
 molname = fchname(1:i-1)//'_NMR.mol'
 call fch2dal_wrap(fchname)
 i = RENAME(TRIM(oldmol), TRIM(molname))

 nclosed = ndb + npair - npair0
 open(newunit=fid,file=TRIM(dalname),status='replace')
 write(fid,'(A)') '**DALTON INPUT'
 write(fid,'(A)') '.RUN PROPERTIES'
 write(fid,'(A)') '**WAVE FUNCTIONS'
 if(scf) then
  write(fid,'(A)') '.MCSCF'
 else
  write(fid,'(A)') '.CI'
 end if
 write(fid,'(A,/,A)') '*CONFIGURATION INPUT', '.SPIN MULTIPLICITY'
 write(fid,'(I0)') mult
 write(fid,'(A,/,I0)') '.INACTIVE ORBITALS', nclosed
 write(fid,'(A,/,I0)') '.CAS SPACE', nacto
 write(fid,'(A,/,I0)') '.ELECTRONS', nacte
 if(scf) then
  write(fid,'(A)') '*OPTIMIZATION'
  write(fid,'(A,/,A)') '.MAX CI','500'
  write(fid,'(A,/,A)') '.MAX MICRO ITERATIONS','200'
 end if
 write(fid,'(A)') '*ORBITAL INPUT'
 write(fid,'(A)') '.MOSTART'

 open(newunit=fid1,file=TRIM(olddal),status='old',position='rewind')
 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:8) == '.MOSTART') exit
 end do ! for while

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:5) == '**END') exit
  write(fid,'(A)') TRIM(buf)
 end do ! for while
 close(fid1,status='delete')

 write(fid,'(A,/,A)') '**PROPERTIES','.SHIELD'
 write(fid,'(A)') '**END OF INPUT'
 close(fid)

 nfile = 0
 if(ICSS) call add_ghost2dalton_mol(dalname, nfile)
end subroutine prt_cas_dalton_nmr_inp

! add ghost atoms into Dalton .mol file
subroutine add_ghost2dalton_mol(inpname, nfile)
 use mol, only: natom, coor
 use icss_param, only: r_ext, interval, ngrid, origin
 implicit none
 integer :: i, j, k, m, n, p, fid, tot_ngrid
 integer, intent(out) :: nfile
 real(kind=8) :: x,y,z, bound(2,3)
 character(len=240) :: buf, proname, molname, molname1, dalname, dalname1
 character(len=240), intent(in) :: inpname

 do i = 1, 3   ! x,y,z
  bound(1,i) = MINVAL(coor(i,:)) - r_ext
  bound(2,i) = MAXVAL(coor(i,:)) + r_ext
  ngrid(i) = INT((bound(2,i) - bound(1,i))/interval)
 end do ! for i
 origin = bound(1,:)
 tot_ngrid = PRODUCT(ngrid)
 write(6,'(A,I0)') 'ICSS rectangular box range: bound atoms + 6.5 Angstrom. To&
                   &tal grid = ', tot_ngrid
 write(6,'(A,F6.3,A,3I5)') 'Interval=', interval, '. Ngrid along x,y,z are', &
                           (ngrid(i),i=1,3)

 i = index(inpname, '.dal', back=.true.)
 j = index(inpname, '_NMR.dal', back=.true.)
 proname = inpname(1:j-1)
 molname = inpname(1:i-1)//'.mol'
 dalname = inpname(1:i-1)//'.dal'
 molname1 = inpname(1:j-1)//'00001.mol'
 dalname1 = inpname(1:j-1)//'00001.dal'
 call copy_file(dalname, dalname1, .false.)
 call copy_mol_and_add_atomtypes(molname, molname1)
 open(newunit=fid,file=TRIM(molname1),status='old',position='append')

 n = 333 - natom
 write(fid,'(A,I0,A)') 'Charge=0. Atoms=', n, ' Basis=INTGRL Ghost'
 m = 0; p = 1

 do i = 1, ngrid(1), 1
  x = bound(1,1) + DBLE(i-1)*interval
  do j = 1, ngrid(2), 1
   y = bound(1,2) + DBLE(j-1)*interval
   do k = 1, ngrid(3), 1
    z = bound(1,3) + DBLE(k-1)*interval
    m = m + 1
    write(fid,'(A,3(1X,F15.8))') 'Bq',x,y,z
    if(m == n) then
     close(fid)
     tot_ngrid = tot_ngrid - n
     if(tot_ngrid == 0) exit
     m = 0; p = p + 1
     write(molname1,'(A,I5.5,A)') TRIM(proname),p,'.mol'
     write(dalname1,'(A,I5.5,A)') TRIM(proname),p,'.dal'
     call copy_file(dalname, dalname1, .false.)
     call copy_mol_and_add_atomtypes(molname, molname1)
     open(newunit=fid,file=TRIM(molname1),status='old',position='append')
     if(tot_ngrid < n) then
      write(fid,'(A,I0,A)') 'Charge=0. Atoms=',tot_ngrid, ' Basis=INTGRL Ghost'
     else
      write(fid,'(A,I0,A)') 'Charge=0. Atoms=', n, ' Basis=INTGRL Ghost'
     end if
    end if
   end do ! for k
  end do ! for j
 end do ! for i

 close(fid)
 nfile = p
end subroutine add_ghost2dalton_mol

! copy Dalton .mol file and add AtomTypes by 1
subroutine copy_mol_and_add_atomtypes(molname, molname1)
 implicit none
 integer :: i, natom, fid, fid1
 character(len=240) :: buf
 character(len=240), intent(in) :: molname, molname1

 open(newunit=fid,file=TRIM(molname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(molname1),status='replace')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:10) == 'AtomTypes=') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 i = index(buf, ' ')
 read(buf(11:i-1),*) natom
 write(fid1,'(A,I0,A)') 'AtomTypes=', natom+1, ' Integrals=1.0D-14 Charge=0 &
                        &NoSymmetry Angstrom'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid)
 close(fid1)
end subroutine copy_mol_and_add_atomtypes

subroutine submit_dalton_job(proname, mem, nproc, sirius, noarch, del_sout)
 implicit none
 integer :: i, system
 integer, intent(in) :: mem, nproc
 character(len=480) :: buf
 character(len=240), intent(in) :: proname
 logical, intent(in) :: sirius, noarch, del_sout

 write(buf,'(2(A,I0))') 'dalton -gb ',MIN(mem,16),' -omp ',nproc
 if(sirius) buf = TRIM(buf)//" -put ""SIRIUS.RST"""
 if(noarch) buf = TRIM(buf)//' -noarch'
 buf = TRIM(buf)//' -ow '//TRIM(proname)//' >'//TRIM(proname)//".sout 2>&1"
 write(6,'(A)') '$'//TRIM(buf)
 i = system(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_cas: Dalton job failed.'
  write(6,'(A)') 'Please open files '//TRIM(proname)//'.out and '//TRIM(proname)&
                  //'.sout and check why.'
  stop
 end if

 if(del_sout) call delete_file(TRIM(proname)//'.sout')
end subroutine submit_dalton_job

subroutine submit_dalton_icss_job(proname, mem, nproc, nfile)
 implicit none
 integer :: i, n, fid, system
 integer, intent(in) :: mem, nproc, nfile
 character(len=240) :: buf, makefile, molname, dalname
 character(len=240), intent(in) :: proname

 if(nproc < 2) then
  write(6,'(/,A)') 'ERROR in subroutine submit_dalton_icss_job: the ICSS job &
                   &requires at least 2 CPU cores.'
  stop
 end if

 makefile = TRIM(proname)//'.make'
 open(newunit=fid,file=TRIM(makefile),status='replace')
 write(fid,'(A,/)') '# Makefile generated by MOKIT to run ICSS batch jobs'

 if(MOD(nproc,4) == 0) then
  write(fid,'(A,I0)') 'mem = ', MIN(16,mem/4)
  n = nproc/4
  write(fid,'(A,I0)') 'np = ', n
 else if(MOD(nproc,3) == 0) then
  write(fid,'(A,I0)') 'mem = ', MIN(16,mem/3)
  n = nproc/3
  write(fid,'(A,I0)') 'np = ', n
 else
  write(fid,'(A,I0)') 'mem = ', MIN(16,CEILING(DBLE(mem)*0.5d0))
  n = nproc/2
  write(fid,'(A,I0)') 'np = ', n
 end if

 write(fid,'(A)') 'pro = '//TRIM(proname)
 write(fid,'(/,A)',advance='no') 'all:'
 do i = 1, nfile, 1
  write(fid,'(A,I0)',advance='no') ' w',i
  if(MOD(i,15) == 0) write(fid,'(A,/,4X)',advance='no') ' \'
 end do ! for i
 write(fid,'(/)',advance='no')

 do i = 1, nfile, 1
  write(fid,'(/,2(A,I0),A)') 'w',i,':'
  write(molname,'(A,I5.5,A)') TRIM(proname),i,'.sout'
  write(fid,'(A1,2(A,I5.5),A)') ACHAR(9),"dalton -gb ${mem} -omp ${np} -put &
   &""SIRIUS.RST"" -noarch -ow ${pro}",i,' >${pro}',i,".sout 2>&1"
 end do ! for i
 close(fid)

 write(buf,'(A,I0,A)') 'make -j',n,' -f '//TRIM(makefile)//' all'
 write(6,'(A)') '$'//TRIM(buf)
 i = system(TRIM(buf))
 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine submit_dalton_icss_job: failed to run&
                & Dalton ICSS jobs.'
  stop
 end if

 do i = 1, nfile, 1
  write(molname,'(A,I5.5,A)') TRIM(proname), i, '.mol'
  write(dalname,'(A,I5.5,A)') TRIM(proname), i, '.dal'
  write(makefile,'(A,I5.5,A)') TRIM(proname), i, '.sout'
  call delete_files(3, [molname, dalname, makefile])
 end do ! for i
end subroutine submit_dalton_icss_job

! Cube format: https://gaussian.com/cubegen; http://sobereva.com/125
subroutine gen_icss_cub(proname, nfile)
 use mol, only: natom, coor, elem
 use fch_content, only: elem2nuc, Bohr_const
 use icss_param, only: ngrid, interval, origin
 implicit none
 integer :: i, j, k, tot_ngrid, fid, cubid
 integer, intent(in) :: nfile
 real(kind=8) :: intv
 real(kind=8), allocatable :: iso(:,:,:)
 character(len=240) :: buf, outname, cubname, bqfile
 character(len=240), intent(in) :: proname

 write(6,'(A)') 'Generating ICSS cube...'
 cubname = TRIM(proname)//'_ICSS.cub'
 bqfile = TRIM(proname)//'_ICSS.txt'
 open(newunit=cubid,file=TRIM(cubname),status='replace')
 write(cubid,'(A)') ' ICSS cube file generated by MOKIT'

 tot_ngrid = PRODUCT(ngrid)
 write(cubid,'(1X,I0,A)') tot_ngrid, ' grid points in total'

 write(cubid,'(I5,3(1X,F11.6),A)') natom, origin/Bohr_const, '    1'
 intv = interval/Bohr_const
 write(cubid,'(I5,1X,F11.6,A)') ngrid(1),intv,'    0.000000    0.000000'
 write(cubid,'(I5,A,1X,F11.6,A)') ngrid(2),'    0.000000',intv,'    0.000000'
 write(cubid,'(I5,A,1X,F11.6)') ngrid(3),'    0.000000    0.000000',intv

 do i = 1, natom, 1
  if(elem(i) == 'Bq') exit
  j = elem2nuc(elem(i))
  write(cubid,'(2I5,A,3(1X,F11.6))') j, j, '.000000', coor(:,i)/Bohr_const
 end do ! for i

 allocate(iso(ngrid(3),ngrid(2),ngrid(1)), source=0d0)
 call merge_shielding_into_one_file(proname, bqfile, nfile)
 open(newunit=fid,file=TRIM(bqfile),status='old',position='rewind')
 do i = 1, ngrid(1), 1
  do j = 1, ngrid(2), 1
   do k = 1, ngrid(3), 1
    read(fid,*) iso(k,j,i)
   end do ! for k
  end do ! for j
 end do ! for i

 close(fid,status='delete')
 do i = 1, ngrid(1), 1
  do j = 1, ngrid(2), 1
   write(cubid,'(6ES13.5)') iso(:,j,i)
  end do ! for j
 end do ! for i

 deallocate(iso)
 close(cubid)
 write(6,'(A)') 'ICSS cube generated.'
end subroutine gen_icss_cub

subroutine merge_shielding_into_one_file(proname, bqfile, nfile)
 implicit none
 integer :: i, fid1, fid2
 integer, intent(in) :: nfile
 character(len=240) :: buf, outname
 character(len=240), intent(in) :: proname, bqfile

 open(newunit=fid1,file=TRIM(bqfile),status='replace')

 do i = 1, nfile, 1
  write(outname,'(A,I5.5,A)') TRIM(proname),i,'.out'
  open(newunit=fid2,file=TRIM(outname),status='old',position='rewind')

  do while(.true.)
   read(fid2,'(A)') buf
   if(buf(1:5) == '@1 Bq') exit
  end do ! for while

  BACKSPACE(fid2)
  do while(.true.)
   read(fid2,'(A)') buf
   if(LEN_TRIM(buf) == 0) exit
   write(fid1,'(A)') buf(8:19)
  end do ! for while

  close(fid2)
 end do ! for i

 close(fid1)
end subroutine merge_shielding_into_one_file

subroutine read_shieldings_from_dalton_out(outname)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(28:57) == 'Summary of chemical shieldings') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(A)') 'ERROR in subroutine read_shieldings_from_dalton_out: no chem&
                 &ical shieldings found in file '//TRIM(outname)
  stop
 end if

 BACKSPACE(fid)
 BACKSPACE(fid)
 write(6,'(/,A)') 'Chemical shieldings copied from Dalton output:'
 do while(.true.)
  read(fid,'(A)') buf
  if(index(buf,'Intera')>0 .or. index(buf,'intera')>0) exit
  write(6,'(A)') TRIM(buf)
 end do ! for while

 close(fid)
end subroutine read_shieldings_from_dalton_out

