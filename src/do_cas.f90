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
 use mr_keyword, only: mem, nproc, casci, dmrgci, casscf, dmrgscf, dmrg_no, &
  ist, hf_fch, datname, nacte_wish, nacto_wish, gvb, casnofch, casci_prog, &
  casscf_prog, dmrgci_prog, dmrgscf_prog, gau_path, gms_path, molcas_omp, &
  molcas_path, orca_path, gms_scr_path, molpro_path, bdf_path, psi4_path, bgchg,&
  chgname, casscf_force, check_gms_path, prt_strategy, RI, nmr, ICSS, on_thres,&
  iroot, xmult, dyn_corr
 use mol, only: nbf, nif, npair, nopen, npair0, ndb, casci_e, casscf_e, nacta, &
  nactb, nacto, nacte, gvb_e, ptchg_e, nuc_pt_e, natom, grad
 use util_wrapper, only: bas_fch2py_wrap, formchk, unfchk, gbw2mkl, mkl2gbw, &
  fch2inp_wrap, mkl2fch_wrap, fch2inporb_wrap, orb2fch_wrap
 implicit none
 integer :: i, j, idx1, idx2, nvir, nfile, system, RENAME
 real(kind=8) :: unpaired_e ! unpaired electrons
 real(kind=8) :: e(2)       ! e(1) is CASCI enery, e(2) is CASSCF energy
 real(kind=8), allocatable :: noon(:)
 character(len=10) :: cas_prog = ' '
 character(len=24) :: data_string = ' '
 character(len=240) :: fchname, pyname, inpname, outname, proname, mklname, &
                       orbname, xmlname, gradname
 character(len=480) :: buf
 logical, intent(in) :: scf
 logical :: alive1, alive2

 if(scf) then
  if((.not. casscf) .and. (.not.dmrgscf)) return
 else
  if((.not. casci) .and. (.not.dmrgci)) return
 end if
 write(6,'(//,A)') 'Enter subroutine do_cas...'

 if(ist == 5) then
  write(6,'(A)') 'Radical index for input NOs:'
  call calc_unpaired_from_fch(hf_fch, 3, .false., unpaired_e)
  ! read nbf, nif, nopen, nacto, ... variables from NO .fch(k) file
  call read_no_info_from_fch(hf_fch, on_thres, nbf, nif, ndb, nopen, nacta, &
                             nactb, nacto, nacte)
  i = nacte; j = nacto
  npair0 = nactb; npair = npair0
 else
  i = 2*npair0 + nopen; j = i
 end if

 alive1 = (nacte_wish>0 .and. nacte_wish/=i)
 alive2 = (nacto_wish>0 .and. nacto_wish/=j)

 if(alive1 .or. alive2) then
  write(6,'(4(A,I0),A)') 'Warning: AutoMR recommends CAS(',i,'e,',j,'o), but&
   & you specify CAS(',nacte_wish,'e,',nacto_wish,'o). Trying to fulfill...'

  ! check the odevity of nacte_wish, in case that the user requires nonsense
  ! number of active electrons
  if(MOD(nacte_wish-nopen,2) /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine do_cas: wrong active space specified.'
   write(6,'(3(A,I0),A)') 'Nopen=',nopen,'. Incompatible with CAS(',nacte_wish,&
                          'e,',nacto_wish,'o)'
   stop
  end if

  if(ist == 5) then
   call prt_active_space_warn(nacte_wish, nacto_wish, nacte, nacto)
   i = nacte_wish; j = nacto_wish
   npair0 = (i - nopen)/2; npair = npair0
   ndb = ndb + nactb - npair
   nacta = npair0 + nopen; nactb = npair0
   nacto = nacto_wish;     nacte = nacto_wish
   write(6,'(A)') 'OK, fulfilled.'

  else ! ist /= 5
   if(2*npair+nopen < nacte_wish) then
    write(6,'(/,A)') 'ERROR in subroutine do_cas: too large space specified. Ca&
                     &nnot be fulfilled.'
    write(6,'(2(A,I0))') '2*npair+nopen=',2*npair+nopen, ', nacte_wish=',nacte_wish
    stop
   else ! 2*npair+nopen >= nacte_wish
    write(6,'(A)') 'OK, fulfilled.'
    npair0 = (nacte_wish - nopen)/2
    i = 2*npair0 + nopen; j = i
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
 write(6,'(A)',advance='no') TRIM(data_string)

 idx1 = ndb + npair - npair0 + 1
 idx2 = idx1 + 2*npair0 + nopen - 1
 nvir = nif - (idx1-1) - 2*npair0 - nopen
 write(6,'(A,2(I0,A))') '(',i,'e,',j,'o) using program '//TRIM(cas_prog)

 if(i*j == 0) then
  write(6,'(/,A)') REPEAT('-', 79)
  write(6,'(A)') 'There is zero active orbital/electron. AutoMR terminated.'
  write(6,'(/,A)') 'The reason is this molecule has little multi-configurationa&
                   &l or multi-reference'
  write(6,'(A)') 'character. This molecule can be well described by single refe&
                 &rence methods, e.g.'
  write(6,'(A)') 'MP2, CCSD(T). Thus no need for multi-reference computation. B&
                 &ut if you still want'
  write(6,'(A)') 'to do it, you can manually specify the size of active space i&
                 &n .gjf file. For ex-'
  write(6,'(A)') 'ample, CASSCF(4,4) for water(H2O), or CASSCF(8,8) for methane&
                 &(CH4). The maximum'
  write(6,'(A)') 'number of active orbitals is 2*npair. You can find npair in G&
                 &VB computations above.'
  write(6,'(A,I0)') 'For this molecule, npair=', npair
  write(6,'(A)') REPEAT('-', 79)
  stop
 end if

 write(6,'(4(A,I0),A,L1)') 'doubly_occ=', idx1-1, ', nvir=', nvir, ', Root=', &
                           iroot, ', Xmult=', xmult, ', RIJK=', RI
 write(6,'(2(A,I0))') 'No. of active alpha/beta e = ', nacta,'/',nactb

 if(nopen+2*npair0 > 15) then
  if(nmr) then
   write(6,'(/,A)') 'ERROR in subroutine do_cas: DMRG invoked, but DMRG-GIAO&
                   & is not supported.'
   stop
  end if
  if(scf) then
   casscf = .false.
   dmrgscf = .true.
   cas_prog = dmrgscf_prog
   write(6,'(A)') 'Remark: CASSCF is switched to DMRG-CASSCF due to active spac&
                  &e larger than (15,15).'
   if(casscf_prog /= 'pyscf') then
    write(6,'(/,A)') 'ERROR in subroutine do_cas: DMRGSCF required. But CASSCF_&
                     &prog='//TRIM(casscf_prog)//'.'
    stop
   end if

  else ! not scf
   casci = .false.
   dmrgci = .true.
   cas_prog = dmrgci_prog
   write(6,'(A)') 'Remark: CASCI is switched to DMRG-CASCI due to active space &
                  &larger than (15,15).'
   if(casci_prog /= 'pyscf') then
    write(6,'(/,A)') 'ERROR in subroutine do_cas: DMRGCI required. But CASCI_pr&
                     &og='//TRIM(casci_prog)//'.'
    stop
   end if
  end if

  write(6,'(A)') 'Strategy updated:'
  call prt_strategy()
 end if

 if(ist<1 .or. ist>6) then
  write(6,'(/,A)') 'ERROR in subroutine do_cas: invalid ist.'
  write(6,'(A,I0)') 'Allowed values are 1~6. But got ist=', ist
  stop
 end if

 if((dmrgci .or. dmrgscf) .and. TRIM(cas_prog)/='pyscf') then
  write(6,'(/,A)') 'ERROR in subroutine do_cas: DMRG-CASCI/CASSCF calculation i&
                   &s only supported by'
  write(6,'(A)') 'PySCF. Wrong CASCI_prog/CASSCF_prog='//TRIM(cas_prog)
  stop
 end if

 select case(ist)
 case(1,3,6)
  i = INDEX(datname, '.dat', back=.true.)
  fchname = datname(1:i-1)//'.fch'
 case(2) ! UHF -> UNO -> CASCI/CASSCF
  i = INDEX(hf_fch, '.fch', back=.true.)
  fchname = hf_fch(1:i-1)//'_uno.fch'
  pyname = hf_fch(1:i-1)//'_uno.py'
  inpname = hf_fch(1:i-1)//'_uno.py2'
  i = RENAME(TRIM(pyname), TRIM(inpname))
  ! bas_fch2py will generate file '_uno.py', so we need to rename it to another filename
 case(5)
  fchname = hf_fch
 end select

 proname = ' '
 i = INDEX(hf_fch, '.fch', back=.true.)
 select case(ist)
 case(1,3,6)
  if(scf) then
   write(proname,'(A,I0,A)') hf_fch(1:i-1)//'_gvb', npair, '_CASSCF'
  else
   write(proname,'(A,I0,A)') hf_fch(1:i-1)//'_gvb', npair, '_CASCI'
  end if
 case(2)
  if(scf) then
   write(proname,'(A)') hf_fch(1:i-1)//'_uno2CASSCF'
  else
   write(proname,'(A)') hf_fch(1:i-1)//'_uno2CASCI'
  end if
 case(5)
  if(scf) then
   proname = hf_fch(1:i-1)//'_CASSCF'
  else
   proname = hf_fch(1:i-1)//'_CASCI'
  end if
 end select
 casnofch = TRIM(proname)//'_NO.fch'

 select case(TRIM(cas_prog))
 case('pyscf')
  inpname = TRIM(proname)//'.py'
  outname = TRIM(proname)//'.out'
  call bas_fch2py_wrap(fchname, .false., inpname)
  call prt_cas_script_into_py(inpname, scf)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  if(casscf_force) call add_force_key2py_script(mem, inpname, .false.)
  call submit_pyscf_job(inpname)

 case('gaussian')
  call check_exe_exist(gau_path)

  inpname = TRIM(proname)//'.gjf'
  outname = TRIM(proname)//'.log'
  mklname = TRIM(proname)//'.chk'
  call prt_cas_gjf(inpname, nacto, nacte, scf, casscf_force)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call unfchk(fchname, mklname)
  call submit_gau_job(gau_path, inpname, .true.)
  call formchk(mklname, casnofch)
  call modify_IROHF_in_fch(casnofch, 0)

 case('gamess')
  call check_gms_path()
  call fch2inp_wrap(fchname, .false., 0, 0)
  i = INDEX(fchname, '.fch', back=.true.)
  outname = fchname(1:i-1)//'.inp'
  inpname = TRIM(proname)//'.inp'
  i = RENAME(TRIM(outname), TRIM(inpname))
  outname = TRIM(proname)//'.gms'
  call prt_cas_gms_inp(inpname, idx1-1, scf)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  if(casscf_force) call add_force_key2gms_inp(inpname)

  call submit_gms_job(gms_path, gms_scr_path, inpname, nproc)

  ! make a copy of the .fch file to save NOs
  if(ist /= 2) then ! datname is a GVB job .dat file
   i = INDEX(datname, '.dat', back=.true.)
   inpname = datname(1:i-1)//'.fch'
  else              ! no GVB job
   i = INDEX(hf_fch, '.fch')
   inpname = hf_fch(1:i-1)//'_uno.fch'
  end if
  call copy_file(inpname, casnofch, .false.)

  ! transfer CASSCF pseudo-canonical MOs from .dat to .fch
  i = INDEX(casnofch, '_NO.fch', back=.true.)
  datname = casnofch(1:i-1)//'.dat'
  if(scf) then
   buf = 'dat2fch '//TRIM(datname)//' '//TRIM(casnofch)
   write(6,'(A)') '$'//TRIM(buf)
   i = SYSTEM(TRIM(buf))
   write(buf,'(A,2(1X,I0))') 'extract_noon2fch '//TRIM(outname)//' '//&
                              TRIM(casnofch), idx1, idx2
   write(6,'(A)') '$'//TRIM(buf)
   i = SYSTEM(TRIM(buf))
  end if

  ! transfer NOs from .dat to .fch
  write(buf,'(A,I0)') 'dat2fch '//TRIM(datname)//' '//TRIM(casnofch)//' -no ',idx2
  write(6,'(A)') '$'//TRIM(buf)
  i = SYSTEM(TRIM(buf))
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine do_cas: failed to call utility dat2fch.'
   write(6,'(A)') 'Related files: '//TRIM(datname)//', '//TRIM(casnofch)//'.'
   stop
  end if

 case('openmolcas')
  call check_exe_exist(molcas_path)
  inpname = TRIM(proname)//'.input'
  outname = TRIM(proname)//'.out'
  call fch2inporb_wrap(fchname, .false., inpname)
  write(orbname,'(A,I0)') TRIM(proname)//'.RasOrb.', iroot+1
  call prt_cas_molcas_inp(inpname, scf)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  if(casscf_force) i = SYSTEM("echo '&ALASKA' >> "//TRIM(inpname))
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)

  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  call orb2fch_wrap(orbname, casnofch, .true.) ! transfer NOs from .dat to .fch
  ! OpenMolcas CASSCF NOs may not be in ascending order, sort them
  call sort_no_by_noon(casnofch, idx1, idx2)

 case('orca')
  call check_exe_exist(orca_path)
  i = SYSTEM('fch2mkl '//TRIM(fchname))

  i = INDEX(fchname, '.fch', back=.true.)
  pyname  = fchname(1:i-1)//'_o.inp'
  orbname = fchname(1:i-1)//'_o.mkl'
  inpname = TRIM(proname)//'.inp'
  mklname = TRIM(proname)//'.mkl'
  outname = TRIM(proname)//'.out'
  gradname = TRIM(proname)//'.engrad'
  i = RENAME(TRIM(orbname),TRIM(mklname))
  i = RENAME(TRIM(pyname),TRIM(inpname))

  call prt_cas_orca_inp(inpname, scf)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(mklname)
  if(casscf_force) call add_force_key2orca_inp(inpname)

  call submit_orca_job(orca_path, inpname, .true.)
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
  call mkl2fch_wrap(mklname, casnofch, 1)

 case('molpro')
  call check_exe_exist(molpro_path)

  i = SYSTEM('fch2com '//TRIM(fchname)) ! generate .com and .txt
  i = INDEX(fchname, '.fch', back=.true.)
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
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  call submit_molpro_job(inpname, mem, nproc)
  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  i = SYSTEM('xml2fch '//TRIM(xmlname)//' '//TRIM(casnofch)//' -no')

 case('bdf')
  call check_exe_exist(bdf_path)

  i = SYSTEM('fch2bdf '//TRIM(fchname)//' -no')
  i = INDEX(fchname, '.fch', back=.true.)
  mklname = fchname(1:i-1)//'_bdf.inp'
  pyname  = fchname(1:i-1)//'_bdf.inporb'
  inpname = TRIM(proname)//'.inp'
  xmlname = TRIM(proname)//'.inporb'
  orbname = TRIM(proname)//'.casorb'
  outname = TRIM(proname)//'.out'
  i = RENAME(TRIM(mklname), TRIM(inpname))
  i = RENAME(TRIM(pyname), TRIM(xmlname))
  call prt_cas_bdf_inp(inpname, scf, casscf_force)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  write(6,'(A)') '$$BDF '//TRIM(proname)
  i = SYSTEM(TRIM(bdf_path)//' '//TRIM(proname))
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine do_cas: BDF CASCI/CASSCF job failed.'
   stop
  end if
  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  i = SYSTEM('bdf2fch '//TRIM(orbname)//' '//TRIM(casnofch)//' -no')

 case('psi4')
  i = SYSTEM('fch2psi '//TRIM(fchname))
  i = INDEX(fchname, '.fch', back=.true.)
  mklname = fchname(1:i-1)//'_psi.inp'
  pyname  = fchname(1:i-1)//'_psi.A'
  inpname = TRIM(proname)//'.inp'
  xmlname = TRIM(proname)//'.A'
  outname = TRIM(proname)//'.out'
  i = RENAME(TRIM(mklname), TRIM(inpname))
  i = RENAME(TRIM(pyname), TRIM(xmlname))
  call prt_cas_psi4_inp(inpname, scf)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  call submit_psi4_job(psi4_path, inpname, nproc)
  write(buf,'(A,2(1X,I0))') 'extract_noon2fch '//TRIM(outname)//' '//&
                             TRIM(casnofch), idx1, idx2
  i = SYSTEM(TRIM(buf))
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine do_cas: failed to call utility extract&
                    &_noon2fch.'
   stop
  end if

 case('dalton')
  write(6,'(/,A)') 'Warning: Currently automr only supports the OpenMP version &
                   &of Dalton. If you'
  write(6,'(A,/)') 'are using MPI version of Dalton, you should kill the job.'
  i = SYSTEM('fch2dal '//TRIM(fchname))
  i = INDEX(fchname, '.fch', back=.true.)
  mklname = fchname(1:i-1)//'.dal'
  pyname  = fchname(1:i-1)//'.mol'
  inpname = TRIM(proname)//'.dal'
  xmlname = TRIM(proname)//'.mol'
  orbname = TRIM(proname)//'.MOPUN'
  outname = TRIM(proname)//'.out'
  i = RENAME(TRIM(mklname), TRIM(inpname))
  i = RENAME(TRIM(pyname), TRIM(xmlname))
  call prt_cas_dalton_inp(inpname, scf, casscf_force)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
                                             ! sirius, noarch, del_sout
  call submit_dalton_job(proname, mem, nproc, .false., .false., .false.)

  ! untar/unzip the compressed package to get DALTON.MOPUN
  i = SYSTEM('tar -xpf '//TRIM(proname)//'.tar.gz DALTON.MOPUN')

  ! make a copy to save NOs
  call copy_file(fchname, casnofch, .false.)

  ! transfer NOs (and NOONs, if any) from DALTON.MOPUN to .fch file
  buf = 'dal2fch DALTON.MOPUN '//TRIM(casnofch)
  if(.not. casscf_force) buf = TRIM(buf)//' -no'
  ! when force=.T., NOONs are not in DALTON.MOPUN, so we only transfer NOs here
  i = SYSTEM(TRIM(buf))
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine do_cas: failed to call utility dal2fch.'
   write(6,'(A)') 'Please open files DALTON.MOPUN and '//TRIM(casnofch)//&
                 &' and check.'
   stop
  end if

  ! Now transfer NOONs when force=.T.
  if(casscf_force) then
   allocate(noon(nif), source=0d0)
   call read_on_from_dalton_out(outname, idx2, noon(1:idx2))
   call write_eigenvalues_to_fch(casnofch, nif, 'a', noon, .true.)
   deallocate(noon)
  end if
  call delete_file('DALTON.MOPUN')

 case default
  write(6,'(/,A)') 'ERROR in subroutine do_cas: allowed programs are Gaussian,&
                   & GAMESS, PySCF,'
  write(6,'(A)') 'OpenMolcas, ORCA, Molpro, BDF, PSI4 and Dalton. But got CAS_&
                 &prog='//TRIM(cas_prog)
  stop
 end select

 ! generate CASCI/CASSCF density from NOs and NOONs
 ! Note: density in .fch of Gaussian CASSCF job is wrong (Gaussian bug), I have
 !  to re-generate density
 select case(TRIM(cas_prog))
 case('pyscf','orca') ! do nothing, density is already correct
 case default
  call update_density_using_no_and_on(casnofch)
 end select

 if(ist == 2) then
  i = INDEX(hf_fch, '.fch', back=.true.)
  pyname = hf_fch(1:i-1)//'_uno.py'
  inpname = hf_fch(1:i-1)//'_uno.py2'
  i = RENAME(TRIM(inpname), TRIM(pyname)) ! rename it back
 end if

 ! read energy, check convergence and check spin
 call read_cas_energy_from_output(cas_prog, outname, e, scf, nacta-nactb, &
                                  (dmrgci.or.dmrgscf), ptchg_e, nuc_pt_e)

 if(gvb .and. 2*npair+nopen==nacto .and. iroot==0 .and. e(1)-gvb_e>2D-6) then
  write(6,'(/,A)') 'ERROR in subroutine do_cas: active space of GVB and CAS&
                  & are equal, but CASCI/CASSCF'
  write(6,'(A)') 'energy is higher than that of GVB. This is probably due to:&
                & (1) CASCI stucks in a higher'
  write(6,'(A)') 'energy local minimum or not pure spin state; (2) GVB MOs are&
                & disordered (GAMESS bug).'
  stop
 end if

 casci_e = e(1)
 write(6,'(/,A,F18.8,1X,A4)') 'E(CASCI)  = ', e(1), 'a.u.'
 if(scf) then
  casscf_e = e(2)
  write(6,'(A,F18.8,1X,A4)') 'E(CASSCF) = ', e(2), 'a.u.'
 end if

 if((dmrgci .and. dmrg_no) .or. (.not. dmrgci)) then
  call calc_unpaired_from_fch(casnofch, 3, .false., unpaired_e)
 end if

 if(casscf_force) then
  allocate(grad(3*natom))
  ! rename outname to use subroutine read_grad_from_output
  select case(TRIM(cas_prog))
  case('gamess')
   outname = datname
  case('orca')
   outname = gradname
  end select
  call read_grad_from_output(cas_prog, outname, natom, grad)
 end if

 if(nmr) then
  call prt_cas_dalton_nmr_inp(casnofch, scf, ICSS, iroot, nfile)
  inpname = TRIM(proname)//'_NMR'
                                             ! sirius, noarch, del_sout
  call submit_dalton_job(inpname, mem, nproc, .false., .false., .false.)
  inpname = TRIM(proname)//'_NMR.out'
  call read_shieldings_from_dalton_out(inpname)
 end if

 if(ICSS) then
  call submit_dalton_icss_job(proname, mem, nproc, nfile)
  call gen_icss_cub(proname, nfile)
 end if

 if(.not. dyn_corr) call delete_file('ss-cas.txt')
 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_cas at '//TRIM(data_string)
end subroutine do_cas

! print CASCI/DMRG-CASCI or CASSCF/DMRG-CASSCF script into a given .py file
subroutine prt_cas_script_into_py(pyname, scf)
 use mr_keyword, only: mem, nproc, dmrgci, dmrgscf, dmrg_no, block_mpi, RI, &
  RIJK_bas, iroot, casnofch
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, pyname1, cmofch
 character(len=240), intent(in) :: pyname
 logical, intent(in) :: scf
 logical :: dmrg

 dmrg = (dmrgci .or. dmrgscf)
 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
 pyname1 = TRIM(pyname)//'.t'
 i = INDEX(pyname, '.py', back=.true.)

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(dmrg) then
  write(fid2,'(A)') 'from pyscf import mcscf, dmrgscf, lib'
 else
  write(fid2,'(A)') 'from pyscf import mcscf, lib'
 end if
 write(fid2,'(A)') 'from mokit.lib.py2fch import py2fch'
 write(fid2,'(A)') 'from shutil import copyfile'
 write(fid2,'(A,/)') 'import numpy as np'
 if(dmrg) then
  write(fid2,'(A)',advance='no') "dmrgscf.settings.MPIPREFIX = '"
  if(block_mpi) write(fid2,'(A,I0)',advance='no') "mpirun -n ", nproc
  write(fid2,'(A)') "'"
 end if
 write(fid2,'(A,I0)') 'nproc = ', nproc
 write(fid2,'(A,/)') 'lib.num_threads(nproc)'

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

 ! For CASCI/CASSCF, always set mc.natorb = True.
 ! For DMRG-CASCI/CASSCF, a two-step job is employed here: 
 !  1) use input MOs to perform (DMRG-)CASCI/CASSCF;
 !  2) perform a CASCI to generate NOs.

 ! For DMRG-CASCI/CASSCF, both the original MOs and NOs will be saved/punched.
 ! Since DMRG is not strictly invariant to unitary rotations of active orbitals,
 ! I hope the original MOs to be used in DMRG-PDFT/NEVPT2/CASPT2 computations.
 ! NOs are usually more delocalized than original MOs.

 ! mem*500 is in fact mem*1000/2. The mc.max_memory and fcisolver.max_memory seem
 ! not to share common memory, they used two memory, so I have to divide them
 if(scf) then
  if(iroot == 0) then
   call prt_gs_casscf_kywrd_py(fid2, RIJK_bas1) ! print ground state CASSCF keywords
  else
   call prt_es_casscf_kywrd_py(fid2, RIJK_bas1) ! print excited state CASSCF keywords
  end if
 else ! (DMRG-)CASCI
  call prt_casci_kywrd_py(fid2, RIJK_bas1, .false.)
 end if

 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.kernel()'

 i = INDEX(casnofch, '_NO', back=.true.)
 cmofch = casnofch(1:i)//'CMO.fch'
 if(dmrg) write(fid2,'(/,A)') "cmofch = '"//TRIM(cmofch)//"'"
 if(dmrgci) then
  write(fid2,'(/,A)') '# backup original MOs to cmofch file'
  write(fid2,'(A)') 'copyfile(hf_fch, cmofch)'
 else if(dmrgscf) then
  write(fid2,'(/,A)') '# save conveged MOs into .fch file'
  write(fid2,'(A)') 'copyfile(hf_fch, cmofch)'
  write(fid2,'(A)') 'noon = np.zeros(nif)'
  write(fid2,'(A)') "py2fch(cmofch, nbf, nif, mc.mo_coeff, 'a', noon, False, False)"
 end if

 if(dmrg) then
  if(dmrg_no) then ! perform DMRG-CASCI to generate NOs
   call prt_casci_kywrd_py(fid2, RIJK_bas1, .true.)
  else
   casnofch = cmofch
  end if
 end if

 if((dmrg .and. dmrg_no) .or. (.not. dmrg)) then
  write(fid2,'(/,A)') '# save NOs into .fch file'
  write(fid2,'(A)') "casnofch = '"//TRIM(casnofch)//"'"
  write(fid2,'(A)') "copyfile(hf_fch, casnofch)"
  write(fid2,'(A)') "py2fch(casnofch, nbf, nif, mc.mo_coeff, 'a', mc.mo_occ, Tr&
                    &ue, True)"
  ! Note: mc.mo_occ is only valid for PySCF >= 1.7.4
 end if

 close(fid2)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
end subroutine prt_cas_script_into_py

! print a CASCI or CASSCF gjf file
subroutine prt_cas_gjf(gjfname, nacto, nacte, scf, force)
 use mr_keyword, only: mem, nproc, dkh2_or_x2c, iroot
 implicit none
 integer :: i, fid
 integer, intent(in) :: nacto, nacte
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: scf, force

 i = INDEX(gjfname, '.gjf', back=.true.)
 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(2(A,I0))',advance='no') '#p CASSCF(',nacte,',',nacto
 if(iroot > 0) write(fid,'(A,I0)',advance='no') ',qc,nroot=',iroot+1
 write(fid,'(A)',advance='no') ') chkbasis nosymm guess=read geom=allcheck'

 if(dkh2_or_x2c) then
  write(fid,'(A)',advance='no') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A)',advance='no') ' int=nobasistransform'
 end if
 if(force) write(fid,'(A)',advance='no') ' force'

 if(scf) then ! CASSCF
  write(fid,'(A)') ' scf(maxcycle=200)'
 else         ! CASCI
  write(fid,'(A)') ' scf(maxcycle=-2)'
  ! to obtain CASCI NOs, we need to use -2, since -1 only calculates CASCI energy
 end if

 write(fid,'(/,A)') '--Link1--'
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

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 if(scf) then   ! CASSCF
  write(fid2,'(A)',advance='no') ' $CONTRL SCFTYP=MCSCF RUNTYP=ENERGY ICHARG='
 else           ! CASCI
  write(fid2,'(A)',advance='no') ' $CONTRL SCFTYP=NONE CITYP=ALDET RUNTYP=ENERG&
                                 &Y ICHARG='
 end if
 write(fid2,'(2(I0,A))') charge, ' MULT=', mult, ' NOSYM=1'

 write(fid2,'(A)',advance='no') '  ICUT=11'
 if(dkh2_or_x2c) write(fid2,'(A)',advance='no') ' RELWFN=DK'

 if(.not. cart) then
  write(fid2,'(A)') ' ISPHER=1 $END'
 else
  write(fid2,'(A)') ' $END'
 end if

 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',CEILING(DBLE(mem*125)/DBLE(nproc)),&
                        ' $END'
 if(scf) then   ! CASSCF
  write(fid2,'(A)',advance='no') ' $DET'
 else
  write(fid2,'(A)',advance='no') ' $CIDET'
 end if

 write(fid2,'(3(A,I0),A)',advance='no') ' NCORE=', ncore, ' NELS=', nacte, &
                                        ' NACT=', nacto, ' ITERMX=500'
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
end subroutine prt_cas_gms_inp

! print CASCI/CASSCF keywords in to a given (Open)Molcas input file
subroutine prt_cas_molcas_inp(inpname, scf)
 use mr_keyword, only: dmrgci, dmrgscf, RI, RIJK_bas, mokit_root
 implicit none
 integer :: i, j, fid1, fid2, RENAME, system
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf
 logical :: dmrg

 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
 dmrg = (dmrgci .or. dmrgscf)

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == 'noCD') exit
  j = INDEX(buf, '...')
  if(j > 0) then
   if(RI) then
    j = INDEX(buf, '.')
    write(fid2,'(A)') buf(1:j)//TRIM(RIJK_bas1)//'..'//TRIM(buf(j+3:))
   else
    write(fid2,'(A)') TRIM(buf)
   end if
  else
   write(fid2,'(A)') TRIM(buf)
  end if
 end do

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine prt_cas_molcas_inp: no 'SEWARD' found i&
                   &n file "//TRIM(inpname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 if(RI) write(fid2,'(A)') 'RIJK'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4)=="&SCF" .or. buf(1:7)=="&RASSCF") exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine prt_cas_molcas_inp: '&SCF'  or '&RASSCF&
                   &' not found in"
  write(6,'(A)') 'file '//TRIM(inpname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 call prt_molcas_cas_para(fid2, dmrg, .false., .false., (.not.scf), inpname)

 write(fid2,'(/)',advance='no')
 close(fid2)
 close(fid1,status='delete')
 i = RENAME(TRIM(inpname1), TRIM(inpname))

 ! if RIJK is on, we need to generate the fitting basis set file for OpenMolcas
 if(RI) then
  i = SYSTEM('cp '//TRIM(mokit_root)//'/mokit/basis/'//TRIM(RIJK_bas1)//' .')
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine prt_cas_molcas_inp: failed to copy fil&
                    &e from'
   write(6,'(A)') TRIM(mokit_root)//'/mokit/basis/'//TRIM(RIJK_bas1)//' to the&
                 & current directory.'
   stop
  end if
  i = SYSTEM('bas_gau2molcas '//TRIM(RIJK_bas1))
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine prt_cas_molcas_inp: failed to call uti&
                    &lity bas_gau2molcas.'
   write(6,'(A)') 'Did you forget to compile it?'
   stop
  end if

  call delete_file(RIJK_bas1)
  call upper(RIJK_bas1)
  i = SYSTEM('mv '//TRIM(RIJK_bas1)//' $MOLCAS/basis_library/jk_Basis/')
 end if
end subroutine prt_cas_molcas_inp

! print CASCI/CASSCF keywords in to a given ORCA .inp file
subroutine prt_cas_orca_inp(inpname, scf)
 use mol, only: nacte, nacto, mult
 use mr_keyword, only: mem, nproc, RI, RIJK_bas, hardwfn, crazywfn, iroot, xmult
 implicit none
 integer :: i, j, fid1, fid2, RENAME
 character(len=3), allocatable :: weight(:)
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf

 if(iroot > 19) then
  write(6,'(A)') 'ERROR in subroutine prt_cas_orca_inp: please contact the MOKI&
                 &T developers to'
  write(6,'(A)') 'modify the printing format.'
  stop
 end if
 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 read(fid1,'(A)') buf   ! skip nproc
 read(fid1,'(A)') buf   ! skip memory
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid2,'(A,I0)') '%maxcore ', CEILING(1d3*DBLE(mem)/DBLE(nproc))

 read(fid1,'(A)') buf   ! skip '!' line
 write(fid2,'(A)',advance='no') '!'
 if(.not. scf) write(fid2,'(A)',advance='no') ' noiter'
 ! RIJK in CASSCF must be combined with CONVentional
 if(RI) write(fid2,'(A)',advance='no') ' RIJK conv '//TRIM(RIJK_bas)
 write(fid2,'(A)') ' TightSCF'

 if(scf) then ! CASSCF
  write(fid2,'(A)') '%casscf'
  write(fid2,'(A,I0)') ' nel ', nacte
  write(fid2,'(A,I0)') ' norb ', nacto
  write(fid2,'(A)') ' maxiter 200'
  write(fid2,'(A)') ' ActOrbs NatOrbs'
  if(iroot > 0) then ! SS-CASSCF
   write(fid2,'(A,I0)') ' mult ', xmult
   j = 2
   if(xmult /= mult) j = 1
   write(fid2,'(A,I0)') ' nroots ', iroot+j+1
   allocate(weight(iroot+j))
   weight = '0.0'
   write(fid2,'(A)',advance='no') ' weights[0]='
   if(xmult == mult) then
    write(fid2,'(19A4)',advance='no') (weight(i)//',',i=1,iroot)
   else ! xmult /= mult
    if(iroot > 1) then
     write(fid2,'(19A4)',advance='no') (weight(i)//',',i=1,iroot-1)
    end if
   end if
   write(fid2,'(A)') '1.0,0.0,0.0'
   deallocate(weight)
  end if
  call prt_hard_or_crazy_casci_orca(fid2, hardwfn, crazywfn)
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
end subroutine prt_cas_orca_inp

! print CASCI/CASSCF keywords into a given Molpro input file
subroutine prt_cas_molpro_inp(inpname, scf, force)
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
  write(6,'(A)') 'ERROR in subroutine prt_cas_molpro_inp: wrong content found&
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

 write(fid,'(2(A,I0))',advance='no') ';closed,', nclosed, ';occ,', nocc
 ! Note: we need 'NoExtra' to completely turn off symmetry.
 ! Otherwise the CASCI energy is slightly different to that of other programs
 if(scf) then
  write(fid,'(A)') ';NoExtra;MaxIter,200}'
 else
  write(fid,'(A)') ';NoExtra;DONT,ORBITAL}'
 end if

 if(force) write(fid,'(A)') 'Forces'
 write(fid,'(A)') TRIM(put)
 close(fid)
end subroutine prt_cas_molpro_inp

! print CASCI/CASSCF keywords into a given BDF input file
subroutine prt_cas_bdf_inp(inpname, scf, force)
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
  write(6,'(A)') "ERROR in subroutine prt_cas_bdf_inp: '$SCF' not found in&
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
end subroutine prt_cas_bdf_inp

! print CASCI/CASSCF keywords into a given PSI4 input file
subroutine prt_cas_psi4_inp(inpname, scf)
 use mol, only: ndb, npair, npair0, nacte
 use mr_keyword, only: mem, hardwfn, crazywfn, RI, RIJK_bas
 implicit none
 integer :: i, nclosed, fid
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, casnofch
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf

 call modify_memory_in_psi4_inp(inpname, mem)
 i = INDEX(inpname, '.inp', back=.true.)
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

 ! The following line can be used with PSI4 1.3.2, but not work for >= 1.4
 !write(fid,'(A)') "fchk(cas_wfn,'"//TRIM(casnofch)//"')"

 ! So I make a workaround: copy CASSCF NOs into scf object, see
 ! http://forum.psicode.org/t/how-to-copy-orbitals-from-cas-wfn-to-scf-wfn/2541/5
 write(fid,'(A)') 'scf_wfn.Ca().copy(cas_wfn.Ca())'
 write(fid,'(A)') "fchk(scf_wfn,'"//TRIM(casnofch)//"')"
end subroutine prt_cas_psi4_inp

! print CASCI/CASSCF keywords into a given Dalton input file
subroutine prt_cas_dalton_inp(inpname, scf, force)
 use mol, only: mult, ndb, npair, npair0, nacto, nacte
 use mr_keyword, only: DKH2
 implicit none
 integer :: i, nclosed, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf, force

 inpname1 = TRIM(inpname)//'.t'
 nclosed = ndb + npair - npair0

 open(newunit=fid1,file=TRIM(inpname1),status='replace')
 write(fid1,'(A)') '**DALTON INPUT'
 if(DKH2) write(fid1,'(A)') '.DOUGLAS-KROLL'
 if(force) then
  write(fid1,'(A)') '.RUN PROPERTIES'
 else
  write(fid1,'(A)') '.RUN WAVE FUNCTIONS'
 end if
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
  write(6,'(/,A)') "ERROR in subroutine prt_cas_dalton_inp: no '**END OF INP' f&
                   &ound in file "//TRIM(inpname)
  close(fid1,status='delete')
  close(fid)
  stop
 end if

 if(force) write(fid1,'(A,/,A)') '**PROPERTIES', '.MOLGRA'
 write(fid1,'(A)') '**END OF INPUT'
 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_cas_dalton_inp

! print CASSCF NMR keywords into a given Dalton input file
subroutine prt_cas_dalton_nmr_inp(fchname, scf, ICSS, iroot, nfile)
 use mol, only: mult, ndb, npair, npair0, nacto, nacte
 use mr_keyword, only: DKH2
 use util_wrapper, only: fch2dal_wrap
 implicit none
 integer :: i, nclosed, fid, fid1, RENAME
 integer, intent(in) :: iroot
 ! 0 for ground state, 1 for the first excited state
 integer, intent(out) :: nfile
 character(len=240) :: buf, dalname, molname, olddal, oldmol
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: scf, ICSS

 i = INDEX(fchname, '.fch', back=.true.)
 olddal = fchname(1:i-1)//'.dal'
 oldmol = fchname(1:i-1)//'.mol'

 i = INDEX(fchname, '_NO.fch', back=.true.)
 if(i == 0) i = INDEX(fchname, '.fch')
 dalname = fchname(1:i-1)//'_NMR.dal'
 molname = fchname(1:i-1)//'_NMR.mol'
 call fch2dal_wrap(fchname)
 i = RENAME(TRIM(oldmol), TRIM(molname))

 nclosed = ndb + npair - npair0
 open(newunit=fid,file=TRIM(dalname),status='replace')
 write(fid,'(A)') '**DALTON INPUT'
 if(DKH2) write(fid,'(A)') '.DOUGLAS-KROLL'
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
  if(iroot > 0) write(fid,'(A,/,I0)') '.STATE',iroot+1
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
 character(len=240) :: proname, molname, molname1, dalname, dalname1
 character(len=240), intent(in) :: inpname

 do i = 1, 3   ! x,y,z
  bound(1,i) = MINVAL(coor(i,:)) - r_ext
  bound(2,i) = MAXVAL(coor(i,:)) + r_ext
  ngrid(i) = INT((bound(2,i) - bound(1,i))/interval)
 end do ! for i
 origin = bound(1,:)
 tot_ngrid = PRODUCT(ngrid)
 write(6,'(/,A,I0)') 'ICSS rectangular box range: bound atoms + 6.5 Angstrom.&
                    & Total grid = ', tot_ngrid
 write(6,'(A,F6.3,A,3I5)') 'Interval=', interval, '. Ngrid along x,y,z are', &
                           (ngrid(i),i=1,3)

 i = INDEX(inpname, '.dal', back=.true.)
 j = INDEX(inpname, '_NMR.dal', back=.true.)
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

 i = INDEX(buf, ' ')
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

subroutine submit_dalton_icss_job(proname, mem, nproc, nfile)
 implicit none
 integer :: i, n, njob, fid, system
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

 do i = 32, 2, -1
  if(MOD(nproc,i) == 0) exit
 end do ! for i
 njob = i
 write(fid,'(A,I0)') 'mem = ', MIN(16,CEILING(DBLE(mem)/DBLE(i)))
 n = nproc/i
 write(fid,'(A,I0)') 'np = ', n

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
  write(fid,'(A1,2(A,I5.5),A)') ACHAR(9),"dalton -gb ${mem} -omp ${np} -noarch&
                              & -ow ${pro}",i,' >${pro}',i,".sout 2>&1"
 end do ! for i
 close(fid)

 write(buf,'(A,I0,A)') 'make -j',njob,' -f '//TRIM(makefile)//' all'
 write(6,'(A)') '$'//TRIM(buf)
 i = SYSTEM(TRIM(buf))
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
 character(len=240) :: cubname, bqfile
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
  if(INDEX(buf,'Intera')>0 .or. INDEX(buf,'intera')>0) exit
  write(6,'(A)') TRIM(buf)
 end do ! for while

 close(fid)
end subroutine read_shieldings_from_dalton_out

! print (Open)Molcas (DMRG-)CASCI/CASSCF input file keywords
subroutine prt_molcas_cas_para(fid, dmrg, nevpt, chemps2, CIonly, inpname)
 use mr_keyword, only: hardwfn, crazywfn, MaxM, iroot, xmult, nstate
 use mol, only: nacte, nacto, charge, mult
 implicit none
 integer :: i, nroots, lroots
 integer, intent(in) :: fid
 integer, allocatable :: weight(:)
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: dmrg, nevpt, chemps2, CIonly

 if(iroot>0 .and. mult/=xmult) then
  write(6,'(/,A)') 'ERROR in subroutine prt_molcas_cas_para: Xmult/=mult!'
  write(6,'(A)') 'When calling OpenMolcas, the ground state spin must be equal &
                 &to the excited'
  write(6,'(A)') 'state spin.'
  stop
 end if

 if(dmrg) then
  write(fid,'(/,A)') "&DMRGSCF"
  write(fid,'(A)') 'ActiveSpaceOptimizer= QCMaquis'
  write(fid,'(A)') 'Fiedler= ON'
  write(fid,'(A)') 'DMRGSettings'
  write(fid,'(A)') ' nsweeps= 5'
  write(fid,'(A,I0)') ' max_bond_dimension= ', MaxM
  write(fid,'(A)') 'EndDMRGSettings'
  write(fid,'(A)') 'OOptimizationSettings'
  if(chemps2) then
   write(fid,'(A,I0)') 'DMRG= ', maxM
   write(fid,'(A)') '3RDM'
  end if
 else
  write(fid,'(/,A)') "&RASSCF"
  write(fid,'(A)') 'CIMX= 200'
  write(fid,'(A)') 'Tight= 5d-8 5d-6'
  if(crazywfn) then
   write(fid,'(A)') 'SDav= 500'
  else if(hardwfn) then
   write(fid,'(A)') 'SDav= 300'
  end if
 end if

 write(fid,'(A,I0)') 'Spin= ', mult
 write(fid,'(A,I0)') 'Charge= ', charge
 write(fid,'(A,I0)') 'nActEl= ', nacte
 write(fid,'(A,I0)') 'RAS2= ', nacto
 if(iroot > 9) then
  close(fid)
  write(6,'(/,A)') 'ERROR in subroutine prt_molcas_cas_para: iroot>9. Please c&
                   &ontact the MOKIT'
  write(6,'(A)') 'developers to modify the input file print format.'
  stop
 end if

 ! SS-CASCI
 if(iroot > 0) then
  nroots = iroot + 1
  lroots = iroot + 3
  if(nacte==2 .and. nacto==2) lroots = 3 ! CAS(2,2)
  allocate(weight(nroots), source=0)
  weight(nroots) = 1
  write(fid,'(A,2(1X,I0),A)',advance='no') 'CIroot=',nroots,lroots,';'
  write(fid,'(9I2)',advance='no') (i,i=1,nroots)
  write(fid,'(A,9I2)') ';', weight
  deallocate(weight)
 end if

 ! multi-root CASCI
 if(nstate > 0) then
  nroots = nstate + 1
  lroots = nstate + 3
  if(nacte==2 .and. nacto==2) lroots = 3 ! CAS(2,2)
  write(fid,'(A,2(1X,I0),A)') 'CIroot=',nroots,lroots,' 1'
 end if

 i = INDEX(inpname, '.input', back=.true.)
 write(fid,'(A)') 'FILEORB= '//inpname(1:i-1)//'.INPORB'

 if(CIonly) write(fid,'(A)') 'CIonly'
 if(nevpt) write(fid,'(A,/,A)') 'NEVPT2Prep','EvRDM'
 if(dmrg) write(fid,'(A)') 'EndOOptimizationSettings'

 if(nstate > 0) then   ! generate hole and particle NTOs
  write(fid,'(/)',advance='no')
  ! these CASCI states share the same basis set and the same MOs, so copy from
  ! one common $Project.JobIph
  do i = 1, nstate+1, 1
   write(fid,'(A,I3.3)') '>> COPY $Project.JobIph JOB',i
  end do ! for i
  write(fid,'(/,A)') "&RASSI"
  write(fid,'(A)') 'NTOC'
  write(fid,'(A,I0,99I2)') 'Nr of JobIphs=',nstate+1,(1,i=1,nstate+1)
  do i = 1, nstate, 1
   write(fid,'(I0,A1)',advance='no') i,';'
  end do ! for i
  write(fid,'(I0)') nstate+1
 end if
end subroutine prt_molcas_cas_para

! print PySCF FCI solver keywords
subroutine prt_hard_or_crazy_casci_pyscf(fid, nopen, hardwfn, crazywfn, skip2)
 implicit none
 integer, intent(in) :: fid, nopen
 real(kind=8) :: ss
 logical, intent(in) :: hardwfn, crazywfn, skip2

 if(skip2) then
  if(hardwfn) then
   write(fid,'(2X,A)') 'mc.fcisolver.pspace_size = 900'
   write(fid,'(2X,A)') 'mc.fcisolver.max_cycle = 400'
  else if(crazywfn) then
   write(fid,'(2X,A)') 'mc.fcisolver.level_shift = 0.2'
   write(fid,'(2X,A)') 'mc.fcisolver.pspace_size = 1400'
   write(fid,'(2X,A)') 'mc.fcisolver.max_space = 100'
   write(fid,'(2X,A)') 'mc.fcisolver.max_cycle = 600'
   ss = DBLE(nopen)*0.5d0
   ss = ss*(ss+1d0)
   write(fid,'(2X,A,F7.3,A)') 'mc.fix_spin_(ss=',ss,')'
  else
   write(fid,'(2X,A)') 'mc.fcisolver.max_cycle = 200'
  end if
 else
  if(hardwfn) then
   write(fid,'(A)') 'mc.fcisolver.pspace_size = 900'
   write(fid,'(A)') 'mc.fcisolver.max_cycle = 400'
  else if(crazywfn) then
   write(fid,'(A)') 'mc.fcisolver.level_shift = 0.2'
   write(fid,'(A)') 'mc.fcisolver.pspace_size = 1400'
   write(fid,'(A)') 'mc.fcisolver.max_space = 100'
   write(fid,'(A)') 'mc.fcisolver.max_cycle = 600'
   ss = DBLE(nopen)*0.5d0
   ss = ss*(ss+1d0)
   write(fid,'(A,F7.3,A)') 'mc.fix_spin_(ss=',ss,')'
  else
   write(fid,'(A)') 'mc.fcisolver.max_cycle = 200'
  end if
 end if
end subroutine prt_hard_or_crazy_casci_pyscf

! print PySCF FCI solver keywords
subroutine prt_hard_or_crazy_casci_orca(fid, hardwfn, crazywfn)
 implicit none
 integer, intent(in) :: fid
 logical, intent(in) :: hardwfn, crazywfn

 write(fid,'(A)') ' CI'
 if(hardwfn) then
  write(fid,'(A)') '  MaxIter 400'
  write(fid,'(A)') '  NGuessMat 1000'
 else if(crazywfn) then
  write(fid,'(A)') '  MaxIter 600'
  write(fid,'(A)') '  NGuessMat 1500'
 else
  write(fid,'(A)') '  MaxIter 200'
 end if

 write(fid,'(A)') ' end'
end subroutine prt_hard_or_crazy_casci_orca

! print ground state CASSCF keywords into a PySCF .py file
subroutine prt_gs_casscf_kywrd_py(fid, RIJK_bas1)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, casscf, dkh2_or_x2c, RI, maxM, hardwfn, &
  crazywfn, block_mpi
 implicit none
 integer :: i
 integer, intent(in) :: fid
 character(len=21), intent(in) :: RIJK_bas1

 if(casscf) then ! CASSCF
  if(dkh2_or_x2c) then
   write(fid,'(A)',advance='no') 'mc = mcscf.CASSCF(mf.x2c1e(),'
  else
   write(fid,'(A)',advance='no') 'mc = mcscf.CASSCF(mf,'
  end if
  write(fid,'(3(I0,A))',advance='no') nacto,',(',nacta,',',nactb,')'
  if(RI) then
   write(fid,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
  else
   write(fid,'(A)') ')'
  end if
  write(fid,'(A,I0,A)') 'mc.max_memory = ', mem*700, ' # MB'
  write(fid,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*300, ' # MB'
  write(fid,'(A,I0)') 'mc.fcisolver.spin = ', nacta-nactb
  call prt_hard_or_crazy_casci_pyscf(fid, nacta-nactb,hardwfn,crazywfn,.false.)
  write(fid,'(A)') 'mc.natorb = True'
 else ! DMRG-CASSCF
  if(dkh2_or_x2c) then
   write(fid,'(A)',advance='no') 'mc = dmrgscf.DMRGSCF(mf.x2c1e(),'
  else
   write(fid,'(A)',advance='no') 'mc = dmrgscf.DMRGSCF(mf,'
  end if
  write(fid,'(3(I0,A))') nacto,',(',nacta,',',nactb,'))'
  if(block_mpi) then
   i = CEILING(0.5*REAL(mem)/REAL(nproc))
   write(fid,'(A,I0,A)') 'mc.max_memory = ', (mem-nproc*i)*1000, ' # MB'
  else
   i = CEILING(0.4*REAL(mem))
   write(fid,'(A,I0,A)') 'mc.max_memory = ', (mem-i)*1000, ' # MB'
   write(fid,'(A)') 'mc.fcisolver.threads = nproc'
  end if
  write(fid,'(A,I0,A)') 'mc.fcisolver.memory = ', i,' # GB'
  write(fid,'(A,I0)') 'mc.fcisolver.maxM = ', maxM
 end if

 write(fid,'(A)') 'mc.max_cycle = 200'
end subroutine prt_gs_casscf_kywrd_py

! print excited state CASSCF keywords into a PySCF .py file
! state tracking is achieved by detecting spin multiplicities of each state
subroutine prt_es_casscf_kywrd_py(fid, RIJK_bas1)
 use mol, only: nacto, nacta, nactb, mult
 use mr_keyword, only: mem, casscf, dkh2_or_x2c, RI, hardwfn, crazywfn, iroot,&
  xmult
 implicit none
 integer, intent(in) :: fid
 integer, parameter :: maxcyc = 250
 real(kind=8) :: spin, xss ! xss: spin square of the target excited state
 character(len=21), intent(in) :: RIJK_bas1

 spin = DBLE(xmult-1)/2
 xss = spin*(spin + 1d0)

 if(casscf) then ! CASSCF
  if(xmult == mult) then
   write(fid,'(A,I0)') 'nroots = ', iroot+6 ! initial nroots
  else
   write(fid,'(A,I0)') 'nroots = ', iroot+3 ! initial nroots
  end if
  write(fid,'(A,I0,A)') 'for i in range(',maxcyc,'):'
  write(fid,'(2X,A)') "print('ITER=',i)"
  if(dkh2_or_x2c) then
   write(fid,'(2X,A)',advance='no') 'mc = mcscf.CASCI(mf.x2c1e(),'
  else
   write(fid,'(2X,A)',advance='no') 'mc = mcscf.CASCI(mf,'
  end if
  write(fid,'(3(I0,A))',advance='no') nacto,',(',nacta,',',nactb,')'
  if(RI) then
   write(fid,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
  else
   write(fid,'(A)') ')'
  end if
  write(fid,'(2X,A,I0,A)') 'mc.max_memory = ', mem*700, ' # MB'
  write(fid,'(2X,A,I0,A)') 'mc.fcisolver.max_memory = ', mem*300, ' # MB'
  write(fid,'(2X,A)') 'mc.fcisolver.nroots = nroots'
  write(fid,'(2X,A,I0)') 'mc.fcisolver.spin = ', nacta-nactb
  call prt_hard_or_crazy_casci_pyscf(fid, nacta-nactb,hardwfn,crazywfn,.true.)
  write(fid,'(2X,A)') 'mc.verbose = 5'
  write(fid,'(2X,A)') 'mc.kernel()'
  if(xmult == mult) then
   write(fid,'(2X,A)') 'iroot = -1'
  else
   write(fid,'(2X,A)') 'iroot = 0'
  end if
  write(fid,'(2X,A)') 'for j in range(nroots):'
  write(fid,'(4X,A)') 'ss = mc.fcisolver.spin_square(mc.ci[j], mc.ncas, mc.nelecas)'
  write(fid,'(4X,A,F7.3,A)') 'if(abs(ss[0]-',xss,') < 1e-4):'
  write(fid,'(6X,A)') 'iroot = iroot + 1'
  write(fid,'(4X,A,I0,A)') 'if(iroot == ',iroot,'):'
  write(fid,'(6X,A)') 'break'
  write(fid,'(2X,A)') 'e = mc.e_tot[j]'
  write(fid,'(2X,A)') 'nroots = j + 3'
  if(dkh2_or_x2c) then
   write(fid,'(2X,A)',advance='no') 'mc = mcscf.CASSCF(mf.x2c1e(),'
  else
   write(fid,'(2X,A)',advance='no') 'mc = mcscf.CASSCF(mf,'
  end if
  write(fid,'(3(I0,A))',advance='no') nacto,',(',nacta,',',nactb,')'
  write(fid,'(A,I0)',advance='no') ').state_specific_(j'
  if(RI) then
   write(fid,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
  else
   write(fid,'(A)') ')'
  end if
  write(fid,'(2X,A,I0,A)') 'mc.max_memory = ', mem*700, ' # MB'
  write(fid,'(2X,A,I0,A)') 'mc.fcisolver.max_memory = ', mem*300, ' # MB'
  write(fid,'(2X,A)') 'mc.fcisolver.nroots = nroots'
  write(fid,'(2X,A,I0)') 'mc.fcisolver.spin = ', nacta-nactb
  call prt_hard_or_crazy_casci_pyscf(fid, nacta-nactb, hardwfn,crazywfn,.true.)
  write(fid,'(2X,A)') 'mc.verbose = 5'
  write(fid,'(2X,A)') 'mc.max_cycle = 1' ! only 1 cycle
  write(fid,'(2X,A)') 'mc.kernel()'
  write(fid,'(2X,A)') 'mf.mo_coeff = mc.mo_coeff.copy()'
  write(fid,'(2X,A)') 'if(abs(e - mc.e_tot) < 1e-7):'
  write(fid,'(4X,A)') 'break'
  ! assume the SS-CASSCF orbital optimization is accomplished, now run a multi
  !  -root CASCI calculation and generate NOs
  ! print a label to show that this is a successful SS-CASSCF job
  write(fid,'(A)') "print('SSS')"
  ! print the target state and nroots for possible use of NEVPT2, etc
  write(fid,'(A)') "f = open('ss-cas.txt', 'w+')"
  write(fid,'(A)') "f.write('nroots=%i\n' %nroots)"
  write(fid,'(A)') "f.write('target_root=%i' %j)"
  write(fid,'(A)') 'f.close()'
  if(dkh2_or_x2c) then
   write(fid,'(A)',advance='no') 'mc = mcscf.CASCI(mf.x2c1e(),'
  else
   write(fid,'(A)',advance='no') 'mc = mcscf.CASCI(mf,'
  end if
  write(fid,'(3(I0,A))',advance='no') nacto,',(',nacta,',',nactb,')'
  write(fid,'(A,I0)',advance='no') ').state_specific_(j'
  if(RI) then
   write(fid,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
  else
   write(fid,'(A)') ')'
  end if
  write(fid,'(A,I0,A)') 'mc.max_memory = ', mem*700, ' # MB'
  write(fid,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*300, ' # MB'
  write(fid,'(A)') 'mc.fcisolver.nroots = nroots'
  write(fid,'(A,I0)') 'mc.fcisolver.spin = ', nacta-nactb
  call prt_hard_or_crazy_casci_pyscf(fid, nacta-nactb,hardwfn,crazywfn,.false.)
  write(fid,'(A)') 'mc.natorb = True'
 else ! DMRG-CASSCF
  write(6,'(/,A)') 'ERROR in subroutine prt_es_casscf_kywrd_py: the excited sta&
                   &te SS-DMRG-CASSCF'
  write(6,'(A)') 'job is not supported by MOKIT currently.'
  stop
 end if
end subroutine prt_es_casscf_kywrd_py

subroutine prt_casci_kywrd_py(fid, RIJK_bas1, natorb)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: dmrgscf, dkh2_or_x2c, iroot, RI, mem, casci, hardwfn, &
  crazywfn, maxM, nproc, block_mpi
 implicit none
 integer :: i
 integer, intent(in) :: fid
 character(len=21), intent(in) :: RIJK_bas1
 logical, intent(in) :: natorb

 if(natorb) then
  write(fid,'(/,A)') '# perform DMRG-CASCI to get NOs'
  write(fid,'(A)') 'mf.mo_coeff = mc.mo_coeff.copy()'
 end if

 if((.not.natorb) .or. (natorb .and. dmrgscf)) then
  if(dkh2_or_x2c) then
   write(fid,'(A)',advance='no') 'mc = mcscf.CASCI(mf.x2c1e(),'
  else
   write(fid,'(A)',advance='no') 'mc = mcscf.CASCI(mf,'
  end if
  write(fid,'(3(I0,A))',advance='no') nacto,',(',nacta,',',nactb,')'
  if(iroot > 0) write(fid,'(A,I0)',advance='no') ').state_specific_(',iroot
  if(RI) then
   write(fid,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
  else
   write(fid,'(A)') ')'
  end if
  if(casci) then ! CASCI
   write(fid,'(A,I0,A)') 'mc.max_memory = ', mem*700, ' # MB'
   write(fid,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*300, ' # MB'
   write(fid,'(A,I0)') 'mc.fcisolver.spin = ', nacta-nactb
   call prt_hard_or_crazy_casci_pyscf(fid,nacta-nactb,hardwfn,crazywfn,.false.)
   write(fid,'(A)') 'mc.natorb = True'
  else           ! DMRG-CASCI
   write(fid,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=', maxM, ')'
   if(block_mpi) then
    i = CEILING(0.5*REAL(mem)/REAL(nproc))
    write(fid,'(A,I0,A)') 'mc.max_memory = ', (mem-nproc*i)*1000, ' # MB'
   else
    i = CEILING(0.4*REAL(mem))
    write(fid,'(A,I0,A)') 'mc.max_memory = ', (mem-i)*1000, ' # MB'
    write(fid,'(A)') 'mc.fcisolver.threads = nproc'
   end if
   write(fid,'(A,I0,A)') 'mc.fcisolver.memory = ', i, ' # GB'
  end if
 end if

 if(natorb) then
  write(fid,'(A,/,A,/,A)') 'mc.natorb = True','mc.verbose = 4','mc.kernel()'
 end if
end subroutine prt_casci_kywrd_py

! print warnings if user-requested active orbitals/electrons are more than
! recommended.
subroutine prt_active_space_warn(nacte_wish, nacto_wish, nacte, nacto)
 implicit none
 integer, intent(in) :: nacte_wish, nacto_wish, nacte, nacto
 logical :: alive1, alive2

 alive1 = (nacte_wish > nacte)
 alive2 = (nacto_wish > nacto)
 if(alive1 .or. alive2) write(6,'(A)') REPEAT('-',79)

 if(alive1) then
  write(6,'(A)') 'Warning from subroutine do_cas: You request more active ele&
                 &ctrons than recommended.'
  write(6,'(A)') 'You should clearly know what you are calculating, otherwise&
                 & nonsense results may'
  write(6,'(A)') 'be obtained.'
 end if

 if(alive2) then
  write(6,'(A)') 'Warning from subroutine do_cas: You request more active orb&
                 &itals than recommended.'
  write(6,'(A)') 'You should clearly know what you are calculating, otherwise&
                 & nonsense results may'
  write(6,'(A)') 'be obtained.'
 end if

 if(alive1 .or. alive2) write(6,'(A)') REPEAT('-',79)
end subroutine prt_active_space_warn

