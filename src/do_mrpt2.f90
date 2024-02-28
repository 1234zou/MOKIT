! written by jxzou at 20210111: moved from subroutine do_mrpt2 in automr.f90
! updated by jxzou at 20210224: add CASPT2 interface with ORCA

! do CASCI/CASPT2(npair<=7) or DMRG-CASCI/DMRG-NEVPT2 (npair>7) when scf=.False.
! do CASSCF/CASPT2(npair<=7) or DMRG-CASSCF/DMRG-NEVPT2 (npair>7) when scf=.True.
subroutine do_mrpt2()
 use mr_keyword, only: casci, casscf, dmrgci, dmrgscf, dmrg_no, CIonly, caspt2,&
  caspt2k, nevpt2, mrmp2, ovbmp2, sdspt2, casnofch, casscf_prog, casci_prog, &
  nevpt2_prog, caspt2_prog, bgchg, chgname, mem, nproc, gms_path, gms_scr_path,&
  check_gms_path, molcas_omp, molcas_path, molpro_path, orca_path, bdf_path, &
  gau_path, FIC, eist, target_root, caspt2_force
 use mol, only: caspt2_e, nevpt2_e, mrmp2_e, sdspt2_e, ovbmp2_e, davidson_e, &
  ptchg_e, nuc_pt_e, natom, grad
 use util_wrapper, only: bas_fch2py_wrap, mkl2gbw, fch2inp_wrap, unfchk, &
  fch2inporb_wrap
 implicit none
 integer :: i, RENAME, system
 character(len=24) :: data_string
 character(len=240) :: string, pyname, outname, inpname, inporb
 character(len=240) :: mklname, cmofch
 real(kind=8) :: ref_e, corr_e
 logical :: alive(5)

 if(eist == 1) return ! excited state calculation
 alive = [caspt2, nevpt2, mrmp2, ovbmp2, sdspt2]
 if(ALL(alive .eqv. .false.)) return
 write(6,'(//,A)') 'Enter subroutine do_mrpt2...'

 if((dmrgci .or. dmrgscf)) then
  if(mrmp2) then
   write(6,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-MRMP2 not supported.'
   stop
  end if
  if(sdspt2) then
   write(6,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-SDSPT2 not supported.'
   stop
  end if
  if(ovbmp2) then
   write(6,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-OVB-MP2 not supported.'
   stop
  end if
  if(nevpt2) then
   select case(TRIM(nevpt2_prog))
   case('molpro', 'orca')
    write(6,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-NEVPT2 is only support&
                     &ed with PySCF or OpenMolcas.'
    write(6,'(A)') 'But you specify NEVPT2_prog='//TRIM(nevpt2_prog)
    stop
   end select
  end if
  if(caspt2 .and. TRIM(caspt2_prog)/='openmolcas') then
   write(6,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-CASPT2 is only supporte&
                    &d with OpenMolcas.'
   write(6,'(A)') 'But you specify CASPT2_prog='//TRIM(caspt2_prog)
   stop
  end if
 end if

 if(.not. CIonly) then
  if(TRIM(casscf_prog) == 'orca') then
   write(6,'(A)') REPEAT('-',79)
   write(6,'(A)') 'Warning: ORCA is used as the CASSCF solver, the NO coefficie&
                  &nts in .mkl file are'
   write(6,'(A)') 'only 7 digits. This will affect the PT2 energy up to 10^-5 a&
                  &.u. Such small error'
   write(6,'(A)') 'is usually not important. If you care about the accuracy, pl&
                  &ease use another'
   write(6,'(A)') 'CASSCF solver.'
   write(6,'(A)') REPEAT('-',79)
  end if

  if(casscf) then
   if(caspt2) then
    string = 'CASPT2 based on CASSCF orbitals.'
   else if(nevpt2) then
    string = 'NEVPT2 based on CASSCF orbitals.'
   else if(mrmp2) then
    string = 'MRMP2 based on CASSCF orbitals.'
   else if(ovbmp2) then
    string = 'OVB-MP2 based on CASSCF orbitals.'
   else
    string = 'SDSMP2 based on CASSCF orbitals.'
   end if
  else ! DMRG-CASSCF
   if(caspt2) then
    string = 'DMRG-CASPT2 based on DMRG-CASSCF orbitals.'
   else if(nevpt2) then
    string = 'DMRG-NEVPT2 based on DMRG-CASSCF orbitals.'
   end if
  end if

 else ! CIonly = .True.
  if(TRIM(casci_prog) == 'orca') then
   write(6,'(A)') REPEAT('-',79)
   write(6,'(A)') 'Warning: ORCA is used as the CASCI solver, the NO coefficien&
                  &ts in .mkl file are'
   write(6,'(A)') 'only 7-digits. This will affect the PT2 energy up to 10^-5 a&
                  &.u. Such small error'
   write(6,'(A)') 'is usually not important. If you care about the accuracy, pl&
                  &ease use another CASCI solver.'
   write(6,'(A)') REPEAT('-',79)
  end if

  if(casci) then
   if(caspt2) then
    string = 'CASPT2 based on CASCI orbitals.'
   else if(nevpt2) then
    string = 'NEVPT2 based on CASCI orbitals.'
   else if(ovbmp2) then
    string = 'OVB-MP2 based on CASCI orbitals.'
   else
    string = 'SDSPT2 based on CASCI orbitals.'
   end if
  else ! DMRG-CASCI
   if(caspt2) then
    string = 'DMRG-CASPT2 based on DMRG-CASCI orbitals.'
   else
    string = 'DMRG-NEVPT2 based on DMRG-CASCI orbitals.'
   end if
  end if
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning: the CASSCF orbital optimization is strongly recommen&
                 &ded to be'
  write(6,'(A)') 'performed before PT2, unless it is too time-consuming.'
  write(6,'(A)') REPEAT('-',79)
 end if

 write(6,'(A)') TRIM(string)
 write(6,'(A)',advance='no') 'Frozen_core = F, '

 if(nevpt2) then
  write(6,'(A)') 'NEVPT2 using program '//TRIM(nevpt2_prog)

  select case(TRIM(nevpt2_prog))
  case('pyscf')
   ! For DMRG-NEVPT2, use CMOs rather than NOs
   if(dmrgci .or. dmrgscf) then ! DMRG-CASCI/CASSCF
    if(dmrg_no) then ! DMRG NOs calculated previously
     i = INDEX(casnofch, '_NO', back=.true.)
     cmofch = casnofch(1:i)//'CMO.fch'
     casnofch = cmofch
    end if
    ! if DMRG NOs not calculated previously, casnofch is set to 'xxx_CMO.fch'
    ! in subroutine do_cas
    i = INDEX(casnofch, '_CMO', back=.true.)
   else                         ! CASCI/CASSCF
    i = INDEX(casnofch, '_NO', back=.true.)
   end if
   pyname = casnofch(1:i)//'NEVPT2.py'
   outname = casnofch(1:i)//'NEVPT2.out'
   call bas_fch2py_wrap(casnofch, .false., pyname)
   call prt_nevpt2_script_into_py(pyname)
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(pyname))
   call submit_pyscf_job(pyname)

  case('molpro')
   call check_exe_exist(molpro_path)

   i = SYSTEM('fch2com '//TRIM(casnofch))
   i = INDEX(casnofch, '.fch', back=.true.)
   pyname = casnofch(1:i-1)//'.com'
   string = casnofch
   call convert2molpro_fname(string, '.a')
   i = INDEX(casnofch, '_NO', back=.true.)
   inpname = casnofch(1:i-1)//'_NEVPT2.com'
   outname = casnofch(1:i-1)//'_NEVPT2.out'
   inporb = inpname
   call convert2molpro_fname(inporb, '.a')
   i = RENAME(TRIM(pyname), TRIM(inpname))
   i = RENAME(TRIM(string), TRIM(inporb))

   call prt_mrpt_molpro_inp(inpname, 1)
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   call submit_molpro_job(inpname, mem, nproc)

  case('openmolcas')
   call check_exe_exist(molcas_path)
   i = INDEX(casnofch, '_NO', back=.true.)
   inpname = casnofch(1:i)//'NEVPT2.input'
   outname = casnofch(1:i)//'NEVPT2.out'
   call fch2inporb_wrap(casnofch, .true., inpname)
   call prt_nevpt2_molcas_inp(inpname)
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   call submit_molcas_job(inpname, mem, nproc, molcas_omp)
   i = 0 ! the value of i will be checked below

  case('orca')
   call check_exe_exist(orca_path)
   i = SYSTEM('fch2mkl '//TRIM(casnofch))
   i = INDEX(casnofch, '.fch', back=.true.)
   inporb = casnofch(1:i-1)//'_o.mkl'
   string = casnofch(1:i-1)//'_o.inp'
   i = INDEX(casnofch, '_NO', back=.true.)
   mklname = casnofch(1:i)//'NEVPT2.mkl'
   inpname = casnofch(1:i)//'NEVPT2.inp'
   outname = casnofch(1:i)//'NEVPT2.out'
   i = RENAME(TRIM(inporb), TRIM(mklname))
   i = RENAME(TRIM(string), TRIM(inpname))

   call prt_nevpt2_orca_inp(inpname)
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   ! if bgchg = .True., .inp and .mkl file will be updated
   call mkl2gbw(mklname)
   call delete_file(mklname)
   call submit_orca_job(orca_path, inpname, .true.)

  case('bdf')
   call check_exe_exist(bdf_path)
   i = SYSTEM('fch2bdf '//TRIM(casnofch)//' -no')
   i = INDEX(casnofch, '.fch', back=.true.)
   inporb = casnofch(1:i-1)//'_bdf.inporb'
   string = casnofch(1:i-1)//'_bdf.inp'
   i = INDEX(casnofch, '_NO', back=.true.)
   mklname = casnofch(1:i)//'NEVPT2.inporb'
   inpname = casnofch(1:i)//'NEVPT2.inp'
   outname = casnofch(1:i)//'NEVPT2.out'
   i = RENAME(TRIM(inporb), TRIM(mklname))
   i = RENAME(TRIM(string), TRIM(inpname))

   call prt_mrpt_bdf_inp(inpname, 2)
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   i = INDEX(inpname, '.inp', back=.true.)
   string = inpname(1:i-1)
   i = SYSTEM(TRIM(bdf_path)//' '//TRIM(string))
  end select

 else if(caspt2) then ! CASPT2
  write(6,'(A)') 'CASPT2 using program '//TRIM(caspt2_prog)

  select case(TRIM(caspt2_prog))
  case('openmolcas')
   call check_exe_exist(molcas_path)
   ! For DMRG-CASPT2, use CMOs rather than NOs
   i = INDEX(casnofch, '_NO', back=.true.)
   if(dmrgci .or. dmrgscf) then
    cmofch = casnofch(1:i)//'CMO.fch'
    casnofch = cmofch
   end if
   inpname = casnofch(1:i-1)//'_CASPT2.input'
   outname = casnofch(1:i-1)//'_CASPT2.out'

   call fch2inporb_wrap(casnofch, .true., inpname)
   call prt_caspt2_molcas_inp(inpname)
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   call submit_molcas_job(inpname, mem, nproc, molcas_omp)
   i = 0 ! the value of i will be checked below

  case('molpro')
   call check_exe_exist(molpro_path)

   i = SYSTEM('fch2com '//TRIM(casnofch))
   i = INDEX(casnofch, '.fch', back=.true.)
   pyname = casnofch(1:i-1)//'.com'
   string = casnofch
   call convert2molpro_fname(string, '.a')
   i = INDEX(casnofch, '_NO', back=.true.)
   inpname = casnofch(1:i-1)//'_CASPT2.com'
   outname = casnofch(1:i-1)//'_CASPT2.out'
   inporb = inpname
   call convert2molpro_fname(inporb, '.a')
   i = RENAME(TRIM(pyname), TRIM(inpname))
   i = RENAME(TRIM(string), TRIM(inporb))

   call prt_mrpt_molpro_inp(inpname, 2)
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   call submit_molpro_job(inpname, mem, nproc)

  case('orca')
   call check_exe_exist(orca_path)

   i = SYSTEM('fch2mkl '//TRIM(casnofch))
   i = INDEX(casnofch, '.fch', back=.true.)
   inporb = casnofch(1:i-1)//'_o.mkl'
   string = casnofch(1:i-1)//'_o.inp'
   i = INDEX(casnofch, '_NO', back=.true.)
   mklname = casnofch(1:i)//'CASPT2.mkl'
   inpname = casnofch(1:i)//'CASPT2.inp'
   outname = casnofch(1:i)//'CASPT2.out'
   i = RENAME(TRIM(inporb), TRIM(mklname))
   i = RENAME(TRIM(string), TRIM(inpname))

   call prt_caspt2_orca_inp(inpname)
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   ! if bgchg = .True., .inp and .mkl file will be updated
   call mkl2gbw(mklname)
   call delete_file(mklname)
   call submit_orca_job(orca_path, inpname, .true.)
  end select

 else if(mrmp2) then ! CASSCF-MRMP2
  write(6,'(A)') 'MRMP2 using program gamess'
  call check_gms_path()

  call fch2inp_wrap(casnofch, .false., 0, 0)
  i = INDEX(casnofch, '.fch')
  pyname = casnofch(1:i-1)//'.inp'
  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i-1)//'_MRMP2.inp'

  outname = casnofch(1:i-1)//'_MRMP2.gms'
  i = RENAME(pyname, inpname)
  call prt_mrmp2_gms_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_gms_job(gms_path, gms_scr_path, inpname, nproc)

 else if(ovbmp2) then ! OVB-MP2
  write(6,'(A)') 'OVB-MP2 using program gaussian'
  call check_exe_exist(gau_path)
  i = INDEX(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'OVBMP2.chk'
  inpname = casnofch(1:i)//'OVBMP2.gjf'
  outname = casnofch(1:i)//'OVBMP2.log'
  call unfchk(casnofch, mklname)
  call prt_ovbmp2_gau_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_gau_job(gau_path, inpname, .true.)

 else ! CASSCF-SDSPT2
  write(6,'(A)') 'SDSPT2 using program bdf'
  call check_exe_exist(bdf_path)

  i = SYSTEM('fch2bdf '//TRIM(casnofch)//' -no')
  i = INDEX(casnofch, '.fch', back=.true.)
  inporb = casnofch(1:i-1)//'_bdf.inporb'
  string = casnofch(1:i-1)//'_bdf.inp'
  i = INDEX(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'SDSPT2.inporb'
  inpname = casnofch(1:i)//'SDSPT2.inp'
  outname = casnofch(1:i)//'SDSPT2.out'
  i = RENAME(TRIM(inporb), TRIM(mklname))
  i = RENAME(TRIM(string), TRIM(inpname))

  call prt_mrpt_bdf_inp(inpname, 1)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  i = INDEX(inpname, '.inp', back=.true.)
  string = inpname(1:i-1)
  i = SYSTEM(TRIM(bdf_path)//' '//TRIM(string))
 end if

 if(nevpt2) then      ! read NEVPT2 energy
  select case(TRIM(nevpt2_prog))
  case('pyscf')
   call read_mrpt_energy_from_pyscf_out(outname, target_root, ref_e, corr_e)
   ref_e = ref_e + ptchg_e
  case('molpro')
   if(FIC) then
    call read_mrpt_energy_from_molpro_out(outname, 2, ref_e, corr_e)
   else
    call read_mrpt_energy_from_molpro_out(outname, 1, ref_e, corr_e)
   end if
   ref_e = ref_e + ptchg_e
  case('openmolcas')
   if(FIC) then
    call read_mrpt_energy_from_molcas_out(outname, 2, ref_e, corr_e)
   else
    call read_mrpt_energy_from_molcas_out(outname, 1, ref_e, corr_e)
   end if
   ref_e = ref_e + ptchg_e
  case('orca')
   if(FIC) then
    call read_mrpt_energy_from_orca_out(outname, 2, ref_e, corr_e)
   else
    call read_mrpt_energy_from_orca_out(outname, 1, ref_e, corr_e)
   end if
  case('bdf')
   call read_mrpt_energy_from_bdf_out(outname, 2, ref_e, corr_e, davidson_e)
   ! here the parameter davidson_e is useless
   ref_e = ref_e + ptchg_e + nuc_pt_e
  end select
 else if(caspt2) then ! read CASPT2 energy
  select case(TRIM(caspt2_prog))
  case('openmolcas')
   call read_mrpt_energy_from_molcas_out(outname, 3, ref_e, corr_e)
   ref_e = ref_e + ptchg_e
  case('molpro')
   call read_mrpt_energy_from_molpro_out(outname, 3, ref_e, corr_e)
   ref_e = ref_e + ptchg_e
  case('orca')
   call read_mrpt_energy_from_orca_out(outname, 3, ref_e, corr_e)
  end select
 else if(mrmp2) then  ! read MRMP2 energy
  call read_mrpt_energy_from_gms_out(outname, ref_e, corr_e)
 else if(ovbmp2) then ! read OVB-MP2 energy
  call read_mrpt_energy_from_gau_out(outname, ref_e, corr_e)
 else                 ! read SDSPT2 energy
  call read_mrpt_energy_from_bdf_out(outname, 1, ref_e, corr_e, davidson_e)
 end if

 write(6,'(/,A,F18.8,1X,A4)') 'E(ref)       = ', ref_e, 'a.u.'

 if(caspt2) then
  if(.not. caspt2k) then
   if(caspt2_force) then
    write(6,'(A)') 'IP-EA shift  =         0.0 (for gradient)'
   else
    write(6,'(A)') 'IP-EA shift  =         0.25 (default)'
   end if
  end if
  caspt2_e = ref_e + corr_e
  write(6,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e, 'a.u.'
  if(caspt2k) then
   write(6,'(A,F18.8,1X,A4)') 'E(CASPT2-K)  = ', caspt2_e, 'a.u.'
  else
   write(6,'(A,F18.8,1X,A4)') 'E(CASPT2)    = ', caspt2_e, 'a.u.'
  end if
 else if(nevpt2) then ! NEVPT2
  nevpt2_e = ref_e + corr_e
  write(6,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,   'a.u.'
  if(FIC) then
   write(6,'(A,F18.8,1X,A4)')'E(FIC-NEVPT2)= ', nevpt2_e, 'a.u.'
  else
   write(6,'(A,F18.8,1X,A4)')'E(SC-NEVPT2) = ', nevpt2_e, 'a.u.'
  end if
 else if(mrmp2) then ! MRMP2
  mrmp2_e = ref_e + corr_e
  write(6,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,  'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(MRMP2)     = ', mrmp2_e, 'a.u.'
 else if(ovbmp2) then ! OVB-MP2
  ovbmp2_e = ref_e + corr_e
  write(6,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,  'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(OVB-MP2)   = ', ovbmp2_e,'a.u.'
 else           ! SDSPT2
  sdspt2_e = ref_e + corr_e + davidson_e
  write(6,'(A,F18.8,1X,A4)') 'Davidson correction=',davidson_e,'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,  'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(SDSPT2)    = ', sdspt2_e, 'a.u.'
 end if

 call delete_file('ss-cas.txt')

 if(caspt2_force) then
  allocate(grad(3*natom))

  select case(TRIM(caspt2_prog))
  case('molpro')
   call read_grad_from_molpro_out(outname, 'RSPT2', natom, grad)
  case('openmolcas')
   call read_grad_from_molcas_out(outname, natom, grad)
  case default
   write(6,'(A)') 'ERROR in subroutine do_mrpt2: program cannot be identified.'
   write(6,'(A)') 'CASPT2_prog='//TRIM(caspt2_prog)
   stop
  end select

  write(6,'(/,A)') 'Cartesian gradients (HARTREE/BOHR):'
  write(6,'(5(1X,ES15.8))') (grad(i),i=1,3*natom)
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mrpt2 at '//TRIM(data_string)
end subroutine do_mrpt2

! print NEVPT2/CASPT2/CASPT3 keywords into Molpro input file
subroutine prt_mrpt_molpro_inp(inpname, itype)
 use mr_keyword, only: CIonly, caspt2_force
 use mol, only: ndb, npair, npair0, nacto
 implicit none
 integer :: fid, nclosed, nocc
 integer, intent(in) :: itype ! 1/2/3 for NEVPT2/CASPT2/CASPT3
 character(len=240) :: put, buf, orbfile
 character(len=240), intent(in) :: inpname

 orbfile = inpname
 call convert2molpro_fname(orbfile, '.a')

 put = ' '
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 BACKSPACE(fid)
 read(fid,'(A)') put

 if(put(1:4) /= '{put') then
  write(6,'(/,A)') 'ERROR in subroutine prt_mrpt_molpro_inp: wrong content foun&
                   &d in the final'
  write(6,'(A)') 'line of file '//TRIM(inpname)
  close(fid)
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

 if(CIonly) then
  write(fid,'(2(A,I0),A)') '{CASSCF;closed,',nclosed,';occ,',nocc,&
                           ';DONT,ORBITAL;CANONICAL;NoExtra}'
 else
  write(fid,'(2(A,I0),A)') '{CASSCF;closed,',nclosed,';occ,',nocc,&
                           ';CANONICAL;NoExtra}'
 end if

 write(fid,'(A)') TRIM(put)

 select case(itype)
 case(1) ! NEVPT2
  write(fid,'(A)') '{NEVPT2;CORE}'
 case(2) ! CASPT2
  if(caspt2_force) then
   write(fid,'(A,/,A)') '{RS2,IPEA=0.0;CORE}','Force'
  else
   write(fid,'(A)') '{RS2C,IPEA=0.25;CORE}'
  end if
 case(3) ! CASPT3
  write(fid,'(A)') '{RS3,IPEA=0.25;CORE}'
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine prt_mrpt_molpro_inp: invalid itype=',&
                       itype
  write(6,'(A,I0)') 'Allowed itype are 1/2/3 for NEVPT2/CASPT2/CASPT3.'
  close(fid)
  stop
 end select

 close(fid)
end subroutine prt_mrpt_molpro_inp

! print SDSPT2/NEVPT2/NEVPT3 keywords in to a given BDF .inp file
subroutine prt_mrpt_bdf_inp(inpname, itype)
 use mol, only: charge, mult, nacte, nacto, npair, npair0, ndb
 use mr_keyword, only: CIonly
 implicit none
 integer :: i, nclosed, fid
 integer, intent(in) :: itype ! 1/2/3 for SDSPT2/NEVPT2/NEVPT3
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 nclosed = ndb + npair - npair0
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == '$SCF') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine prt_mrpt_bdf_inp: no '$SCF' found in fi&
                   &le "//TRIM(inpname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 write(fid,'(A)') '$MCSCF'
 write(fid,'(A,/,I0)') 'Charge', charge
 write(fid,'(A,/,I0)') 'Spin', mult
 write(fid,'(A,/,I0)') 'Close', nclosed
 write(fid,'(A,/,I0)') 'Active', nacto
 write(fid,'(A,/,I0)') 'Actel', nacte
 write(fid,'(A)') 'Guess'
 write(fid,'(A)') 'read'
 if(CIonly) write(fid,'(A)') 'CASCI'
 write(fid,'(A)') '$END'

 write(fid,'(/,A)') '$TRAINT'
 write(fid,'(A,/,A1)') 'Frozen', '0'
 write(fid,'(A)') 'Orbital'
 write(fid,'(A)') ' mcorb'
 if(itype == 3) then
  write(fid,'(A)') 'MRPT3'
 else
  write(fid,'(A)') 'MRPT2'
 end if

 write(fid,'(A)') '$END'
 write(fid,'(/,A)') '$XIANCI'

 select case(itype)
 case(1)
  write(fid,'(A)') 'SDSPT2'
 case(2)
  write(fid,'(A)') 'NEVPT2'
 case(3)
  write(fid,'(A)') 'NEVPT3'
 case default
  write(6,'(A,I0)') 'ERROR in subroutine prt_mrpt_bdf_inp: invalid itype=',itype
  write(6,'(A)') 'itype=1/2/3 for SDSPT2/NEVPT2/NEVPT3.'
  stop
 end select

 write(fid,'(A)') '$END'
 close(fid)
end subroutine prt_mrpt_bdf_inp

! print NEVPT2 keywords in to a given ORCA .inp file
subroutine prt_nevpt2_orca_inp(inpname)
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc, X2C, CIonly, RI, RIJK_bas, F12, F12_cabs, &
  FIC, DLPNO, hardwfn, crazywfn
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 read(fid1,'(A)') buf   ! skip nproc
 read(fid1,'(A)') buf   ! skip memory
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid2,'(A,I0,A)') '%maxcore ', CEILING(1d3*DBLE(mem)/DBLE(nproc))

 read(fid1,'(A)') buf   ! skip '!' line
 write(fid2,'(A)',advance='no') '!'
 if(CIonly) write(fid2,'(A)',advance='no') ' noiter'
 if(RI) then
  ! RIJK in CASSCF must be combined with CONVentional
  ! /C basis set cannot be used with /JK in CASSCF
  write(fid2,'(A)',advance='no') ' RIJK conv '//TRIM(RIJK_bas)
 end if
 if(F12) write(fid2,'(A)',advance='no') ' '//TRIM(F12_cabs)
 if(DLPNO) write(fid2,'(A)',advance='no') ' TightPNO'
 write(fid2,'(A)') ' TightSCF'

 if(X2C) then
  write(6,'(/,A)') 'ERROR in subroutine prt_nevpt2_orca_inp: NEVPT2 with X2C is&
                   & not supported in ORCA.'
  write(6,'(A)') 'You can specify NEVPT2_prog=Molpro/OpenMolcas.'
  stop
 end if

 write(fid2,'(A)') '%casscf'
 write(fid2,'(A,I0)') ' nel ', nacte
 write(fid2,'(A,I0)') ' norb ', nacto
 write(fid2,'(A,I0)') ' maxiter 1'
 call prt_hard_or_crazy_casci_orca(fid2, hardwfn, crazywfn)

 if(FIC) then
  if(DLPNO) then
   write(fid2,'(A)') ' PTMethod DLPNO_NEVPT2'
  else
   write(fid2,'(A)') ' PTMethod FIC_NEVPT2'
  end if
 else
  write(fid2,'(A)') ' PTMethod SC_NEVPT2'
 end if
 if(F12) then
  write(fid2,'(A)') ' PTSettings'
  write(fid2,'(A)') '  F12 true'
  write(fid2,'(A)') ' end'
 end if
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
end subroutine prt_nevpt2_orca_inp

! print CASPT2 keywords in to a given ORCA .inp file
subroutine prt_caspt2_orca_inp(inpname)
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc, caspt2k, X2C, CIonly, RI, RIJK_bas, &
  hardwfn, crazywfn
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 read(fid1,'(A)') buf   ! skip nproc
 read(fid1,'(A)') buf   ! skip memory
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid2,'(A,I0,A)') '%maxcore ', CEILING(1d3*DBLE(mem)/DBLE(nproc))

 read(fid1,'(A)') buf   ! skip '!' line
 write(fid2,'(A)',advance='no') '!'
 if(CIonly) write(fid2,'(A)',advance='no') ' noiter'
 if(RI) then
  ! RIJK in CASSCF must be combined with CONVentional
  ! /C basis set cannot be used with /JK in CASSCF
  write(fid2,'(A)',advance='no') ' RIJK conv '//TRIM(RIJK_bas)
 end if
 write(fid2,'(A)') ' TightSCF'

 if(X2C) then
  write(6,'(/,A)') 'ERROR in subroutine prt_nevpt2_orca_inp: CASPT2 with X2C is&
                   & not supported in ORCA.'
  write(6,'(A)') 'You can specify CASPT2_prog=Molpro or OpenMolcas.'
  close(fid2,status='delete')
  close(fid1)
  stop
 end if

 write(fid2,'(A)') '%casscf'
 write(fid2,'(A,I0)') ' nel ', nacte
 write(fid2,'(A,I0)') ' norb ', nacto
 write(fid2,'(A,I0)') ' maxiter 1'
 call prt_hard_or_crazy_casci_orca(fid2, hardwfn, crazywfn)

 if(caspt2k) then
  write(fid2,'(A)') ' PTMethod FIC_CASPT2K'
 else
  write(fid2,'(A)') ' PTMethod FIC_CASPT2'
 end if
 write(fid2,'(A)') ' PTSettings'
 if(caspt2k) then
  write(fid2,'(A)') '  TReg 2e-2'
  ! According to my experience, assuming we are use the same converged CASSCF
  !  orbitals but with different TReg, the default Treg=1e-2 sometimes makes
  !  the CASCI ref weight too small, e.g. for a transition state geometry.
  ! But 2e-2, 5e-2 or 1e-1 give similar and numerically stable results, and
  !  the calculated CASCI ref weight is close to 0.8, so 2e-2 is used here.
 else
  write(fid2,'(A)') '  CASPT2_IPEAshift 0.25'
 end if
 write(fid2,'(A)') '  MaxIter 200'
 write(fid2,'(A)') ' end'
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
end subroutine prt_caspt2_orca_inp

! print NEVPT2 script into a given .py file
subroutine prt_nevpt2_script_into_py(pyname)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, casci, casscf, maxM, X2C, RI, RIJK_bas, &
  hardwfn, crazywfn, iroot, target_root, block_mpi
 implicit none
 integer :: i, nroots, fid1, fid2, RENAME
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, pyname1
 character(len=240), intent(in) :: pyname

 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
 pyname1 = TRIM(pyname)//'.t'

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:17) == 'from pyscf import') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(casci .or. casscf) then
  buf = TRIM(buf)//', mcscf, mrpt, lib'
 else
  buf = TRIM(buf)//', mcscf, dmrgscf, mrpt, lib'
 end if
 write(fid2,'(A)') TRIM(buf)

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(.not. (casci .or. casscf)) then
  write(fid2,'(/,A)',advance='no') "dmrgscf.settings.MPIPREFIX = '"
  if(block_mpi) write(fid2,'(A,I0)',advance='no') 'mpirun -n ', nproc
  write(fid2,'(A)') "'"
 end if
 write(fid2,'(A,I0)') 'nproc = ', nproc
 write(fid2,'(A)') 'lib.num_threads(nproc)'

 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:12) == 'mf.max_cycle') exit
 end do ! for while
 write(fid2,'(A,I0,A)') 'mf.max_memory = ', mem*1000, ' # MB'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')

 write(fid2,'(A)') '# generate CASCI wfn'
 if(X2C) then
  write(fid2,'(A)',advance='no') 'mc = mcscf.CASCI(mf.x2c1e(),'
 else
  write(fid2,'(A)',advance='no') 'mc = mcscf.CASCI(mf,'
 end if
 write(fid2,'(3(I0,A))',advance='no') nacto,',(',nacta,',',nactb,')'
 if(RI) then
  write(fid2,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
 else
  write(fid2,'(A)') ')'
 end if

 if(casci .or. casscf) then
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*700, ' # MB'
  write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*300, ' # MB'
 else
  write(fid2,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=', maxM, ')'
  if(block_mpi) then
   i = CEILING(0.5*REAL(mem)/REAL(nproc))
   write(fid2,'(A,I0,A)') 'mc.max_memory = ', (mem-nproc*i)*1000, ' # MB'
  else
   i = CEILING(0.4*REAL(mem))
   write(fid2,'(A,I0,A)') 'mc.max_memory = ', (mem-i)*1000, ' # MB'
   write(fid2,'(A)') 'mc.fcisolver.threads = nproc'
  end if
  write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ', i, ' # GB'
 end if

 if(iroot > 0) then ! NEVPT2 based on SS-CASSCF
  call read_ss_root_from_txt(nroots, target_root)
  write(fid2,'(A,I0)') 'mc.fcisolver.nroots = ', nroots
 end if

 if(casci .or. casscf) then
  call prt_hard_or_crazy_casci_pyscf(fid2,nacta-nactb,hardwfn,crazywfn,.false.)
 end if
 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.kernel()'

 if(casci .or. casscf) then
  write(fid2,'(/,A)',advance='no') 'mrpt.NEVPT(mc'
  if(iroot > 0) write(fid2,'(A,I0)',advance='no') ', root=', target_root
  write(fid2,'(A)') ').kernel()'
 else
  call prt_dmrg_nevpt2_setting(fid2)
 end if

 close(fid2)
 i = RENAME(pyname1, pyname)
end subroutine prt_nevpt2_script_into_py

! print PySCF DMRG-NEVPT2 settings
subroutine prt_dmrg_nevpt2_setting(fid)
 use mr_keyword, only: mem, nproc, maxM, iroot, target_root
 use mol, only: nif, nacto
 implicit none
 integer :: i, real_nproc, nthread
 integer, intent(in) :: fid

 ! The number of MPI processors should be less than or equal to ndb+nvir. This
 ! is required by dmrgscf/pyscf/dmrgscf/nevpt_mpi.py.
 real_nproc = min(nif-nacto, nproc)

 ! MPI/OpenMP hybrid can be used if the number of MPI processors is small
 nthread = 1
 if(real_nproc < nproc) nthread = nproc/real_nproc

 ! In the NEVPT2 step, always use MPI or MPI/OpenMP hybrid parallelism since
 ! they are faster than pure OpenMP.
 write(fid,'(/,A,I0,A)') "dmrgscf.settings.MPIPREFIX = 'mpirun -n ",real_nproc,&
                         "'"
 write(fid,'(A)') 'lib.num_threads(1)'
 i = CEILING(0.5*REAL(mem)/REAL(real_nproc))
 write(fid,'(A,I0,A)') 'mc.max_memory = ', (mem-real_nproc*i)*1000, ' # MB'
 write(fid,'(A,I0)') 'mc.fcisolver.threads = ', nthread
 write(fid,'(A,I0,A)') 'mc.fcisolver.memory = ', i, ' # GB'
 write(fid,'(A,I0,A)') "mc.fcisolver.mpiprefix = 'mpirun -n ", real_nproc, &
                       " --bind-to none'"
 write(fid,'(A)',advance='no') 'mrpt.NEVPT(mc'
 if(iroot > 0) write(fid,'(A,I0)',advance='no') ', root=', target_root
 write(fid,'(A,I0,A)') ').compress_approx(maxM=', maxM, ').kernel()'
end subroutine prt_dmrg_nevpt2_setting

! print NEVPT2 keywords into OpenMolcas .input file
! It seems that OpenMolcas does not support CASSCF-NEVPT2. So I have to use
! DMRG-NEVPT2.
subroutine prt_nevpt2_molcas_inp(inpname)
 use mr_keyword, only: CIonly, RI, RIJK_bas
 implicit none
 integer :: i, j, fid1, fid2, RENAME
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:7) == 'SEWARD') exit
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
  write(6,'(A)') "ERROR in subroutine prt_nevpt2_molcas_inp: no 'SEWARD'&
                & found in file "//TRIM(inpname)
  stop
 end if

 if(RI) write(fid2,'(A)') 'RIJK'
 write(fid2,'(A)') "&SEWARD"

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4)=="&SCF" .or. buf(1:7)=="&RASSCF") exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine prt_nevpt2_molcas_inp: no 'SCF'&
                & found in file "//TRIM(inpname)
  stop
 end if
 close(fid1,status='delete')

 call prt_molcas_cas_para(fid2, .true., .true., .false., CIonly, inpname)

 write(fid2,'(/,A)') "&MOTRA"
 write(fid2,'(A)') 'Frozen = 0'
 write(fid2,'(A)') 'HDF5'

 write(fid2,'(/,A)') "&NEVPT2"
 write(fid2,'(A)') 'Frozen = 0'
 close(fid2)
 i = RENAME(inpname1, inpname)
end subroutine prt_nevpt2_molcas_inp

! print CASTP2 keywords into OpenMolcas .input file
subroutine prt_caspt2_molcas_inp(inputname)
 use mr_keyword, only: CIonly, dmrgci, dmrgscf, caspt2_force
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inputname1
 character(len=240), intent(in) :: inputname
 logical :: dmrg

 dmrg = (dmrgci .or. dmrgscf)
 inputname1 = TRIM(inputname)//'.t'

 open(newunit=fid1,file=TRIM(inputname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inputname1),status='replace')
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(2:7) == 'SEWARD') exit
 end do

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine prt_caspt2_molcas_inp: no 'SEWARD'&
                & found in file "//TRIM(inputname)
  stop
 end if

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  if(buf(1:4)=="&SCF" .or. buf(1:7)=="&RASSCF") exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid1,status='delete')

 call prt_molcas_cas_para(fid2, dmrg, .false., dmrg, CIonly, inputname)

 if(CIonly) then
  if(dmrgci) then
   write(6,'(A)') 'ERROR in subroutine prt_caspt2_molcas_inp: CIonly is not all&
                  &owed in DMRG-CASPT2.'
   write(6,'(A)') 'It must be based on converged (DMRG-)CASSCF orbitals.'
   stop
  end if
 end if

 if(caspt2_force) write(fid2,'(/,A)') '>>> Export MOLCAS_NPROCS=1'
 ! currently MPI is not supported for the CASPT2 gradient in OpenMolcas-v23.02

 write(fid2,'(/,A)') "&CASPT2"
 if(dmrgscf) write(fid2,'(A)') ' CheMPS2'
 if(.not. caspt2_force) write(fid2,'(A)') ' IPEA= 0.25'
 ! OpenMolcas has an environment variable $MOLCAS_NEW_DEFAULTS, which will
 !  affect the default value of IPEA shift. So here we write it explicitly.
 write(fid2,'(A)') ' Frozen= 0'
 write(fid2,'(A,/,A)') ' MaxIter',' 500'

 if(caspt2_force) then
  write(fid2,'(A)') ' GRDT'
  write(fid2,'(A)') ' IPEA= 0.0'
  write(fid2,'(/,A,/)') "&ALASKA"
 else
  write(fid2,'(/)',advance='no')
 end if

 close(fid2)
 i = RENAME(inputname1, inputname)
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

 inpname1 = TRIM(inpname)//'.t'
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
end subroutine prt_mrmp2_gms_inp

! print OVB-MP2 keywords into Gaussian .gjf file
! Note: OVB-MP2 should be based on the converged CASSCF, not CASCI
subroutine prt_ovbmp2_gau_inp(gjfname)
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc, DKH2
 implicit none
 integer :: i, fid
 character(len=240), intent(in) :: gjfname

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 i = INDEX(gjfname, '.gjf', back=.true.)
 write(fid,'(A,I0)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A5,I0,A2)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(4(A,I0),A)',advance='no') '#p CASSCF(',nacte,',',nacto,') MP2&
 & chkbasis nosymm guess=read geom=allcheck iop(5/52=100)'
 ! IOp(5/52): configuration cutoff for mp2
 ! i     float(1/i)

 if(DKH2) then
  write(fid,'(A,/)') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A,/)') ' int=nobasistransform'
 end if

 close(fid)
end subroutine prt_ovbmp2_gau_inp

! read CASSCF OVB-MP2 energy from a Gaussian output file
subroutine read_mrpt_energy_from_gau_out(outname, ref_e, corr_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, corr_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(13:22) == 'EIGENVALUE') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_gau_out: no '&
                      &EIGENVALUE' found in file "//TRIM(outname)
  close(fid)
 end if
 read(buf(23:),*) ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(28:34) == 'EUMP2 =') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_gau_out: no '&
                      &EUMP2 =' found in file "//TRIM(outname)
 end if

 read(buf(35:),*) corr_e
 corr_e = corr_e - ref_e
end subroutine read_mrpt_energy_from_gau_out

! read nroots and target_root from a plain text file
subroutine read_ss_root_from_txt(nroots, target_root)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nroots, target_root
 character(len=50) :: buf
 character(len=10), parameter :: txtname = 'ss-cas.txt'

 nroots = 3; target_root = 1 ! initialization
 open(newunit=fid,file=txtname,status='old',position='rewind')

 read(fid,'(A)') buf
 i = INDEX(buf, '=')
 read(buf(i+1:),*) nroots

 read(fid,'(A)') buf
 i = INDEX(buf, '=')
 read(buf(i+1:),*) target_root

 close(fid)
end subroutine read_ss_root_from_txt

