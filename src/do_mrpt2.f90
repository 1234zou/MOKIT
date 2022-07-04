! written by jxzou at 20210111: moved from subroutine do_mrpt2 in automr.f90
! updated by jxzou at 20210224: add CASPT2 interface with ORCA

! do CASCI/CASPT2(npair<=7) or DMRG-CASCI/DMRG-NEVPT2 (npair>7) when scf=.False.
! do CASSCF/CASPT2(npair<=7) or DMRG-CASSCF/DMRG-NEVPT2 (npair>7) when scf=.True.
subroutine do_mrpt2()
 use print_id, only: iout
 use mr_keyword, only: casci, casscf, dmrgci, dmrgscf, CIonly, caspt2, caspt2k,&
  nevpt2, mrmp2, ovbmp2, sdspt2, casnofch, casscf_prog, casci_prog, nevpt2_prog, &
  caspt2_prog, bgchg, chgname, mem, nproc, gms_path, gms_scr_path, check_gms_path,&
  molcas_path, molpro_path, orca_path, bdf_path, gau_path, FIC, dryrun, eist
 use mol, only: caspt2_e, nevpt2_e, mrmp2_e, sdspt2_e, ovbmp2_e, davidson_e, &
  ptchg_e, nuc_pt_e
 use util_wrapper, only: mkl2gbw, fch2inp_wrap, unfchk
 implicit none
 integer :: i, mem0, RENAME, system
 character(len=24) :: data_string
 character(len=240) :: string, pyname, outname, inpname, inporb
 character(len=240) :: mklname, cmofch
 real(kind=8) :: ref_e, corr_e
 logical :: alive(5)

 if(eist == 1) return ! excited state calculation
 alive = [caspt2, nevpt2, mrmp2, ovbmp2, sdspt2]
 if(ALL(alive .eqv. .false.)) return
 write(iout,'(//,A)') 'Enter subroutine do_mrpt2...'
 mem0 = CEILING(DBLE(mem*125)/DBLE(nproc))

 if((dmrgci .or. dmrgscf)) then
  if(mrmp2) then
   write(iout,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-MRMP2 not supported.'
   stop
  end if
  if(sdspt2) then
   write(iout,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-SDSPT2 not supported.'
   stop
  end if
  if(ovbmp2) then
   write(iout,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-OVB-MP2 not supported.'
   stop
  end if
  if(nevpt2 .and. (nevpt2_prog=='molpro' .or. nevpt2_prog=='orca')) then
   write(iout,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-NEVPT2 is only supported&
                    & with PySCF or OpenMolcas.'
   write(iout,'(A)') 'But you specify NEVPT2_prog='//TRIM(nevpt2_prog)
   stop
  end if
  if(caspt2 .and. caspt2_prog/='openmolcas') then
   write(iout,'(/,A)') 'ERROR in subroutine do_mrpt2: DMRG-CASPT2 is only supported&
                    & with OpenMolcas.'
   write(iout,'(A)') 'But you specify CASPT2_prog='//TRIM(caspt2_prog)
   stop
  end if
 end if

 if(.not. CIonly) then
  if(casscf_prog == 'orca') then
   write(iout,'(A)') 'Warning: ORCA is used as the CASSCF solver, the NO&
                    & coefficients in .mkl file are'
   write(iout,'(A)') 'only 7 digits. This will affect the PT2 energy up to 10^-5&
                    & a.u. Such small error'
   write(iout,'(A)') 'is usually not important. If you care about the accuracy,&
                    & please use another CASSCF solver.'
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
  if(casci_prog == 'orca') then
   write(iout,'(A)') 'Warning: ORCA is used as the CASCI solver, the NO&
                    & coefficients in .mkl file are'
   write(iout,'(A)') 'only 7-digits. This will affect the PT2 energy up to 10^-5&
                    & a.u. Such small error'
   write(iout,'(A)') 'is usually not important. If you care about the accuracy,&
                    & please use another CASCI solver.'
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
  write(iout,'(A)') 'Warning: the CASSCF orbital optimization is strongly&
                     & recommended'
  write(iout,'(A)') 'to be performed before PT2, unless it is too time-consuming.'
 end if

 write(iout,'(A)') TRIM(string)
 if(dryrun) return
 write(iout,'(A)',advance='no') 'Frozen_core = F, '

 if(nevpt2) then
  write(iout,'(A)') 'NEVPT2 using program '//TRIM(nevpt2_prog)

  select case(TRIM(nevpt2_prog))
  case('pyscf')
   ! For DMRG-NEVPT2, use CMOs rather than NOs
   if(dmrgci .or. dmrgscf) then
    i = index(casnofch, '_NO', back=.true.)
    cmofch = casnofch(1:i)//'CMO.fch'
    casnofch = cmofch
   end if

   i = system('bas_fch2py '//TRIM(casnofch))
   i = index(casnofch, '.fch', back=.true.)
   inpname = casnofch(1:i-1)//'.py'
   if(dmrgci .or. dmrgscf) then
    i = index(casnofch, '_CMO', back=.true.)
   else
    i = index(casnofch, '_NO', back=.true.)
   end if
   pyname  = casnofch(1:i)//'NEVPT2.py'
   outname = casnofch(1:i)//'NEVPT2.out'
   i = RENAME(TRIM(inpname), TRIM(pyname))
   call prt_nevpt2_script_into_py(pyname)
   if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(pyname))
   i = system('python '//TRIM(pyname)//' >'//TRIM(outname)//" 2>&1")

  case('molpro')
   call check_exe_exist(molpro_path)

   i = system('fch2com '//TRIM(casnofch))
   i = index(casnofch, '.fch', back=.true.)
   pyname = casnofch(1:i-1)//'.com'
   string = casnofch
   call convert2molpro_fname(string, '.a')
   i = index(casnofch, '_NO', back=.true.)
   inpname = casnofch(1:i-1)//'_NEVPT2.com'
   outname = casnofch(1:i-1)//'_NEVPT2.out'
   inporb = inpname
   call convert2molpro_fname(inporb, '.a')
   i = RENAME(TRIM(pyname), TRIM(inpname))
   i = RENAME(TRIM(string), TRIM(inporb))

   call prt_mrpt_molpro_inp(inpname, 1)
   if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   write(string,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -t 1 -m ',mem0,&
                            'm '//TRIM(inpname)
   i = system(TRIM(string))

  case('openmolcas')
   call check_exe_exist(molcas_path)
   i = system('fch2inporb '//TRIM(casnofch)//' -no') ! generate .input and .INPORB
   i = index(casnofch, '.fch', back=.true.)
   pyname  = casnofch(1:i-1)//'.input'
   mklname = casnofch(1:i-1)//'.INPORB'
   i = index(casnofch, '_NO', back=.true.)
   inpname = casnofch(1:i)//'NEVPT2.input'
   inporb  = casnofch(1:i)//'NEVPT2.INPORB'
   outname = casnofch(1:i)//'NEVPT2.out'
   i = RENAME(TRIM(pyname), TRIM(inpname))
   i = RENAME(TRIM(mklname), TRIM(inporb))
 
   call prt_nevpt2_molcas_inp(inpname)
   if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   i = system(TRIM(molcas_path)//' '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1")

  case('orca')
   call check_exe_exist(orca_path)
   i = system('fch2mkl '//TRIM(casnofch))
   i = index(casnofch, '.fch', back=.true.)
   inporb = casnofch(1:i-1)//'_o.mkl'
   string = casnofch(1:i-1)//'_o.inp'
   i = index(casnofch, '_NO', back=.true.)
   mklname = casnofch(1:i)//'NEVPT2.mkl'
   inpname = casnofch(1:i)//'NEVPT2.inp'
   outname = casnofch(1:i)//'NEVPT2.out'
   i = RENAME(TRIM(inporb), TRIM(mklname))
   i = RENAME(TRIM(string), TRIM(inpname))

   call prt_nevpt2_orca_inp(inpname)
   if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   ! if bgchg = .True., .inp and .mkl file will be updated
   call mkl2gbw(mklname)
   call delete_file(mklname)
   call submit_orca_job(orca_path, inpname)

  case('bdf')
   call check_exe_exist(bdf_path)
   i = system('fch2bdf '//TRIM(casnofch)//' -no')
   i = index(casnofch, '.fch', back=.true.)
   inporb = casnofch(1:i-1)//'_bdf.inporb'
   string = casnofch(1:i-1)//'_bdf.inp'
   i = index(casnofch, '_NO', back=.true.)
   mklname = casnofch(1:i)//'NEVPT2.inporb'
   inpname = casnofch(1:i)//'NEVPT2.inp'
   outname = casnofch(1:i)//'NEVPT2.out'
   i = RENAME(TRIM(inporb), TRIM(mklname))
   i = RENAME(TRIM(string), TRIM(inpname))

   call prt_mrpt_bdf_inp(inpname, 2)
   if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   i = index(inpname, '.inp', back=.true.)
   string = inpname(1:i-1)
   i = system(TRIM(bdf_path)//' '//TRIM(string))
  end select

 else if(caspt2) then ! CASPT2
  write(iout,'(A)') 'CASPT2 using program '//TRIM(caspt2_prog)

  select case(TRIM(caspt2_prog))
  case('openmolcas')
   call check_exe_exist(molcas_path)
   ! For DMRG-CASPT2, use CMOs rather than NOs
   if(dmrgci .or. dmrgscf) then
    i = index(casnofch, '_NO', back=.true.)
    cmofch = casnofch(1:i)//'CMO.fch'
    casnofch = cmofch
   end if

   i = system('fch2inporb '//TRIM(casnofch)//' -no') ! generate .input and .INPORB
   i = index(casnofch, '.fch', back=.true.)
   pyname = casnofch(1:i-1)//'.input'
   string = casnofch(1:i-1)//'.INPORB'
   if(dmrgci .or. dmrgscf) then
    i = index(casnofch, '_CMO', back=.true.)
   else
    i = index(casnofch, '_NO', back=.true.)
   end if
   inpname = casnofch(1:i-1)//'_CASPT2.input'
   inporb  = casnofch(1:i-1)//'_CASPT2.INPORB'
   outname = casnofch(1:i-1)//'_CASPT2.out'
   i = RENAME(pyname, inpname)
   i = RENAME(string, inporb)
 
   call prt_caspt2_molcas_inp(inpname)
   if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   i = system(TRIM(molcas_path)//' '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1")

  case('molpro')
   call check_exe_exist(molpro_path)

   i = system('fch2com '//TRIM(casnofch))
   i = index(casnofch, '.fch', back=.true.)
   pyname = casnofch(1:i-1)//'.com'
   string = casnofch
   call convert2molpro_fname(string, '.a')
   i = index(casnofch, '_NO', back=.true.)
   inpname = casnofch(1:i-1)//'_CASPT2.com'
   outname = casnofch(1:i-1)//'_CASPT2.out'
   inporb = inpname
   call convert2molpro_fname(inporb, '.a')
   i = RENAME(TRIM(pyname), TRIM(inpname))
   i = RENAME(TRIM(string), TRIM(inporb))

   call prt_mrpt_molpro_inp(inpname, 2)
   if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   write(string,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -t 1 -m ',mem0,&
                            'm '//TRIM(inpname)
   i = system(TRIM(string))

  case('orca')
   call check_exe_exist(orca_path)

   i = system('fch2mkl '//TRIM(casnofch))
   i = index(casnofch, '.fch', back=.true.)
   inporb = casnofch(1:i-1)//'_o.mkl'
   string = casnofch(1:i-1)//'_o.inp'
   i = index(casnofch, '_NO', back=.true.)
   mklname = casnofch(1:i)//'CASPT2.mkl'
   inpname = casnofch(1:i)//'CASPT2.inp'
   outname = casnofch(1:i)//'CASPT2.out'
   i = RENAME(TRIM(inporb), TRIM(mklname))
   i = RENAME(TRIM(string), TRIM(inpname))

   call prt_caspt2_orca_inp(inpname)
   if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   ! if bgchg = .True., .inp and .mkl file will be updated
   call mkl2gbw(mklname)
   call delete_file(mklname)
   call submit_orca_job(orca_path, inpname)
  end select

 else if(mrmp2) then ! CASSCF-MRMP2
  write(iout,'(A)') 'MRMP2 using program gamess'
  call check_gms_path()

  call fch2inp_wrap(casnofch, .false., 0, 0)
  i = index(casnofch, '.fch')
  pyname = casnofch(1:i-1)//'.inp'
  i = index(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i-1)//'_MRMP2.inp'

  mklname = casnofch(1:i-1)//'_MRMP2.dat'
  mklname = TRIM(gms_scr_path)//'/'//TRIM(mklname) ! delete the possible .dat file
  call delete_file(mklname)

  outname = casnofch(1:i-1)//'_MRMP2.gms'
  i = RENAME(pyname, inpname)
  call prt_mrmp2_gms_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  write(string,'(A,I0,A)') TRIM(gms_path)//' '//TRIM(inpname)//' 01 ',nproc,&
                           ' >'//TRIM(outname)//" 2>&1"
  i = system(TRIM(string))
  i = system('mv '//TRIM(mklname)//' .')

 else if(ovbmp2) then ! OVB-MP2
  write(iout,'(A)') 'OVB-MP2 using program gaussian'
  call check_exe_exist(gau_path)
  i = index(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'OVBMP2.chk'
  inpname = casnofch(1:i)//'OVBMP2.gjf'
  outname = casnofch(1:i)//'OVBMP2.log'
  call unfchk(casnofch, mklname)
  call prt_ovbmp2_gau_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  string = TRIM(gau_path)//' '//TRIM(inpname)
  write(iout,'(A)') '$'//TRIM(string)
  i = system(TRIM(string))

 else ! CASSCF-SDSPT2
  write(iout,'(A)') 'SDSPT2 using program bdf'
  call check_exe_exist(bdf_path)

  i = system('fch2bdf '//TRIM(casnofch)//' -no')
  i = index(casnofch, '.fch', back=.true.)
  inporb = casnofch(1:i-1)//'_bdf.inporb'
  string = casnofch(1:i-1)//'_bdf.inp'
  i = index(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'SDSPT2.inporb'
  inpname = casnofch(1:i)//'SDSPT2.inp'
  outname = casnofch(1:i)//'SDSPT2.out'
  i = RENAME(TRIM(inporb), TRIM(mklname))
  i = RENAME(TRIM(string), TRIM(inpname))

  call prt_mrpt_bdf_inp(inpname, 1)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  i = index(inpname, '.inp', back=.true.)
  string = inpname(1:i-1)
  i = system(TRIM(bdf_path)//' '//TRIM(string))
 end if

 if(i /= 0) then
  if(nevpt2) then
   write(iout,'(/,A)') 'ERROR in subroutine do_mrpt2: NEVPT2 computation failed.'
  else if(caspt2) then
   write(iout,'(/,A)') 'ERROR in subroutine do_mrpt2: CASPT2 computation failed.'
  else if(ovbmp2) then
   write(iout,'(/,A)') 'ERROR in subroutine do_mrpt2: OVB-MP2 computation failed.'
  else if(sdspt2) then
   write(iout,'(/,A)') 'ERROR in subroutine do_mrpt2: SDSPT2 computation failed.'
  end if
  write(iout,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
  stop
 end if

 if(nevpt2) then      ! read NEVPT2 energy
  select case(TRIM(nevpt2_prog))
  case('pyscf')
   call read_mrpt_energy_from_pyscf_out(outname, ref_e, corr_e)
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

 write(iout,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e, 'a.u.'

 if(caspt2) then
  if(.not.caspt2k) write(iout,'(A)') 'IP-EA shift  =         0.25 (default)'
  caspt2_e = ref_e + corr_e
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e, 'a.u.'
  if(caspt2k) then
   write(iout,'(A,F18.8,1X,A4)') 'E(CASPT2-K)  = ', caspt2_e, 'a.u.'
  else
   write(iout,'(A,F18.8,1X,A4)') 'E(CASPT2)    = ', caspt2_e, 'a.u.'
  end if
 else if(nevpt2) then ! NEVPT2
  nevpt2_e = ref_e + corr_e
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,   'a.u.'
  if(FIC) then
   write(iout,'(A,F18.8,1X,A4)')'E(FIC-NEVPT2)= ', nevpt2_e, 'a.u.'
  else
   write(iout,'(A,F18.8,1X,A4)')'E(SC-NEVPT2) = ', nevpt2_e, 'a.u.'
  end if
 else if(mrmp2) then ! MRMP2
  mrmp2_e = ref_e + corr_e
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,  'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(MRMP2)     = ', mrmp2_e, 'a.u.'
 else if(ovbmp2) then ! OVB-MP2
  ovbmp2_e = ref_e + corr_e
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,  'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(OVB-MP2)   = ', ovbmp2_e,'a.u.'
 else           ! SDSPT2
  sdspt2_e = ref_e + corr_e + davidson_e
  write(iout,'(A,F18.8,1X,A4)') 'Davidson correction=',davidson_e,'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,  'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(SDSPT2)    = ', sdspt2_e, 'a.u.'
 end if

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_mrpt2 at '//TRIM(data_string)
 return
end subroutine do_mrpt2

! print NEVPT2/CASPT2/CASPT3 keywords into Molpro input file
subroutine prt_mrpt_molpro_inp(inpname, itype)
 use print_id, only: iout
 use mr_keyword, only: CIonly
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
  close(fid)
  write(iout,'(A)') 'ERROR in subroutine prt_mrpt_molpro_inp: wrong content found&
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
  write(fid,'(A)') '{RS2C,IPEA=0.25;CORE}'
 case(3) ! CASPT3
  write(fid,'(A)') '{RS3,IPEA=0.25;CORE}'
 case default
  write(iout,'(A,I0)') 'ERROR in subroutine prt_mrpt_molpro_inp: invalid&
                      & itype=',itype
  write(iout,'(A,I0)') 'Valid itype=1/2/3 for NEVPT2/CASPT2/CASPT3.'
  close(fid)
  stop
 end select

 close(fid)
 return
end subroutine prt_mrpt_molpro_inp

! print SDSPT2/NEVPT2/NEVPT3 keywords in to a given BDF .inp file
subroutine prt_mrpt_bdf_inp(inpname, itype)
 use print_id, only: iout
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
  close(fid)
  write(iout,'(A)') "ERROR in subroutine prt_mrpt_bdf_inp: no '$SCF' found&
                   & in file "//TRIM(inpname)
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
  write(iout,'(A,I0)') 'ERROR in subroutine prt_mrpt_bdf_inp: invalid itype=',itype
  write(iout,'(A)') 'itype=1/2/3 for SDSPT2/NEVPT2/NEVPT3.'
  stop
 end select

 write(fid,'(A)') '$END'
 close(fid)
 return
end subroutine prt_mrpt_bdf_inp

! print NEVPT2 keywords in to a given ORCA .inp file
subroutine prt_nevpt2_orca_inp(inpname)
 use print_id, only: iout
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc, X2C, CIonly, RI, RIJK_bas, F12, F12_cabs, &
  FIC, DLPNO, hardwfn, crazywfn
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
 write(fid2,'(A)',advance='no') '!'
 if(CIonly) write(fid2,'(A)',advance='no') ' noiter'
 if(RI) then
  ! RIJK in CASSCF must be combined with CONVentional
  ! /C basis set cannot be used with /JK in CASSCF
  write(fid2,'(A)',advance='no') ' RIJK conv '//TRIM(RIJK_bas)
 end if
 if(F12) write(fid2,'(A)',advance='no') ' '//TRIM(F12_cabs)
 write(fid2,'(A)') ' TightSCF'

 if(X2C) then
  write(iout,'(A)') 'ERROR in subroutine prt_nevpt2_orca_inp: NEVPT2 with X2C&
                   & is not supported in ORCA.'
  write(iout,'(A)') 'You can specify NEVPT2_prog=Molpro or OpenMolcas.'
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
 return
end subroutine prt_nevpt2_orca_inp

! print CASPT2 keywords in to a given ORCA .inp file
subroutine prt_caspt2_orca_inp(inpname)
 use print_id, only: iout
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc, caspt2k, DKH2, X2C, CIonly, RI, &
  RIJK_bas, hardwfn, crazywfn
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
 write(fid2,'(A,I0,A)') '%maxcore ', CEILING(1000.0d0*DBLE(mem)/DBLE(nproc))

 read(fid1,'(A)') buf   ! skip '!' line
 write(fid2,'(A)',advance='no') '!'
 if(CIonly) write(fid2,'(A)',advance='no') ' noiter'
 if(RI) then
  ! RIJK in CASSCF must be combined with CONVentional
  ! /C basis set cannot be used with /JK in CASSCF
  write(fid2,'(A)',advance='no') ' RIJK conv '//TRIM(RIJK_bas)
 end if
 write(fid2,'(A)') ' TightSCF'

 if(DKH2) then
  write(fid2,'(A)') '%rel'
  write(fid2,'(A)') ' method DKH'
  write(fid2,'(A)') ' order 2'
  write(fid2,'(A)') 'end'
 else if(X2C) then
  write(iout,'(A)') 'ERROR in subroutine prt_nevpt2_orca_inp: CASPT2 with X2C&
                   & is not supported in ORCA.'
  write(iout,'(A)') 'You can specify CASPT2_prog=Molpro or OpenMolcas.'
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
 if(.not. caspt2k) write(fid2,'(A)') '  CASPT2_IPEAshift 0.25'
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
 return
end subroutine prt_caspt2_orca_inp

! print NEVPT2 script into a given .py file
subroutine prt_nevpt2_script_into_py(pyname)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, casci, casscf, maxM, X2C, RI, RIJK_bas, &
  hardwfn, crazywfn
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, pyname1
 character(len=240), intent(in) :: pyname

 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
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

 write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*500, ' # MB'

 if(casci .or. casscf) then
  write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ', mem*500, ' # MB'
 else
  write(fid2,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=', maxM, ')'
  write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ', CEILING(DBLE(mem)/2.0d0), ' # GB'
 end if

 call prt_hard_or_crazy_casci_pyscf(fid2, nacta-nactb, hardwfn, crazywfn)
 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.kernel()'
 if(casci .or. casscf) then
  write(fid2,'(/,A)') 'mrpt.NEVPT(mc).kernel()'
 else
  write(fid2,'(/,A,I0,A)') 'mrpt.NEVPT(mc).compress_approx(maxM=',maxM,').kernel()'
 end if
 close(fid2)

 i = RENAME(pyname1, pyname)
 return
end subroutine prt_nevpt2_script_into_py

! print NEVPT2 keywords into OpenMolcas .input file
! It seems that OpenMolcas does not support CASSCF-NEVPT2. So I have to use
! DMRG-NEVPT2.
subroutine prt_nevpt2_molcas_inp(inpname)
 use mr_keyword, only: CIonly, maxM, RI, RIJK_bas
 use mol, only: nacte, nacto, charge, mult
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
  write(6,'(A)') "ERROR in subroutine prt_nevpt2_molcas_inp: no 'SEWARD'&
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
 use mr_keyword, only: CIonly, maxM, dmrgci, dmrgscf, hardwfn, crazywfn
 use mol, only: nacte, nacto, charge, mult
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
  if(buf(1:4) == "&SCF") exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while
 close(fid1,status='delete')

 call prt_molcas_cas_para(fid2, dmrg, .false., dmrg, CIonly, inputname)
 !write(fid2,'(A,I0,A)') 'nActEl= ', nacte, ' 0 0'

 if(CIonly) then
  if(dmrgci) then
   write(6,'(A)') 'ERROR in subroutine prt_caspt2_molcas_inp: CIonly is not&
                    & allowed in DMRG-CASPT2.'
   write(6,'(A)') 'It must be based on converged (DMRG-)CASSCF orbitals.'
   stop
  end if
 end if

 write(fid2,'(/,A)') "&CASPT2"
 write(fid2,'(A)') 'MultiState= 1 1'
 if(dmrgscf) write(fid2,'(A)') 'CheMPS2'
 write(fid2,'(A)') 'Frozen= 0'
 write(fid2,'(A,/,A,/)') 'MaxIter','300'

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

! print OVB-MP2 keywords into Gaussian .gjf file
! Note: OVB-MP2 should be based on the converged CASSCF, not CASCI
subroutine prt_ovbmp2_gau_inp(gjfname)
 use mol, only: nacte, nacto
 use mr_keyword, only: mem, nproc, DKH2
 implicit none
 integer :: i, fid
 character(len=240), intent(in) :: gjfname

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 i = index(gjfname, '.gjf', back=.true.)
 write(fid,'(A,I0)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A5,I0,A2)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(4(A,I0),A)',advance='no') '#p CAS(',nacte,',',nacto,') MP2&
 & chkbasis nosymm guess=read geom=allcheck iop(5/52=100)'

 if(DKH2) then
  write(fid,'(A,/)') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A,/)') ' int=nobasistransform'
 end if

 close(fid)
 return
end subroutine prt_ovbmp2_gau_inp

! read CASSCF OVB-MP2 energy from a Gaussian output file
subroutine read_mrpt_energy_from_gau_out(outname, ref_e, corr_e)
 use print_id, only: iout
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
  write(iout,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_gau_out: no '&
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
  write(iout,'(/,A)') "ERROR in subroutine read_mrpt_energy_from_gau_out: no '&
                      &EUMP2 =' found in file "//TRIM(outname)
 end if

 read(buf(35:),*) corr_e
 corr_e = corr_e - ref_e
 return
end subroutine read_mrpt_energy_from_gau_out

