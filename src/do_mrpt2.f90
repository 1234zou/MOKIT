! written by jxzou at 20210111: moved from subroutine do_mrpt2 in automr.f90

! do CASCI/CASPT2(npair<=7) or DMRG-CASCI/DMRG-NEVPT2 (npair>7) when scf=.False.
! do CASSCF/CASPT2(npair<=7) or DMRG-CASSCF/DMRG-NEVPT2 (npair>7) when scf=.True.
subroutine do_mrpt2()
 use print_id, only: iout
 use mr_keyword, only: casci, casscf, dmrgci, dmrgscf, CIonly, caspt2, &
  nevpt2, mrmp2, sdspt2, casnofch, casscf_prog, casci_prog, nevpt2_prog, &
  caspt2_prog, bgchg, chgname, mem, nproc, gms_path, gms_scr_path, check_gms_path,&
  molcas_path, molpro_path, orca_path, bdf_path
 use mol, only: casci_e, casscf_e, caspt2_e, nevpt2_e, mrmp2_e, sdspt2_e, &
                davidson_e, ptchg_e, nuc_pt_e
 use util_wrapper, only: mkl2gbw
 implicit none
 integer :: i, mem0, RENAME, system
 character(len=24) :: data_string
 character(len=240) :: string, pyname, outname, inpname, inporb
 character(len=240) :: mklname, fchname, cmofch
 real(kind=8) :: ref_e, corr_e
 logical :: alive(4)

 alive = [caspt2, nevpt2, mrmp2, sdspt2]
 if(ALL(alive .eqv. .false.)) return
 write(iout,'(//,A)') 'Enter subroutine do_mrpt2...'
 mem0 = CEILING(DBLE(mem*125)/DBLE(nproc))

 if((dmrgci .or. dmrgscf)) then
  if(mrmp2) then
   write(iout,'(A)') 'ERROR in subroutine do_mrpt2: DMRG-MRMP2 not supported.'
   stop
  end if
  if(sdspt2) then
   write(iout,'(A)') 'ERROR in subroutine do_mrpt2: DMRG-SDSPT2 not supported.'
   stop
  end if
  if(nevpt2 .and. (nevpt2_prog=='molpro' .or. nevpt2_prog=='orca')) then
   write(iout,'(A)') 'ERROR in subroutine do_mrpt2: DMRG-NEVPT2 is only supported&
                    & with PySCF or OpenMolcas.'
   write(iout,'(A)') 'But you specify NEVPT2_prog='//TRIM(nevpt2_prog)
   stop
  end if
  if(caspt2 .and. caspt2_prog/='openmolcas') then
   write(iout,'(A)') 'ERROR in subroutine do_mrpt2: DMRG-CASPT2 is only supported&
                    & with OpenMolcas.'
   write(iout,'(A)') 'But you specify CASPT2_prog='//TRIM(caspt2_prog)
   stop
  end if
 end if

 if(.not. CIonly) then
  if(casscf_prog == 'orca') then
   write(iout,'(A)') 'Warning: ORCA is used as the CASSCF solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(iout,'(A)') 'This will affect the PT2 energy up to 10^-5 a.u.&
                    & Such small error is usually not important.'
   write(iout,'(A)') 'If you care about the accuracy, please use another CASSCF solver.'
  end if

  if(casscf) then
   if(caspt2) then
    string = 'CASPT2 based on optimized CASSCF orbitals.'
   else if(nevpt2) then
    string = 'NEVPT2 based on optimized CASSCF orbitals.'
   else if(mrmp2) then
    string = 'MRMP2 based on optimized CASSCF orbitals.'
   else
    string = 'SDSMP2 based on optimized CASSCF orbitals.'
   end if
  else ! DMRG-CASSCF
   if(caspt2) then
    string = 'DMRG-CASPT2 based on optimized DMRG-CASSCF orbitals.'
   else if(nevpt2) then
    string = 'DMRG-NEVPT2 based on optimized DMRG-CASSCF orbitals.'
   end if
  end if

 else ! CIonly = .True.
  if(casci_prog == 'orca') then
   write(iout,'(A)') 'Warning: ORCA is used as the CASCI solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(iout,'(A)') 'This will affect the PT2 energy up to 10^-5 a.u.&
                    & Such small error is usually not important.'
   write(iout,'(A)') 'If you care about the accuracy, please use another CASCI solver.'
  end if

  if(casci) then
   if(caspt2) then
    string = 'CASPT2 based on CASCI orbitals.'
   else if(nevpt2) then
    string = 'NEVPT2 based on CASCI orbitals.'
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
 write(iout,'(A)',advance='no') 'Frozen_core = F, '

 if(nevpt2) then
  write(iout,'(A)') 'NEVPT2 using program '//TRIM(nevpt2_prog)

  select case(TRIM(nevpt2_prog))
  case('pyscf')
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
   i = system('python '//TRIM(pyname)//" >& "//TRIM(outname))

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
   write(string,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -m ',mem0, &
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
   i = system(TRIM(molcas_path)//' '//TRIM(inpname)//" >& "//TRIM(outname))

  case('orca')
   call check_exe_exist(orca_path)
   i = system('fch2mkl '//TRIM(casnofch))
   i = index(casnofch, '.fch', back=.true.)
   inporb = casnofch(1:i-1)//'.mkl'
   string = casnofch(1:i-1)//'.inp'
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
   i = system(TRIM(orca_path)//' '//TRIM(inpname)//" >& "//TRIM(outname))

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
   i = system('fch2inporb '//TRIM(casnofch)//' -no') ! generate .input and .INPORB
   i = index(casnofch, '.fch', back=.true.)
   pyname = casnofch(1:i-1)//'.input'
   outname = casnofch(1:i-1)//'.INPORB'
   i = index(casnofch, '_NO', back=.true.)
   inpname = casnofch(1:i-1)//'_CASPT2.input'
   inporb = casnofch(1:i-1)//'_CASPT2.INPORB'
   i = RENAME(pyname, inpname)
   i = RENAME(outname, inporb)
 
   i = index(casnofch, '_NO', back=.true.)
   outname = casnofch(1:i-1)//'_CASPT2.out'
   call prt_caspt2_molcas_inp(inpname)
   if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   i = system(TRIM(molcas_path)//' '//TRIM(inpname)//" >& "//TRIM(outname))

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
   write(string,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -m ',mem0, &
                            'm '//TRIM(inpname)
   i = system(TRIM(string))
  end select

 else if(mrmp2) then ! CASSCF-MRMP2
  write(iout,'(A)') 'MRMP2 using program gamess'
  call check_gms_path()

  i = system('fch2inp '//TRIM(casnofch))
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
  write(string,'(A,I0,A)') TRIM(gms_path)//' '//TRIM(inpname)//' 01 ', nproc, " >&"//TRIM(outname)
  i = system(TRIM(string))

  i = system('mv '//TRIM(mklname)//' .')

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
   write(iout,'(A)') 'ERROR in subroutine do_mrpt2: NEVPT2 computation failed.'
  else if(caspt2) then
   write(iout,'(A)') 'ERROR in subroutine do_mrpt2: CASPT2 computation failed.'
  else if(sdspt2) then
   write(iout,'(A)') 'ERROR in subroutine do_mrpt2: SDSPT2 computation failed.'
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
   call read_mrpt_energy_from_molpro_out(outname, 1, ref_e, corr_e)
   ref_e = ref_e + ptchg_e
  case('openmolcas')
   call read_mrpt_energy_from_molcas_out(outname, 1, ref_e, corr_e)
   ref_e = ref_e + ptchg_e
  case('orca')
   call read_mrpt_energy_from_orca_out(outname, 1, ref_e, corr_e)
  case('bdf')
   call read_mrpt_energy_from_bdf_out(outname, 2, ref_e, corr_e, davidson_e)
   ! here the parameter davidson_e is useless
   ref_e = ref_e + ptchg_e + nuc_pt_e
  end select
 else if(caspt2) then ! read CASPT2 energy
  select case(TRIM(caspt2_prog))
  case('openmolcas')
   call read_mrpt_energy_from_molcas_out(outname, 2, ref_e, corr_e)
  case('molpro')
   call read_mrpt_energy_from_molpro_out(outname, 2, ref_e, corr_e)
  end select
  ref_e = ref_e + ptchg_e
 else if(mrmp2) then  ! read MRMP2 energy
  call read_mrpt_energy_from_gms_out(outname, ref_e, corr_e)
 else                 ! read SDSPT2 energy
  call read_mrpt_energy_from_bdf_out(outname, 1, ref_e, corr_e, davidson_e)
 end if

 write(iout,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e, 'a.u.'

 if(caspt2) then
  write(iout,'(A)')             'IP-EA shift  =         0.25 (default)'
  caspt2_e = ref_e + corr_e
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(CASPT2)    = ', caspt2_e, 'a.u.'
 else if(nevpt2) then ! NEVPT2
  nevpt2_e = ref_e + corr_e
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,   'a.u.'
  if(nevpt2_prog == 'bdf') then
   write(iout,'(A,F18.8,1X,A4)')'E(FIC-NEVPT2)= ', nevpt2_e, 'a.u.'
  else
   write(iout,'(A,F18.8,1X,A4)')'E(SC-NEVPT2) = ', nevpt2_e, 'a.u.'
  end if
 else if(mrmp2) then ! MRMP2
  mrmp2_e = ref_e + corr_e
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)      = ', corr_e,  'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(MRMP2)     = ', mrmp2_e, 'a.u.'
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
 integer :: i, fid, nclosed, nocc
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

