! written by jxzou at 20210120: subroutine do_mcpdft from automr.f90 is moved here
! TODO: interfaces of SS-/SA-/CMS- MC-PDFT variants

! do MC-PDFT, valid for CASCI/CASSCF based MC-PDFT, and DMRG-PDFT
subroutine do_mcpdft()
 use mr_keyword, only: mem, nproc, casci, dmrgci, dmrgscf, mcpdft, mcpdft_prog,&
  casnofch, molcas_omp, molcas_path, gms_path, bgchg, chgname, check_gms_path, &
  gms_scr_path, n_otpdf, otpdf, mcpdft_force, eist
 use mol, only: nacte, nacto, ptchg_e, mcpdft_e, natom, grad
 use util_wrapper, only: bas_fch2py_wrap, fch2inp_wrap, fch2inporb_wrap
 implicit none
 integer :: i, SYSTEM, RENAME
 real(kind=8) :: ref_e
 real(kind=8), allocatable :: otpdf_e(:)
 character(len=10), allocatable :: split_otpdf(:)
 character(len=24) :: data_string
 character(len=240) :: fname(2), inpname, outname, cmofch
 logical :: dmrg

 if(eist == 1) return ! excited state calculation
 if(.not. mcpdft) return
 write(6,'(//,A)') 'Enter subroutine do_mcpdft...'
 dmrg = (dmrgci .or. dmrgscf)

 if(dmrg) then
  if(mcpdft_prog == 'gamess') then
   write(6,'(/,A)') 'ERROR in subroutine do_mcpdft: DMRG-PDFT has not been impl&
                    &emented in GAMESS.'
   write(6,'(A)') 'You can set MCPDFT_prog=PySCF or OpenMolcas in mokit{}.'
   stop
  end if
  if(dmrgci) then
   write(6,'(A)') 'DMRG-PDFT based on DMRG-CASCI orbitals.'
  else
   write(6,'(A)') 'DMRG-PDFT based on DMRG-CASSCF orbitals.'
  end if
  write(6,'(A,2(I0,A))') 'DMRG-PDFT(', nacte, 'e,', nacto, &
                         'o) using program OpenMolcas'
 else
  if(casci) then
   write(6,'(A)') 'MC-PDFT based on CASCI orbitals.'
  else
   write(6,'(A)') 'MC-PDFT based on CASSCF orbitals.'
  end if
  write(6,'(A,2(I0,A))') 'MC-PDFT(',nacte,'e,',nacto,'o) using program '//&
                         TRIM(mcpdft_prog)
 end if

 if(dmrgci .or. casci) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning: orbital optimization is strongly recommended to be&
                 & performed before'
  write(6,'(A)') 'PDFT, unless it is too time-consuming or the input orbitals&
                 & have been optimized.'
  write(6,'(A)') REPEAT('-',79)
 end if

 ! For DMRG-PDFT, use input orbitals or pseudo-canonical MOs rather than NOs
 if(dmrg) then
  i = INDEX(casnofch, '_NO', back=.true.)
  cmofch = casnofch(1:i)//'CMO.fch'
  casnofch = cmofch
 end if

 ! split multiple on-top pair density functionals, n_otpdf >= 1
 allocate(otpdf_e(n_otpdf), split_otpdf(n_otpdf))
 call split_str_otpdf(otpdf, n_otpdf, split_otpdf)

 select case(TRIM(mcpdft_prog))
 case('pyscf')
  call find_specified_suffix(casnofch, '_', i)
  inpname = casnofch(1:i)//'MCPDFT.py'
  outname = casnofch(1:i)//'MCPDFT.out'
  call bas_fch2py_wrap(casnofch, .false., inpname)
  call prt_mcpdft_script_into_py(inpname, n_otpdf, split_otpdf)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_pyscf_job(inpname, .true.)
  call read_mcpdft_e_from_pyscf_out(outname, n_otpdf, ref_e, otpdf_e)

 case('openmolcas')
  call check_exe_exist(molcas_path)
  call find_specified_suffix(casnofch, '_', i)
  inpname = casnofch(1:i)//'MCPDFT.input'
  outname = casnofch(1:i)//'MCPDFT.out'
  call fch2inporb_wrap(casnofch, .false., inpname)
  call prt_mcpdft_molcas_inp(inpname, n_otpdf, split_otpdf)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)
  call read_mcpdft_e_from_molcas_out(outname, n_otpdf, ref_e, otpdf_e)

 case('gamess')
  if(n_otpdf > 1) then
   write(6,'(/,A)') 'ERROR in subroutine do_mcpdft: multiple on-top pair densit&
                    &y functional is'
   write(6,'(A)') 'not supported when MCPDFT_prog=GAMESS. You can use MCPDFT_pr&
                  &og=PySCF or OpenMolcas.'
   stop
  end if
  call check_gms_path()
  call fch2inp_wrap(casnofch, .false., 0, 0, .false.)
  call find_specified_suffix(casnofch, '.fch', i)
  fname(1) = casnofch(1:i-1)//'.inp'
  i = INDEX(casnofch, '_NO', back=.true.)
  fname(2) = casnofch(1:i)//'MCPDFT.dat'
  inpname = casnofch(1:i)//'MCPDFT.inp'
  outname = casnofch(1:i)//'MCPDFT.gms'
  i = RENAME(TRIM(fname(1)), TRIM(inpname))

  call prt_mcpdft_gms_inp(inpname, split_otpdf(1))
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! MC-PDFT in GAMESS cannot run in parallel currently, use 1 core
  call submit_gms_job(gms_path, gms_scr_path, inpname, 1)
  call read_mcpdft_e_from_gms_gms(outname, ref_e, otpdf_e(1))

 case default
  write(6,'(/,A)') 'ERROR in subroutine do_mcpdft: wrong MCPDFT_prog='//&
                   TRIM(mcpdft_prog)
  stop
 end select

 if(bgchg .and. TRIM(mcpdft_prog)/='gamess') then
  ref_e = ref_e + ptchg_e
  otpdf_e = otpdf_e + ptchg_e
 end if
 mcpdft_e = otpdf_e(1)

 write(6,'(/,A,F18.8,A)')'E(ref)      = ',    ref_e, ' a.u.'
 ! the energy of the 1st on-top pair density function is printed here
 if(dmrg) then
  write(6,'(A,F18.8,A)') 'E(DMRG-PDFT)= ', mcpdft_e, ' a.u.'
 else
  write(6,'(A,F18.8,A)') 'E(MC-PDFT)  = ', mcpdft_e, ' a.u.'
 end if

 do i = 1, n_otpdf, 1
  call upper(split_otpdf(i))
  write(6,'(A,F18.8,A)') 'E('//TRIM(split_otpdf(i))//') = ',otpdf_e(i),' a.u.'
 end do ! for i

 if(mcpdft_force) then
  allocate(grad(3*natom))

  select case(TRIM(mcpdft_prog))
  case('pyscf')
   call read_grad_from_pyscf_out(outname, natom, grad)
  case('gamess')
   call read_grad_from_dat(fname(2), natom, grad)
  case('openmolcas')
   call read_grad_from_molcas_out(outname, natom, grad)
  end select

  write(6,'(A)') 'Cartesian gradient (HARTREE/BOHR):'
  write(6,'(5(1X,ES15.8))') (grad(i),i=1,3*natom)
 end if

 deallocate(otpdf_e, split_otpdf)
 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mcpdft at '//TRIM(data_string)
end subroutine do_mcpdft

subroutine split_str_otpdf(otpdf, n_otpdf, split_otpdf)
 implicit none
 integer :: i, k
 integer, intent(in) :: n_otpdf
 character(len=60) :: buf
 character(len=60), intent(in) :: otpdf
 character(len=10), intent(out) :: split_otpdf(n_otpdf)

 if(n_otpdf == 1) then
  split_otpdf(1) = TRIM(otpdf)
  return
 end if

 buf = otpdf
 k = 0
 i = INDEX(otpdf, ';')

 do while(i > 0)
  k = k + 1
  split_otpdf(k) = buf(1:i-1)
  buf = TRIM(buf(i+1:))
  i = INDEX(buf, ';')
 end do ! while

 k = k + 1
 split_otpdf(k) = TRIM(buf)
end subroutine split_str_otpdf

! print MC-PDFT or DMRG-PDFT keywords into PySCF .py file
subroutine prt_mcpdft_script_into_py(inpname, n_otpdf, split_otpdf)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, dmrgci, dmrgscf, mcpdft_force
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: n_otpdf
 character(len=10), intent(in) :: split_otpdf(n_otpdf)
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical :: dmrg

 dmrg = (dmrgci .or. dmrgscf)
 call find_specified_suffix(inpname, '.py', i)
 inpname1 = inpname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do i = 1, 3
  read(fid,'(A)') buf
  if(buf(1:17) == 'from pyscf import') then
   buf = TRIM(buf)//', mcpdft'
   write(fid1,'(A)') TRIM(buf)
   exit
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for i

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(/,A,I0)') 'nproc = ', nproc
 write(fid1,'(A)') 'lib.num_threads(nproc)'

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:15) == 'lib.num_threads') cycle
  if(buf(1:13) == 'mf.max_memory') cycle
  if(buf(1:9) == 'mf.kernel') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(A,I0,A)') 'mf.max_memory = ', mem*1000, ' # MB'
 write(fid1,'(A)') TRIM(buf)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while
 close(fid,status='delete')

 write(fid1,'(3(A,I0),A)') "mc = mcpdft.CASCI(mf,'"//TRIM(split_otpdf(1))//&
                           "',",nacto,',(',nacta,',',nactb,'))'
 write(fid1,'(A)') 'mc.grids.atom_grid = (99,590) # ultrafine'
 write(fid1,'(A,I0,A)') 'mc.max_memory = ',mem*1000,' # MB'
 write(fid1,'(A)') 'mc.verbose = 5'
 write(fid1,'(A)') 'mc.kernel()'

 do i = 2, n_otpdf, 1
  write(fid1,'(A)') "mc.compute_pdft_energy_(otxc='"//TRIM(split_otpdf(i))//"')"
 end do ! for i

 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 if(mcpdft_force) call add_force_key2py_script(mem, inpname, .false.)
end subroutine prt_mcpdft_script_into_py

! print MC-PDFT or DMRG-PDFT keywords into OpenMolcas .input file
subroutine prt_mcpdft_molcas_inp(inpname, n_otpdf, split_otpdf)
 use mr_keyword, only: RI, dmrgci, dmrgscf, DKH2, mcpdft_force
 implicit none
 integer :: i, fid1, fid2, RENAME
 integer, intent(in) :: n_otpdf
 character(len=6) :: key
 character(len=10), intent(in) :: split_otpdf(n_otpdf)
 character(len=10), allocatable :: new_otpdf(:)
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical :: dmrg

 allocate(new_otpdf(n_otpdf), source=split_otpdf)
 call detect_and_update_otpdf(n_otpdf, new_otpdf, key)

 if(RI) call add_RI_kywd_into_molcas_inp(inpname, .true.)
 dmrg = (dmrgci .or. dmrgscf)

 call find_specified_suffix(inpname, '.input', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
  if(buf(2:7) == 'SEWARD') exit
 end do ! for while
 close(fid1,status='delete')

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine prt_mcpdft_molcas_inp: no 'SEWARD' foun&
                   &d in file "//TRIM(inpname)
  close(fid2)
  stop
 end if

 write(fid2,'(A)') 'Grid input'
 write(fid2,'(A)') ' grid=ultrafine'
 write(fid2,'(A)') 'End of grid input'
 if(DKH2) write(fid2,'(A)') 'Relativistic = R02O'

 ! I used to use &RASSCF, DMRG, RGinput for DMRGCI calculations. I found it
 ! seems that RGinput in &RASSCF had many disadvantages and often converged
 ! to higher energy local minimum. So I changed in into &DMRGSCF, and it
 ! seems to get more reasonable results

 call prt_molcas_cas_para(fid2, dmrg, .false., .false., .true., inpname)

 do i = 1, n_otpdf, 1
  write(fid2,'(/,A)') "&MCPDFT"
  write(fid2,'(A)') TRIM(key)//TRIM(new_otpdf(i))
  if(i==1 .and. mcpdft_force) then
   write(fid2,'(A)') 'Grad'
   write(fid2,'(/,A)') '&ALASKA'
  end if
 end do ! for i

 deallocate(new_otpdf)
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_mcpdft_molcas_inp

! print MC-PDFT keywords into GAMESS .inp file
subroutine prt_mcpdft_gms_inp(inpname, otpdf)
 use mol, only: charge, mult, ndb, nacte, nacto
 use mr_keyword, only: mem, cart, DKH2, hardwfn, crazywfn, CIonly, mcpdft_force
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=10), intent(in) :: otpdf
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'.t'

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

 ! MC-PDFT in GAMESS cannot run in parallel currently. All memory is assigned to
 ! 1 core.
 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',mem*125,' $END'

 if(CIonly) then
  write(fid2,'(A)',advance='no') ' $CIDET'
 else
  write(fid2,'(A)',advance='no') ' $DET'
 end if
 write(fid2,'(3(A,I0))',advance='no') ' NCORE=',ndb,' NELS=',nacte,' NACT=',nacto

 if(hardwfn) then
  write(fid2,'(A)',advance='no') ' NSTATE=5'
 else if(crazywfn) then
  write(fid2,'(A)',advance='no') ' NSTATE=10'
 end if
 write(fid2,'(A)') ' ITERMX=500 $END'
 write(fid2,'(A)') ' $DFT NRAD=99 NLEB=590 $END'

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
 if(mcpdft_force) call add_force_key2gms_inp(inpname)
end subroutine prt_mcpdft_gms_inp

! detect the OpenMolcas version and update otpdf if needed
! OpenMolcas-v22.xx: T:PBE, FT:PBE
! OpenMolcas-v21.xx: tPBE, ftPBE
subroutine detect_and_update_otpdf(n_otpdf, split_otpdf, key)
 implicit none
 integer :: i, j, k, iv, fid, SYSTEM
 integer, intent(in) :: n_otpdf
 character(len=6), intent(out) :: key
 character(len=10) :: str
 character(len=10), intent(inout) :: split_otpdf(n_otpdf)
 character(len=30) :: ftmp
 character(len=240) :: buf

 key = 'KSDFT='
 ! Find the OpenMolcas version. Note: `pymolcas --banner` cannot be used at 1st,
 ! Apr each year. So we use `pymolcas version` instead.
 write(ftmp,'(A,I0)') 'molcas.ver.', i
 i = SYSTEM('pymolcas version >'//TRIM(ftmp)//" 2>&1")

 iv = 0
 open(newunit=fid,file=TRIM(ftmp),status='old',position='rewind')
 read(fid,'(A)',iostat=i) buf
 close(fid,status='delete')
 if(i /= 0) return
 ! if something like 'version: v22.06' not found, just return
 buf = ADJUSTL(buf)
 i = INDEX(buf, 'v')
 j = INDEX(buf, '.')
 read(buf(i+1:j-1),fmt=*,iostat=k) iv ! 22
 if(k /= 0) return
 ! if the user uses some strange version of OpenMolcas, just return

 if(iv >= 22) then ! tPBE, ftPBE -> T:PBE, FT:PBE
  key = 'FUNC= '

  do i = 1, n_otpdf, 1
   str = split_otpdf(i)
   if(LEN_TRIM(str)==0 .or. TRIM(str)=='NONE') then
    write(6,'(/,A)') 'ERROR in subroutine detect_and_update_otpdf: invalid OtPD&
                     &F='//TRIM(str)
    write(6,'(A)') 'Something must be wrong.'
    stop
   end if

   call lower(str)
   if(str(1:1) == 't') then
    split_otpdf(i) = 'T:'//TRIM(split_otpdf(i)(2:))
   else if(str(1:2) == 'ft') then
    split_otpdf(i) = 'FT:'//TRIM(split_otpdf(i)(3:))
   else
    write(6,'(/,A)') 'ERROR in subroutine detect_and_update_otpdf: wrong OtPDF=&
                     &'//TRIM(str)
    stop
   end if
  end do ! for i
 end if
end subroutine detect_and_update_otpdf

! read MC-PDFT energy from a given PySCF/OpenMolcas/GAMESS output file
subroutine read_mcpdft_e_from_pyscf_out(outname, n_otpdf, ref_e, pdft_e)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: n_otpdf
 real(kind=8), intent(out) :: ref_e, pdft_e(n_otpdf)
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0; pdft_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:7) == 'CASCI E') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mcpdft_e_from_pyscf_out: no 'CASCI&
                   & E' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) ref_e

 do i = 1, n_otpdf, 1
  do while(.true.)
   read(fid,'(A)',iostat=k) buf
   if(k /= 0) exit
   if(buf(1:9)=='MC-PDFT E' .or. buf(1:17)=='MC-PDFT state 0 E') exit
  end do ! for while

  if(k /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine read_mcpdft_e_from_pyscf_out: there is&
                    & no enough number'
   write(6,'(A)') 'of MC-PDFT energies in file '//TRIM(outname)
   close(fid)
   stop
  end if

  k = INDEX(buf, '=')
  read(buf(k+1:),*) pdft_e(i)
 end do ! for i

 close(fid)
end subroutine read_mcpdft_e_from_pyscf_out

! read MC-PDFT energy from an OpenMolcas output file
subroutine read_mcpdft_e_from_molcas_out(outname, n_otpdf, ref_e, pdft_e)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: n_otpdf
 real(kind=8), intent(out) :: ref_e, pdft_e(n_otpdf)
 character(len=240) :: buf
 character(len=9) :: str(3)
 character(len=240), intent(in) :: outname

 ref_e = 0d0; pdft_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'MCSCF reference e') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mcpdft_e_from_molcas_out: no 'MCSC&
                   &F reference e'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if
 read(buf,*) str(1), str(2), str(3), ref_e

 do i = 1, n_otpdf, 1
  do while(.true.)
   read(fid,'(A)',iostat=k) buf
   if(k /= 0) exit
   if(INDEX(buf,'Total MC-PDFT')>0 .or. INDEX(buf,'Total MCPDFT')>0) exit
  end do ! for while

  if(k /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine read_mcpdft_e_from_molcas_out: there i&
                    &s no enough number'
   write(6,'(A)') 'of MC-PDFT energies in file '//TRIM(outname)
   close(fid)
   stop
  end if

  read(buf(60:),*) pdft_e(i)
 end do ! for i

 close(fid)
end subroutine read_mcpdft_e_from_molcas_out

! read MC-PDFT energy from a GAMESS output file (.gms)
subroutine read_mcpdft_e_from_gms_gms(gmsname, ref_e, pdft_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, pdft_e
 character(len=240) :: buf
 character(len=240), intent(in) :: gmsname

 ref_e = 0d0; pdft_e = 0d0
 open(newunit=fid,file=TRIM(gmsname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(INDEX(buf,'TOTAL MC-PDFT') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mcpdft_e_from_gms_gms: no 'Total M&
                   &C-PDFT' found in"
  write(6,'(A)') 'file '//TRIM(gmsname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) pdft_e

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:15) == 'STATE=   1   E') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mcpdft_e_from_gms_gms: no 'STATE= &
                   &  1   E' found in"
  write(6,'(A)') 'file '//TRIM(gmsname)
  close(fid)
  stop
 end if

 i = INDEX(buf, 'Y=', back=.true.)
 read(buf(i+2:),*) ref_e
 close(fid)
end subroutine read_mcpdft_e_from_gms_gms

