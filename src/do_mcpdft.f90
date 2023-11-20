! written by jxzou at 20210120: subroutine do_mcpdft from automr.f90 is moved here

! do MC-PDFT, valid for CASCI/CASSCF based MC-PDFT, and DMRG-PDFT
subroutine do_mcpdft()
 use mr_keyword, only: mem, nproc, casci, dmrgci, dmrgscf, mcpdft, mcpdft_prog,&
  casnofch, molcas_omp, molcas_path, gms_path, bgchg, chgname, check_gms_path, &
  gms_scr_path, mcpdft_force, eist
 use mol, only: ptchg_e, mcpdft_e, natom, grad
 use util_wrapper, only: bas_fch2py_wrap, fch2inp_wrap, fch2inporb_wrap
 implicit none
 integer :: i, system, RENAME
 real(kind=8) :: ref_e
 character(len=24) :: data_string
 character(len=240) :: fname(2), inpname, outname, cmofch
 logical :: dmrg

 if(eist == 1) return ! excited state calculation
 if(.not. mcpdft) return
 write(6,'(//,A)') 'Enter subroutine do_mcpdft...'

 dmrg = (dmrgci .or. dmrgscf)

 if(dmrg) then
  if(mcpdft_prog == 'gamess') then
   write(6,'(A)') 'ERROR in subroutine do_mcpdft: DMRG-PDFT has not been&
                  & implemented in GAMESS.'
   write(6,'(A)') 'You can set MCPDFT_prog=OpenMolcas in mokit{}.'
   stop
  end if
  if(dmrgci) then
   write(6,'(A)') 'DMRG-PDFT based on DMRG-CASCI orbitals.'
  else
   write(6,'(A)') 'DMRG-PDFT based on DMRG-CASSCF orbitals.'
  end if
  write(6,'(A)') 'DMRG-PDFT using program openmolcas'
 else
  if(casci) then
   write(6,'(A)') 'MC-PDFT based on CASCI orbitals.'
  else
   write(6,'(A)') 'MC-PDFT based on CASSCF orbitals.'
  end if
  write(6,'(A)') 'MC-PDFT using program '//TRIM(mcpdft_prog)
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

 select case(TRIM(mcpdft_prog))
 case('pyscf')
  i = INDEX(casnofch, '_', back=.true.)
  inpname = casnofch(1:i)//'MCPDFT.py'
  outname = casnofch(1:i)//'MCPDFT.out'
  call bas_fch2py_wrap(casnofch, .false., inpname)
  call prt_mcpdft_script_into_py(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_pyscf_job(inpname)

 case('openmolcas')
  call check_exe_exist(molcas_path)
  i = INDEX(casnofch, '_', back=.true.)
  inpname = casnofch(1:i)//'MCPDFT.input'
  outname = casnofch(1:i)//'MCPDFT.out'
  call fch2inporb_wrap(casnofch, .false., inpname)
  call prt_mcpdft_molcas_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)

 case('gamess')
  call check_gms_path()
  call fch2inp_wrap(casnofch, .false., 0, 0)
  i = INDEX(casnofch, '.fch', back=.true.)
  fname(1) = casnofch(1:i-1)//'.inp'
  i = INDEX(casnofch, '_NO', back=.true.)
  fname(2) = casnofch(1:i)//'MCPDFT.dat'
  inpname = casnofch(1:i)//'MCPDFT.inp'
  outname = casnofch(1:i)//'MCPDFT.gms'
  i = RENAME(TRIM(fname(1)), TRIM(inpname))

  call prt_mcpdft_gms_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_gms_job(gms_path, gms_scr_path, inpname, 1)
  ! MC-PDFT in GAMESS cannot run in parallel currently, use 1 core
 end select

 ! read MC-PDFT/DMRG-PDFT energy from output file
 call read_mcpdft_e_from_output(mcpdft_prog, outname, ref_e, mcpdft_e)
 if(bgchg .and. TRIM(mcpdft_prog)/='gamess') then
  ref_e = ref_e + ptchg_e
  mcpdft_e = mcpdft_e + ptchg_e
 end if

 write(6,'(/,A,F18.8,1X,A4)')'E(ref)      = ',    ref_e, 'a.u.'
 if(dmrg) then
  write(6,'(A,F18.8,1X,A4)') 'E(DMRG-PDFT)= ', mcpdft_e, 'a.u.'
 else
  write(6,'(A,F18.8,1X,A4)') 'E(MC-PDFT)  = ', mcpdft_e, 'a.u.'
 end if

 if(mcpdft_force) then
  allocate(grad(3*natom))

  select case(mcpdft_prog)
  case('pyscf')
   call read_grad_from_pyscf_out(outname, natom, grad)
  case('gamess')
   call read_grad_from_gms_dat(fname(2), natom, grad)
  case('openmolcas')
   call read_grad_from_molcas_out(outname, natom, grad)
  case default
   write(6,'(A)') 'ERROR in subroutine do_cas: program cannot be identified.'
   write(6,'(A)') 'MCPDFT_prog='//TRIM(mcpdft_prog)
   stop
  end select

  write(6,'(A)') 'Cartesian gradient (HARTREE/BOHR):'
  write(6,'(5(1X,ES15.8))') (grad(i),i=1,3*natom)
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mcpdft at '//TRIM(data_string)
end subroutine do_mcpdft

! print MC-PDFT or DMRG-PDFT keywords into PySCF .py file
subroutine prt_mcpdft_script_into_py(inpname)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: dmrgci, dmrgscf, otpdf, mem, mcpdft_force
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical :: dmrg

 dmrg = (dmrgci .or. dmrgscf)
 inpname1 = TRIM(inpname)//'.t'
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
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 write(fid1,'(3(A,I0),A)') "mc = mcpdft.CASCI(mf,'"//TRIM(otpdf)//"',",nacto, &
                           ',(',nacta,',',nactb,'))'
 write(fid1,'(A)') 'mc.grids.atom_grid = (99,590)'
 write(fid1,'(A,I0,A)') 'mc.max_memory = ',mem*1000,' # MB'
 write(fid1,'(A)') 'mc.kernel()'
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))

 if(mcpdft_force) call add_force_key2py_script(mem, inpname, .false.)
end subroutine prt_mcpdft_script_into_py

! print MC-PDFT or DMRG-PDFT keywords into OpenMolcas .input file
subroutine prt_mcpdft_molcas_inp(inpname)
 use mr_keyword, only: dmrgci, dmrgscf, otpdf, DKH2, mcpdft_force
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical :: dmrg

 ! detect the OpenMolcas version and update otpdf if needed
 call detect_and_change_otpdf(otpdf)

 inpname1 = TRIM(inpname)//'.t'
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
  write(6,'(A)') "ERROR in subroutine prt_mcpdft_molcas_inp: no 'SEWARD'&
                 & found in file "//TRIM(inpname)
  stop
 end if

 write(fid2,'(A)') 'Grid input'
 write(fid2,'(A)') ' grid=ultrafine'
 write(fid2,'(A)') 'End of grid input'
 if(DKH2) write(fid2,'(A)') ' Relativistic = R02O'

 ! I used to use &RASSCF, DMRG, RGinput for DMRGCI calculations. I found it
 ! seems that RGinput in &RASSCF had many disadvantages and often converged
 ! to higher energy local minimum. So I changed in into &DMRGSCF, and it
 ! seems to get more reasonable results

 dmrg = (dmrgci .or. dmrgscf)
 call prt_molcas_cas_para(fid2, dmrg, .false., .false., .true., inpname)

 write(fid2,'(/,A)') "&MCPDFT"
 write(fid2,'(A)') 'KSDFT='//TRIM(otpdf)

 if(mcpdft_force) then
  write(fid2,'(A)') 'Grad'
  write(fid2,'(/,A,/)') '&ALASKA'
 end if

 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_mcpdft_molcas_inp

! print MC-PDFT keywords into GAMESS .inp file
subroutine prt_mcpdft_gms_inp(inpname)
 use mol, only: charge, mult, ndb, nacte, nacto, npair, npair0
 use mr_keyword, only: mem, cart, otpdf, DKH2, hardwfn, crazywfn, CIonly, &
  mcpdft_force
 implicit none
 integer :: i, ncore, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 ncore = ndb + npair - npair0
 inpname1 = TRIM(inpname)//'.t'

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
subroutine detect_and_change_otpdf(otpdf)
 implicit none
 integer :: i, k, iv, fid, system
 character(len=3) :: sv
 character(len=9), intent(inout) :: otpdf
 character(len=10), parameter :: ftmp = 'molcas.ver'
 character(len=240) :: buf

 if(LEN_TRIM(otpdf)==0 .or. otpdf=='NONE') then
  write(6,'(/,A)') 'ERROR in subroutine detect_and_change_otpdf: OtPDF=None.'
  write(6,'(A)') 'Something must be wrong.'
  stop
 end if

 ! find the OpenMolcas version
 i = SYSTEM('pymolcas --banner >'//ftmp)
 open(newunit=fid,file=ftmp,status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  k = INDEX(buf,'version:')
  if(k > 0) exit
 end do ! for while

 close(fid,status='delete')
 if(i /= 0) return
 ! if something like 'version: v22.06' not found, just return

 i = INDEX(buf, '.')
 read(buf(k+8:i-1),*) sv     ! 'v22'
 read(sv(2:3),*,iostat=i) iv ! 22
 if(i /= 0) return
 ! if the user uses some strange version of OpenMolcas, just return

 call lower(otpdf)
 if(iv >= 22) then ! T:PBE, FT:PBE
  i = INDEX(otpdf, ':')
  if(i == 0) then
   select case(otpdf(1:1))
   case('t')
    otpdf = 't:'//TRIM(otpdf(2:))
   case('f')
    otpdf = 'ft:'//TRIM(otpdf(3:))
   case default
    write(6,'(A)') 'ERROR in subroutine detect_and_change_otpdf: wrong OtPDF='&
                   //TRIM(OtPDF)
    stop
   end select
  end if

 else ! tPBE, ftPBE
  i = INDEX(otpdf, ':')
  if(i > 0) then
   select case(otpdf(1:1))
   case('t')
    otpdf = 't'//TRIM(otpdf(i+1:))
   case('f')
    otpdf = 'ft'//TRIM(otpdf(i+1:))
   case default
    write(6,'(A)') 'ERROR in subroutine detect_and_change_otpdf: wrong OtPDF='&
                   //TRIM(OtPDF)
   end select
  end if
 end if
end subroutine detect_and_change_otpdf

