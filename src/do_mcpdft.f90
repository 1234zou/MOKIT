! written by jxzou at 20210120: subroutine do_mcpdft from automr.f90 is moved here

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
  write(fid2,'(A)') ' conv_thresh = 1E-8'
  write(fid2,'(A)') ' nsweeps = 20'
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
 return
end subroutine prt_mcpdft_gms_inp

! do MC-PDFT for npair<=7, or <=CAS(14,14)
subroutine do_mcpdft()
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, casci, casscf, dmrgci, dmrgscf, mcpdft, &
  mcpdft_prog, casnofch, molcas_path, gms_path, bgchg, chgname, check_gms_path,&
  gms_scr_path
 use mol, only: casci_e, casscf_e, ptchg_e, mcpdft_e
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, system, RENAME
 real(kind=8) :: ref_e, corr_e
 character(len=24) :: data_string
 character(len=240) :: fname(3), inpname, outname, cmofch
 logical :: dmrg

 if(.not. mcpdft) return
 write(iout,'(//,A)') 'Enter subroutine do_mcpdft...'

 dmrg = (dmrgci .or. dmrgscf)

 if(dmrg) then
  if(mcpdft_prog=='gamess') then
   write(iout,'(A)') 'ERROR in subroutine do_mcpdft: DMRG-PDFT has not been&
                    & implemented in GAMESS.'
   write(iout,'(A)') 'You can set MCPDFT_prog=OpenMolcas in mokit{}.'
   stop
  end if
  if(dmrgci) then
   write(iout,'(A)') 'DMRG-PDFT based on DMRG-CASCI orbitals.'
  else
   write(iout,'(A)') 'DMRG-PDFT based on optimized DMRG-CASSCF orbitals.'
  end if
  write(iout,'(A)') 'DMRG-PDFT using program openmolcas'
 else
  if(casci) then
   write(iout,'(A)') 'MC-PDFT based on CASCI orbitals.'
  else
   write(iout,'(A)') 'MC-PDFT based on optimized CASSCF orbitals.'
  end if
  write(iout,'(A)') 'MC-PDFT using program '//TRIM(mcpdft_prog)
 end if

 if(dmrgci .or. casci) then
  write(iout,'(A)') 'Warning: orbital optimization is strongly recommended to&
                   & be performed'
  write(iout,'(A)') 'before PDFT, unless it is too time-consuming.'
 end if

 select case(TRIM(mcpdft_prog))
 case('openmolcas')
  call check_exe_exist(molcas_path)
  ! For DMRG-PDFT, use CMOs rather than NOs
  if(dmrg) then
   i = index(casnofch, '_NO', back=.true.)
   cmofch = casnofch(1:i)//'CMO.fch'
   casnofch = cmofch
  end if

  i = system('fch2inporb '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  fname(1) = casnofch(1:i-1)//'.INPORB'
  fname(2) = casnofch(1:i-1)//'.input'
  if(dmrg) then
   i = index(casnofch, '_CMO', back=.true.)
  else
   i = index(casnofch, '_NO', back=.true.)
  end if
  fname(3) = casnofch(1:i)//'MCPDFT.INPORB'
  inpname  = casnofch(1:i)//'MCPDFT.input'
  outname  = casnofch(1:i)//'MCPDFT.out'
  i = RENAME(TRIM(fname(1)), TRIM(fname(3)))
  i = RENAME(TRIM(fname(2)), TRIM(inpname))
 
  call prt_mcpdft_molcas_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  i = system(TRIM(molcas_path)//' '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1")

 case('gamess')
  call check_gms_path()
  call fch2inp_wrap(casnofch, .false., 0, 0)
  i = index(casnofch, '.fch', back=.true.)
  fname(1) = casnofch(1:i-1)//'.inp'
  i = index(casnofch, '_NO', back=.true.)
  fname(2) = casnofch(1:i)//'MCPDFT.dat'
  inpname = casnofch(1:i)//'MCPDFT.inp'
  outname = casnofch(1:i)//'MCPDFT.gms'
  i = RENAME(TRIM(fname(1)), TRIM(inpname))
  fname(2) = TRIM(gms_scr_path)//'/'//TRIM(fname(2)) ! delete the possible .dat file
  call delete_file(fname(2))

  call prt_mcpdft_gms_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  write(fname(3),'(A,I0,A)') TRIM(gms_path)//' '//TRIM(inpname)//" 01 1 >"//&
                             TRIM(outname)//" 2>&1"
  i = system(TRIM(fname(3)))
  ! MC-PDFT in GAMESS cannot run in parallel currently, use 1 core

  ! move the .dat file into current directory
  i = system('mv '//TRIM(fname(2))//' .')
 end select

 ! read MC-PDFT/DMRG-PDFT energy from OpenMolcas output file
 call read_mcpdft_e_from_output(mcpdft_prog, outname, ref_e, corr_e)
 if(TRIM(mcpdft_prog) == 'openmolcas') ref_e = ref_e + ptchg_e
 mcpdft_e = ref_e + corr_e

 write(iout,'(/,A,F18.8,1X,A4)')'E(ref)      = ',    ref_e, 'a.u.'
 if(dmrg) then
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)     = ',   corr_e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(DMRG-PDFT)= ', mcpdft_e, 'a.u.'
 else
  write(iout,'(A,F18.8,1X,A4)') 'E(corr)     = ',   corr_e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(MC-PDFT)  = ', mcpdft_e, 'a.u.'
 end if

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_mcpdft at '//TRIM(data_string)
 return
end subroutine do_mcpdft

