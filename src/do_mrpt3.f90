! written by jxzou at 20210111: 3rd order MRPT interfaces

subroutine do_mrpt3()
 use mr_keyword, only: dmrgci, dmrgscf, CIonly, caspt3, nevpt3, casnofch, &
  casscf_prog, casci_prog, bgchg, chgname, mem, nproc, molpro_path, bdf_path, &
  eist
 use mol, only: nacte, nacto, caspt2_e, nevpt2_e, caspt3_e, nevpt3_e, ptchg_e, &
  nuc_pt_e
 implicit none
 integer :: i, mem0, RENAME, system
 character(len=24) :: data_string
 character(len=240) :: pyname, mklname
 character(len=240) :: string, outname, inpname, inporb
 real(kind=8) :: ref_e, corr2_e, corr3_e
 logical :: alive(2)

 if(eist == 1) return ! excited state calculation
 alive = [caspt3, nevpt3]
 if(ALL(alive .eqv. .false.)) return

 write(6,'(//,A)') 'Enter subroutine do_mrpt3...'
 mem0 = CEILING(DBLE(mem*125)/DBLE(nproc))

 if((dmrgci .or. dmrgscf)) then
  write(6,'(/,A)') 'ERROR in subroutine do_mrpt3: CASPT3/NEVPT3 based on DMRG r&
                   &eference is not supported.'
  stop
 end if

 if(.not. CIonly) then
  if(casscf_prog == 'orca') then
   write(6,'(A)') 'Warning: ORCA is used as the CASSCF solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(6,'(A)') 'This will affect the PT3 energy up to 10^-5 a.u.&
                    & Such small error is usually not important.'
   write(6,'(A)') 'If you care about the accuracy, please use another CASSCF solver.'
  end if

  if(caspt3) then
   string = 'CASPT3 based on CASSCF orbitals.'
  else
   string = 'NEVPT3 based on CASSCF orbitals.'
  end if

 else ! CIonly = .True.
  if(casci_prog == 'orca') then
   write(6,'(A)') 'Warning: ORCA is used as the CASCI solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(6,'(A)') 'This will affect the PT3 energy up to 10^-5 a.u.&
                    & Such small error is usually not important.'
   write(6,'(A)') 'If you care about the accuracy, please use another CASCI solver.'
  end if

  string = 'CASPT3 based on CASCI orbitals.'
  write(6,'(A)') 'Warning: the CASSCF orbital optimization is strongly&
                     & recommended'
  write(6,'(A)') 'to be performed before PT3, unless it is too time-consuming.'
 end if

 write(6,'(A)') TRIM(string)
 write(6,'(A)',advance='no') 'Frozen_core = F, '

 if(caspt3) then ! CASPT3
  write(6,'(A,2(I0,A))') 'CASPT3(',nacte,'e,',nacto,'o) using program molpro'
  call check_exe_exist(molpro_path)

  i = SYSTEM('fch2com '//TRIM(casnofch))
  i = INDEX(casnofch, '.fch', back=.true.)
  pyname = casnofch(1:i-1)//'.com'
  string = casnofch
  call convert2molpro_fname(string, '.a')
  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i-1)//'_CASPT3.com'
  outname = casnofch(1:i-1)//'_CASPT3.out'
  inporb = inpname
  call convert2molpro_fname(inporb, '.a')
  i = RENAME(TRIM(pyname), TRIM(inpname))
  i = RENAME(TRIM(string), TRIM(inporb))

  call prt_mrpt_molpro_inp(inpname, 3) ! 1/2/3 for NEVPT2/CASPT2/CASPT3
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  write(string,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -t 1 -m ',mem0, &
                           'm '//TRIM(inpname)
  i = SYSTEM(TRIM(string))

 else ! FIC-NEVPT3
  write(6,'(A,2(I0,A))') 'NEVPT3(',nacte,'e,',nacto,'o) using program BDF'
  call check_exe_exist(bdf_path)

  i = SYSTEM('fch2bdf '//TRIM(casnofch)//' -no')
  i = INDEX(casnofch, '.fch', back=.true.)
  inporb = casnofch(1:i-1)//'_bdf.inporb'
  string = casnofch(1:i-1)//'_bdf.inp'
  i = INDEX(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'NEVPT3.inporb'
  inpname = casnofch(1:i)//'NEVPT3.inp'
  outname = casnofch(1:i)//'NEVPT3.out'
  i = RENAME(TRIM(inporb), TRIM(mklname))
  i = RENAME(TRIM(string), TRIM(inpname))

  call prt_mrpt_bdf_inp(inpname, 3) ! 1/2/3 for SDSPT2/NEVPT2/NEVPT3
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  i = INDEX(inpname, '.inp', back=.true.)
  string = inpname(1:i-1)
  i = SYSTEM(TRIM(bdf_path)//' '//TRIM(string))
 end if

 if(i /= 0) then
  if(caspt3) then
   write(6,'(/,A)') 'ERROR in subroutine do_mrpt3: CASPT3 computation failed.'
  else
   write(6,'(/,A)') 'ERROR in subroutine do_mrpt3: NEVPT3 computation failed.'
  end if
  write(6,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
  stop
 end if

 if(caspt3) then ! read CASPT3 energy
  call read_caspt3_energy_from_molpro_out(outname, ref_e, corr2_e, corr3_e)
  ref_e = ref_e + ptchg_e
 else            ! read NEVPT3 energy
  call read_nevpt3_energy_from_bdf_out(outname, ref_e, corr2_e, corr3_e)
  ! here the parameter davidson_e is useless
  ref_e = ref_e + ptchg_e + nuc_pt_e
 end if

 write(6,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e, 'a.u.'

 if(caspt3) then ! CASPT2
  write(6,'(A)')             'IP-EA shift  =         0.25 (default)'
  caspt2_e = ref_e + corr2_e
  caspt3_e = ref_e + corr3_e
  write(6,'(A,F18.8,1X,A4)') 'E(corr2)     = ', corr2_e, 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(CASPT2)    = ', caspt2_e,'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(corr3)     = ', corr3_e, 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(CASPT3)    = ', caspt3_e,'a.u.'
 else            ! FIC-NEVPT3
  nevpt2_e = ref_e + corr2_e
  nevpt3_e = ref_e + corr3_e
  write(6,'(A,F18.8,1X,A4)') 'E(corr2)     = ', corr2_e, 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(FIC-NEVPT2)= ', nevpt2_e,'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(corr3)     = ', corr3_e, 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(FIC-NEVPT3)= ', nevpt3_e,'a.u.'
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mrpt3 at '//TRIM(data_string)
end subroutine do_mrpt3

! read CASPT3 electronic energy from Molpro .out file
subroutine read_caspt3_energy_from_molpro_out(outname, ref_e, corr2_e, corr3_e)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: ref_e, corr2_e, corr3_e
 character(len=19), parameter :: key(6) = ['!MCSCF STATE  1.1 E', &
 '!MCSCF STATE 1.1 En', '!RSPT2 STATE  1.1 E', '!RSPT2 STATE 1.1 En', &
 '!RSPT3 STATE  1.1 E', '!RSPT3 STATE 1.1 En']
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0
 corr2_e = 0d0
 corr3_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:20)==key(1) .or. buf(2:20)==key(2)) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_caspt3_energy_from_molpro_out:'
  write(6,'(A)') 'CASSCF energy not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:20)==key(3) .or. buf(2:20)==key(4)) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_caspt3_energy_from_molpro_out:'
  write(6,'(A)') 'CASPT2 energy not found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) corr2_e
 corr2_e = corr2_e - ref_e

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:20)==key(5) .or. buf(2:20)==key(6)) exit
 end do ! for while
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_caspt3_energy_from_molpro_out:'
  write(6,'(A)') 'CASPT3 energy not found in file '//TRIM(outname)
  stop
 end if

 i = INDEX(buf,'ergy')
 read(buf(i+4:),*) corr3_e
 corr3_e = corr3_e - ref_e
end subroutine read_caspt3_energy_from_molpro_out

! read NEVPT3 electronic energy from BDF .out file
subroutine read_nevpt3_energy_from_bdf_out(outname, ref_e, corr2_e, corr3_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: ref_e, corr2_e, corr3_e
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ref_e = 0d0
 corr2_e = 0d0
 corr3_e = 0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:6) == 'NROOT') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mrpt_nevpt3_from_bdf_out: no&
                    & 'NROOT' found in file "//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,*) i, ref_e, corr2_e, rtmp, corr3_e
 close(fid)
 corr2_e = corr2_e - ref_e
 corr3_e = corr3_e - ref_e
end subroutine read_nevpt3_energy_from_bdf_out

