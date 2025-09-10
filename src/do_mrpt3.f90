! written by jxzou at 20210111: 3rd order MRPT interfaces

subroutine do_mrpt3()
 use mr_keyword, only: dmrgci, dmrgscf, CIonly, caspt3, nevpt3, nevpt_prog, &
  casnofch, bgchg, chgname, mem, nproc, molpro_path, orca_path, bdf_path, eist
 use mol, only: nacte, nacto, caspt2_e, nevpt2_e, caspt3_e, nevpt3_e, ptchg_e, &
  nuc_pt_e
 use util_wrapper, only: fch2mkl_wrap, mkl2gbw
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=24) :: data_string
 character(len=240) :: pyname, mklname
 character(len=240) :: string, outname, inpname, inporb
 real(kind=8) :: ref_e, corr_e(3)
 logical :: alive(2)

 if(eist == 1) return ! excited state calculation
 alive = [caspt3, nevpt3]
 if(ALL(alive .eqv. .false.)) return
 write(6,'(//,A)') 'Enter subroutine do_mrpt3...'

 if((dmrgci .or. dmrgscf)) then
  write(6,'(/,A)') 'ERROR in subroutine do_mrpt3: CASPT3/NEVPT3 based on DMRG r&
                   &eference is not supported.'
  stop
 end if

 if(.not. CIonly) then
  if(caspt3) then
   string = 'CASPT3 based on CASSCF orbitals.'
  else
   string = 'NEVPT3 based on CASSCF orbitals.'
  end if
 else ! CIonly = .True.
  string = 'CASPT3 based on CASCI orbitals.'
  write(6,'(A)') 'Warning: the CASSCF orbital optimization is strongly recommen&
                 &ded to be'
  write(6,'(A)') 'performed before PT3, unless it is too time-consuming.'
 end if

 write(6,'(A)') TRIM(string)
 write(6,'(A)',advance='no') 'Frozen_core = F, '

 if(caspt3) then ! CASPT3
  write(6,'(A,2(I0,A))') 'CASPT3(',nacte,'e,',nacto,'o) using program molpro'
  call check_exe_exist(molpro_path)
  i = SYSTEM('fch2com '//TRIM(casnofch))
  call find_specified_suffix(casnofch, '.fch', i)
  pyname = casnofch(1:i-1)//'.com'
  string = casnofch
  call convert2molpro_fname(string, '.a')
  call find_specified_suffix(casnofch, '_NO', i)
  inpname = casnofch(1:i-1)//'_CASPT3.com'
  outname = casnofch(1:i-1)//'_CASPT3.out'
  inporb = inpname
  call convert2molpro_fname(inporb, '.a')
  i = RENAME(TRIM(pyname), TRIM(inpname))
  i = RENAME(TRIM(string), TRIM(inporb))
  call prt_mrpt_molpro_inp(inpname, 3) ! 1/2/3 for NEVPT2/CASPT2/CASPT3
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_molpro_job(inpname, mem, nproc)

 else ! FIC-NEVPT3
  write(6,'(2(A,I0),A)') 'NEVPT3(', nacte, 'e,', nacto, 'o) using program '//&
                         TRIM(nevpt_prog)
  select case(TRIM(nevpt_prog))
  case('orca')
   call check_exe_exist(orca_path)
   call fch2mkl_wrap(casnofch)
   call find_specified_suffix(casnofch, '.fch', i)
   inporb = casnofch(1:i-1)//'_o.mkl'
   string = casnofch(1:i-1)//'_o.inp'
   call find_specified_suffix(casnofch, '_NO', i)
   mklname = casnofch(1:i)//'NEVPT3.mkl'
   inpname = casnofch(1:i)//'NEVPT3.inp'
   outname = casnofch(1:i)//'NEVPT3.out'
   i = RENAME(TRIM(inporb), TRIM(mklname))
   i = RENAME(TRIM(string), TRIM(inpname))
   call prt_nevpt34_orca_inp(inpname, .false.)
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   ! if bgchg = .True., .inp and .mkl file will be updated
   call mkl2gbw(mklname)
   call delete_file(mklname)
   call submit_orca_job(orca_path, inpname, .true., .false., .false.)

  case('bdf')
   call check_exe_exist(bdf_path)
   i = SYSTEM('fch2bdf '//TRIM(casnofch)//' -no')
   call find_specified_suffix(casnofch, '.fch', i)
   inporb = casnofch(1:i-1)//'_bdf.inporb'
   string = casnofch(1:i-1)//'_bdf.inp'
   call find_specified_suffix(casnofch, '_NO', i)
   mklname = casnofch(1:i)//'NEVPT3.inporb'
   inpname = casnofch(1:i)//'NEVPT3.inp'
   outname = casnofch(1:i)//'NEVPT3.out'
   i = RENAME(TRIM(inporb), TRIM(mklname))
   i = RENAME(TRIM(string), TRIM(inpname))
   call prt_mrpt_bdf_inp(inpname, 3) ! 1/2/3 for SDSPT2/NEVPT2/NEVPT3
   if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
   call find_specified_suffix(inpname, '.inp', i)
   string = inpname(1:i-1)
   i = SYSTEM(TRIM(bdf_path)//' '//TRIM(string))
   if(i /= 0) then
    write(6,'(/,A)') 'ERROR in subroutine do_mrpt3: NEVPT3 computation failed i&
                     &n BDF.'
    write(6,'(A)') 'You can open file '//TRIM(outname)//' and check why.'
    stop
   end if
  case default
   write(6,'(/,A)') 'ERROR in subroutine do_mrpt3: NEVPT_prog cannot be recogn&
                    &ized. Allowed'
   write(6,'(A)') 'NEVPT_prog are ORCA/BDF. But currently NEVPT_prog='//&
                  TRIM(nevpt_prog)
   stop
  end select
 end if

 if(caspt3) then ! read CASPT3 energy
  call read_caspt3_energy_from_molpro_out(outname, ref_e, corr_e(1), corr_e(2))
  ref_e = ref_e + ptchg_e
 else            ! read NEVPT3 energy
  if(TRIM(nevpt_prog) == 'orca') then
   call read_nevpt34_e_from_orca_out(outname, .false., ref_e, corr_e)
  else
   call read_nevpt3_energy_from_bdf_out(outname, ref_e, corr_e(1), corr_e(2))
  end if
  ref_e = ref_e + ptchg_e + nuc_pt_e
 end if
 ! corr_e(3) is not used in this subroutine

 write(6,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e, 'a.u.'

 if(caspt3) then ! CASPT2
  write(6,'(A)')             'IP-EA shift  =         0.25 (default)'
  caspt2_e = ref_e + corr_e(1)
  caspt3_e = ref_e + corr_e(2)
  write(6,'(A,F18.8,1X,A4)') 'E(corr2)     = ', corr_e(1), 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(CASPT2)    = ', caspt2_e,'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(corr3)     = ', corr_e(2), 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(CASPT3)    = ', caspt3_e,'a.u.'
 else            ! FIC-NEVPT3
  nevpt2_e = ref_e + corr_e(1)
  nevpt3_e = ref_e + corr_e(2)
  write(6,'(A,F18.8,1X,A4)') 'E(corr2)     = ', corr_e(1), 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(FIC-NEVPT2)= ', nevpt2_e,'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(corr3)     = ', corr_e(2), 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(FIC-NEVPT3)= ', nevpt3_e,'a.u.'
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mrpt3 at '//TRIM(data_string)
end subroutine do_mrpt3

! print FIC-NEVPT3/FIC-NEVPT4(SD) keywords in to a given ORCA .inp file
subroutine prt_nevpt34_orca_inp(inpname, nevpt4)
 use mol, only: mult, nacte, nacto
 use mr_keyword, only: mem, nproc, xmult, iroot, RI, RIJK_bas, RIC_bas, F12, &
  F12_cabs, DLPNO, hardwfn, crazywfn
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: nevpt4

 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 read(fid1,'(A)') buf   ! skip nproc
 read(fid1,'(A)') buf   ! skip memory
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid2,'(A,I0)') '%maxcore ', FLOOR(1d3*DBLE(mem)/DBLE(nproc))

 read(fid1,'(A)') buf   ! skip '!' line
 if(nevpt4) then
  write(fid2,'(A)',advance='no') '! FIC-NEVPT4(SD)'
 else
  write(fid2,'(A)',advance='no') '! FIC-NEVPT3'
 end if
 ! TODO: use only RIJCOSX in ORCA
 if(RI) then
  write(fid2,'(A)',advance='no') ' RIJK '//TRIM(RIJK_bas)//' '//TRIM(RIC_bas)
 end if
 if(DLPNO) write(fid2,'(A)',advance='no') ' TightPNO'
 if(F12) write(fid2,'(A)',advance='no') ' '//TRIM(F12_cabs)
 write(fid2,'(A)') ' NoIter'

 write(fid2,'(A)') '%casscf'
 write(fid2,'(1X,A,I0)') 'nel ', nacte
 write(fid2,'(1X,A,I0)') 'norb ', nacto
 if(iroot > 0) call prt_orca_ss_cas_weight(fid2, mult, xmult, iroot)
 call prt_hard_or_crazy_casci_orca(1, fid2, hardwfn, crazywfn)
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
end subroutine prt_nevpt34_orca_inp

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

 ref_e = 0d0; corr2_e = 0d0; corr3_e = 0d0
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

 ref_e = 0d0; corr2_e = 0d0; corr3_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:6) == 'NROOT') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_mrpt_nevpt3_from_bdf_out: no 'NROO&
                   &T' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,*) i, ref_e, corr2_e, rtmp, corr3_e
 close(fid)
 corr2_e = corr2_e - ref_e
 corr3_e = corr3_e - ref_e
end subroutine read_nevpt3_energy_from_bdf_out

! read FIC-NEVPT3/FIC-NEVPT4(SD) electronic energy from ORCA output file
subroutine read_nevpt34_e_from_orca_out(outname, pt4, ref_e, corr_e)
 implicit none
 integer :: i, fid
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: ref_e, corr_e(3)
 real(kind=8), parameter :: thres = 1d-5
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(in) :: pt4

 ref_e = 0d0; corr_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:5) == 'E(0)') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nevpt34_e_from_orca_out: no 'E(0)'&
                   & found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if
 call get_dpv_after_flag(buf, '=', .true., ref_e) ! E(0)

 read(fid,'(A)') buf
 call get_dpv_after_flag(buf, '=', .true., rtmp) ! E(1)
 if(DABS(rtmp) > thres) then
  write(6,'(/,A)') 'ERROR in subroutine read_nevpt34_e_from_orca_out: E(1) is n&
                   &ot zero.'
  write(6,'(A)') 'Unexpected case.'
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 call get_dpv_after_flag(buf, '=', .true., corr_e(1)) ! E(2)

 read(fid,'(A)') buf
 close(fid)
 call get_dpv_after_flag(buf, '=', .true., corr_e(2)) ! E(3)
 corr_e(2) = corr_e(2) + corr_e(1)

 if(pt4) then
  call get_dpv_after_flag(buf, '=', .true., corr_e(3)) ! E(3)
  corr_e(3) = corr_e(3) + corr_e(2)
 end if
end subroutine read_nevpt34_e_from_orca_out

! Temporarily put in this file. Maybe moved to do_mrpt4.f90 in the future.
subroutine do_mrpt4()
 use mr_keyword, only: dmrgci, dmrgscf, CIonly, nevpt4, nevpt_prog, &
  casnofch, bgchg, chgname, orca_path, eist
 use mol, only: nacte, nacto, nevpt2_e, nevpt3_e, nevpt4_e, ptchg_e, &
  nuc_pt_e
 use util_wrapper, only: fch2mkl_wrap, mkl2gbw
 implicit none
 integer :: i, SYSTEM, RENAME
 character(len=24) :: data_string
 character(len=240) :: string, outname, inpname, inporb, mklname
 real(kind=8) :: ref_e, corr_e(3)

 if(eist == 1) return ! excited state calculation
 if(.not. nevpt4) return
 write(6,'(//,A)') 'Enter subroutine do_mrpt4...'

 if((dmrgci .or. dmrgscf)) then
  write(6,'(/,A)') 'ERROR in subroutine do_mrpt4: NEVPT4 based on DMRG referenc&
                   &e is not supported.'
  stop
 end if

 if(.not. CIonly) then
  string = 'NEVPT4(SD) based on CASSCF orbitals.'
 else ! CIonly = .True.
  string = 'NEVPT4(SD) based on CASCI orbitals.'
  write(6,'(A)') 'Warning: the CASSCF orbital optimization is strongly recommen&
                 &ded to be'
  write(6,'(A)') 'performed before PT4, unless it is too time-consuming.'
 end if

 write(6,'(A)') TRIM(string)
 write(6,'(A)',advance='no') 'Frozen_core = F, '

 write(6,'(2(A,I0),A)') 'FIC-NEVPT4(SD) (', nacte, 'e,', nacto, &
                        'o) using program '//TRIM(nevpt_prog)
 call check_exe_exist(orca_path)
 call fch2mkl_wrap(casnofch)
 call find_specified_suffix(casnofch, '.fch', i)
 inporb = casnofch(1:i-1)//'_o.mkl'
 string = casnofch(1:i-1)//'_o.inp'
 call find_specified_suffix(casnofch, '_NO', i)
 mklname = casnofch(1:i)//'NEVPT4.mkl'
 inpname = casnofch(1:i)//'NEVPT4.inp'
 outname = casnofch(1:i)//'NEVPT4.out'
 i = RENAME(TRIM(inporb), TRIM(mklname))
 i = RENAME(TRIM(string), TRIM(inpname))
 call prt_nevpt34_orca_inp(inpname, .true.)
 if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
 ! if bgchg = .True., .inp and .mkl file will be updated
 call mkl2gbw(mklname)
 call delete_file(mklname)
 call submit_orca_job(orca_path, inpname, .true., .false., .false.)

 call read_nevpt34_e_from_orca_out(outname, .true., ref_e, corr_e)
 ref_e = ref_e + ptchg_e + nuc_pt_e

 write(6,'(/,A,F18.8,1X,A4)') 'E(ref)           = ', ref_e, 'a.u.'
 nevpt2_e = ref_e + corr_e(1)
 nevpt3_e = ref_e + corr_e(2)
 nevpt4_e = ref_e + corr_e(3)
 write(6,'(A,F18.8,1X,A4)') 'E(corr2)         = ',corr_e(1), 'a.u.'
 write(6,'(A,F18.8,1X,A4)') 'E(FIC-NEVPT2)    = ', nevpt2_e, 'a.u.'
 write(6,'(A,F18.8,1X,A4)') 'E(corr3)         = ',corr_e(2), 'a.u.'
 write(6,'(A,F18.8,1X,A4)') 'E(FIC-NEVPT3)    = ', nevpt3_e, 'a.u.'
 write(6,'(A,F18.8,1X,A4)') 'E(corr4)         = ',corr_e(3), 'a.u.'
 write(6,'(A,F18.8,1X,A4)') 'E(FIC-NEVPT4(SD))= ', nevpt4_e, 'a.u.'

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mrpt4 at '//TRIM(data_string)
end subroutine do_mrpt4

