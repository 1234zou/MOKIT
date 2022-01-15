! written by jxzou at 20210727

! do MRCC based on CASCI/CASSCF, npair<=8
subroutine do_mrcc()
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, bgchg, ficmrcc, bccc2b, casci, casscf, &
  CIonly, mrcc_prog, casnofch, orca_path, chgname
 use mol, only: nevpt2_e, mrcc_e
 use util_wrapper, only: mkl2gbw
 implicit none
 integer :: i, system, RENAME
 real(kind=8) :: ref_e, corr_e(2)
 character(len=24) :: data_string
 character(len=240) :: string, chkname, inpname, mklname, outname

 if(.not. (ficmrcc .or. bccc2b)) return
 write(iout,'(//,A)') 'Enter subroutine do_mrcc...'

 select case(TRIM(mrcc_prog))
 case('orca')
  if(CIonly) then
   write(iout,'(A)') 'FIC-MRCC based on CASCI orbitals.'
  else
   write(iout,'(A)') 'FIC-MRCC based on optimized CASSCF orbitals.'
  end if
  write(iout,'(A)') 'Frozen_core = F, FIC-MRCC using program orca'
  write(iout,'(A)') 'Note: 1) this is actually an approximate FIC-MRCC method since&
                   & the H_bar'
  write(iout,'(A)') '         operator is truncated after the quadratic terms.'
  write(iout,'(A)') '      2) this calculation will output the FIC-NEVPT2 energy&
                   & as a byproduct.'
  write(iout,'(A)') '      3) FIC-MRCC is supported since ORCA 5.0. Do not use&
                   & any older version.'
  call check_exe_exist(orca_path)
  i = system('fch2mkl '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  chkname = casnofch(1:i-1)//'_o.mkl'
  string  = casnofch(1:i-1)//'_o.inp'
  i = index(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'MRCC.mkl'
  inpname = casnofch(1:i)//'MRCC.inp'
  outname = casnofch(1:i)//'MRCC.out'
  i = RENAME(TRIM(chkname), TRIM(mklname))
  i = RENAME(TRIM(string), TRIM(inpname))
  chkname = ' '
  call prt_mrcc_orca_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(mklname)
  string = TRIM(inpname)//' >'//TRIM(outname)//" 2>&1"
  write(iout,'(A)') '$ORCA '//TRIM(string)

  i = system(TRIM(orca_path)//' '//TRIM(string))
  if(i /= 0) then
   write(iout,'(A)') 'ERROR in subroutine do_mrcc: FIC-MRCC job failed.'
   write(iout,'(A)') 'Please find the error information in file '//TRIM(outname)
   stop
  end if

 case('gvb_bccc2b')
  if(CIonly) then
   write(iout,'(A)') 'ERROR in subroutine do_mrcc: CIonly=.True. found. It&
                    & should not be activated.'
   write(iout,'(A)') 'Because this option has nothing to do with GVB-BCCC2b.'
   stop
  end if

 case default
  write(iout,'(/,A)') 'ERROR in subroutine do_mrcc: MRCC_prog='//TRIM(mrcc_prog)
  write(iout,'(A)') 'Currently only MRCC_prog=ORCA or GVB_BCCC2b is supported.'
  stop
 end select

 ! read MRCC energy from ORCA output file
 call read_mrcc_energy_from_output(mrcc_prog, outname, ref_e, corr_e)

 select case(TRIM(mrcc_prog))
 case('orca')
  nevpt2_e = ref_e + corr_e(1)
  mrcc_e = ref_e + corr_e(2)
  write(iout,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e,    'a.u.'
  write(iout,'(A,F18.8,1X,A4)')  'E(corr_PT2)  = ', corr_e(1),'a.u.'
  write(iout,'(A,F18.8,1X,A4)')  'E(FIC-NEVPT2)= ', nevpt2_e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)')  'E(corr_CC)   = ', corr_e(2),'a.u.'
  write(iout,'(A,F18.8,1X,A4)')  'E(FIC-MRCC)  = ', mrcc_e,   'a.u.'
 case default
  ! GVB-BCCC2b
 end select

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_mrcc at '//TRIM(data_string)
 return
end subroutine do_mrcc

! print FIC-MRCC keywords into ORCA .inp file
subroutine prt_mrcc_orca_inp(inpname1)
 use print_id, only: iout
 use mol, only: nopen, nacta, nactb, npair0, mult
 use mr_keyword, only: mem, nproc, DKH2
 implicit none
 integer :: i, fid1, fid2
 integer :: RENAME
 character(len=240), intent(in) :: inpname1
 character(len=240) :: buf, inpname2

 inpname2 = TRIM(inpname1)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname2),status='replace')
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid2,'(A,I0,A)') '%maxcore ', CEILING(1d3*DBLE(mem)/DBLE(nproc))
 write(fid2,'(A)') '! NoIter'

 if(DKH2) then
  write(fid2,'(A)') '%rel'
  write(fid2,'(A)') ' method DKH'
  write(fid2,'(A)') ' order 2'
  write(fid2,'(A)') 'end'
 end if

 write(fid2,'(A)') '%autoci'
 write(fid2,'(A)') ' CItype FICMRCC'
 write(fid2,'(A,I0)') ' nel ',  nacta+nactb
 write(fid2,'(A,I0)') ' norb ', 2*npair0+nopen
 write(fid2,'(A,I0)') ' mult ', mult
 write(fid2,'(A)') ' nroots 1'
 write(fid2,'(A)') ' MaxIter 100'
 write(fid2,'(A)') 'end'
 write(fid2,'(A)') '%method'
 write(fid2,'(A)') ' FrozenCore FC_NONE'
 write(fid2,'(A)') 'end'
 write(fid2,'(A)') '%coords'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:6) == '%coord') exit
 end do ! for while

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname2), TRIM(inpname1))
 return
end subroutine prt_mrcc_orca_inp

subroutine read_mrcc_energy_from_output(mrcc_prog, outname, ref_e, corr_e)
 use print_id, only: iout
 implicit none
 integer :: i, fid
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: ref_e, corr_e(2)
 character(len=240) :: buf
 character(len=10), intent(in) :: mrcc_prog
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr_e = 0d0

 select case(mrcc_prog)
 case('orca')
  ! corr_e(1): FIC-NEVPT2 energy
  ! corr_e(2): FIC-MRCC energy
  open(newunit=fid,file=TRIM(outname),status='old',position='append')
  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(1:17) == ' Summary of multi') exit
  end do ! for while

  if(i /= 0) then
   write(iout,'(/,A)') "ERROR in subroutine read_mrcc_energy_from_output: ' Summary&
                      & of multi' not found"
   write(iout,'(A)') 'in file '//TRIM(outname)
   close(fid)
   stop
  end if

  read(fid,'(A)') buf
  read(fid,'(A)') buf
  read(fid,*) i,i,rtmp,rtmp,ref_e,corr_e(2)

  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(1:5) == 'MULT=') exit
  end do ! for while
  close(fid)

  if(i /= 0) then
   write(iout,'(/,A)') "ERROR in subroutine read_mrcc_energy_from_output: 'MULT='&
                      & not found in file "//TRIM(outname)
   stop
  end if

  i = index(buf,'EC=')
  read(buf(i+3:),*) corr_e(1)

 case('gvb_bccc2b') ! GVB-BCCC2b

 case default
  write(iout,'(/,A)') 'ERROR in subroutine read_mrcc_energy_from_output: invalid&
                     & mrcc_prog='//TRIM(mrcc_prog)
  stop
 end select

 return
end subroutine read_mrcc_energy_from_output

