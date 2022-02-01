! written by jxzou at 20210727

! do MRCC based on CASCI/CASSCF, npair<=8
subroutine do_mrcc()
 use print_id, only: iout
 use mr_keyword, only: bgchg, mrcc, CIonly, mrcc_type, mrcc_prog, casnofch, &
  orca_path, chgname
 use mol, only: nevpt2_e, mrcc_e
 use util_wrapper, only: mkl2gbw
 implicit none
 integer :: i, system, RENAME
 real(kind=8) :: ref_e, corr_e(2)
 character(len=12), parameter :: method(8) = ['FIC-MRCCSD  ','Mk-MRCCSD   ',&
  'Mk-MRCCSD(T)','BW-MRCCSD   ','BW-MRCCSD(T)','BCCC2b      ','BCCC3b      ',&
  'BCCC4b      ']
 character(len=24) :: data_string
 character(len=240) :: string, chkname, inpname, mklname, outname

 if(.not. mrcc) return
 write(iout,'(//,A)') 'Enter subroutine do_mrcc...'
 write(iout,'(A)',advance='no') TRIM(method(mrcc_type))//' based on CAS'
 if(CIonly) then
  write(iout,'(A)') 'CI orbitals.'
 else
  write(iout,'(A)') 'SCF orbitals.'
 end if
 write(iout,'(A)') 'Frozen_core = F. Using program '//TRIM(mrcc_prog)

 select case(TRIM(mrcc_prog))
 case('orca')
  write(iout,'(A)') 'Note: 1) this is actually an approximate FIC-MRCCSD metho&
                   &d since the'
  write(iout,'(A)') '         H_bar operator is truncated after the quadratic terms.'
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
 case('nwchem')
  ! call prt_mrcc_nwchem_inp(inpname)
 case default
  write(iout,'(/,A)') 'ERROR in subroutine do_mrcc: MRCC_prog='//TRIM(mrcc_prog)
  write(iout,'(A)') 'Currently only MRCC_prog=ORCA/NWChem is supported.'
  stop
 end select

 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine do_mrcc: MRCC job failed.'
  write(iout,'(A)') 'Please find the error information in file '//TRIM(outname)
  stop
 end if

 ! read MRCC energy from ORCA output file
 call read_mrcc_energy_from_output(mrcc_prog, mrcc_type, outname, ref_e, corr_e)
 mrcc_e = ref_e + corr_e(2)

 select case(TRIM(mrcc_prog))
 case('orca')
  select case(mrcc_type)
  case(1)
   nevpt2_e = ref_e + corr_e(1)
   write(iout,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e,    'a.u.'
   write(iout,'(A,F18.8,1X,A4)')  'E(corr_PT2)  = ', corr_e(1),'a.u.'
   write(iout,'(A,F18.8,1X,A4)')  'E(FIC-NEVPT2)= ', nevpt2_e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)')  'E(corr_CCSD) = ', corr_e(2),'a.u.'
   write(iout,'(A,F18.8,1X,A4)')  'E(FIC-MRCCSD)= ', mrcc_e,   'a.u.'
  case(2) ! Mk-MRCCSD
   write(iout,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e,    'a.u.'
   write(iout,'(A,F18.8,1X,A4)')  'E(corr_CCSD) = ', corr_e(2),'a.u.'
   write(iout,'(A,F18.8,1X,A4)')  'E(Mk-MRCCSD) = ', mrcc_e,   'a.u.'
  case(4) ! BW-MRCCSD
   write(iout,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e,    'a.u.'
   write(iout,'(A,F18.8,1X,A4)')  'E(corr_CCSD) = ', corr_e(2),'a.u.'
   write(iout,'(A,F18.8,1X,A4)')  'E(BW-MRCCSD) = ', mrcc_e,   'a.u.'
  case default
   write(iout,'(A,I0)') 'ERROR in subroutine do_mrcc: invalid mrcc_type=',mrcc_type
   stop
  end select
 case('nwchem')
  select case(mrcc_type)
  case(2) ! Mk-MRCCSD
  case(3) ! Mk-MRCCSD(T)
  case(4) ! BW-MRCCSD
  case(5) ! BW-MRCCSD(T)
  case default
   write(iout,'(A,I0)') 'ERROR in subroutine do_mrcc: invalid mrcc_type=',mrcc_type
   stop
  end select
 end select

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_mrcc at '//TRIM(data_string)
 return
end subroutine do_mrcc

! print FIC-MRCCSD keywords into ORCA .inp file
subroutine prt_mrcc_orca_inp(inpname1)
 use mol, only: ndb, nopen, nacta, nactb, npair, npair0, mult
 use mr_keyword, only: mem, nproc, DKH2, mrcc_type, RI
 implicit none
 integer :: i, fid1, fid2, RENAME
 character(len=240), intent(in) :: inpname1
 character(len=240) :: buf, inpname2

 inpname2 = TRIM(inpname1)//'.tmp'
 open(newunit=fid1,file=TRIM(inpname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname2),status='replace')
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid2,'(A,I0,A)') '%maxcore ', CEILING(1d3*DBLE(mem)/DBLE(nproc))
 if(mrcc_type == 1) then
  write(fid2,'(A)') '! NoIter'
 else
  write(fid2,'(A)') '! NoIter UHF CCSD'
 end if

 if(DKH2) then
  write(fid2,'(A)') '%rel'
  write(fid2,'(A)') ' method DKH'
  write(fid2,'(A)') ' order 2'
  write(fid2,'(A)') 'end'
 end if

 select case(mrcc_type)
 case(1) ! FIC-MRCCSD
  write(fid2,'(A)') '%autoci'
  write(fid2,'(A)') ' CItype FICMRCC'
  write(fid2,'(A,I0)') ' nel ',  nacta+nactb
  write(fid2,'(A,I0)') ' norb ', 2*npair0+nopen
  write(fid2,'(A,I0)') ' mult ', mult
  write(fid2,'(A)') ' nroots 1'
  write(fid2,'(A)') ' MaxIter 100'
 case(2,4) ! Mk-MRCCSD, BWMRCCSD
  write(fid2,'(A)') '%mdci'
  write(fid2,'(A)') ' STol 1e-6'
  write(fid2,'(A)') ' mrcc on'
  write(fid2,'(A)',advance='no') ' mrcctype '
  if(mrcc_type == 2) then
   write(fid2,'(A)') 'mkcc'
  else
   write(fid2,'(A)') 'bwcc'
  end if
  write(fid2,'(A,I0)') ' n_docc ', ndb+npair-npair0
  write(fid2,'(A)') " refs ""2200,2020,2002,0220,0202,0022"""
 case default
  write(fid2,'(A,I0)') 'ERROR in subroutine prt_mrcc_orca_inp: invalid mrcc_typ&
                       &e=',mrcc_type
  stop
 end select
 write(fid2,'(A)') 'end'
 write(fid2,'(A)') '%method'
 write(fid2,'(A)') ' FrozenCore FC_NONE'
 if(RI .and. mrcc_type>1) write(fid2,'(A)') ' RIJKSinglesFock 1'
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

subroutine read_mrcc_energy_from_output(mrcc_prog, mrcc_type, outname, ref_e, &
                                        corr_e)
 use print_id, only: iout
 implicit none
 integer :: i, fid
 integer, intent(in) :: mrcc_type ! 1/2/3 for FIC-/Mk-/BW-MRCCSD
 ! 4/5 for Mk-/BW-MRCCSD(T), 6~8 for GVB-BCCC2b/3b/4b
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: ref_e, corr_e(2)
 character(len=240) :: buf
 character(len=10), intent(in) :: mrcc_prog
 character(len=240), intent(in) :: outname

 ref_e = 0d0; corr_e = 0d0
 ! corr_e(1): FIC-NEVPT2 energy
 ! corr_e(2): FIC-/Mk-/BW- MRCCSD energy

 select case(mrcc_prog)
 case('orca')
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

 case('nwchem')

 case default
  write(iout,'(/,A)') 'ERROR in subroutine read_mrcc_energy_from_output: invalid&
                     & mrcc_prog='//TRIM(mrcc_prog)
  stop
 end select

 return
end subroutine read_mrcc_energy_from_output

