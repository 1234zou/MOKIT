! written by jxzou at 20210727

! do MRCC based on CASCI/CASSCF, npair<=8
subroutine do_mrcc()
 use mr_keyword, only: bgchg, mrcc, CIonly, mrcc_type, mrcc_prog, casnofch, &
  orca_path, chgname, eist
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

 if(eist == 1) return ! excited state calculation
 if(.not. mrcc) return
 if(mrcc_type > 5) then
  call do_gvb_bccc()
  return
 end if

 write(6,'(//,A)') 'Enter subroutine do_mrcc...'
 write(6,'(A)',advance='no') TRIM(method(mrcc_type))//' based on CAS'
 if(CIonly) then
  write(6,'(A)') 'CI orbitals.'
 else
  write(6,'(A)') 'SCF orbitals.'
 end if
 write(6,'(A)') 'Frozen_core = F. Using program '//TRIM(mrcc_prog)

 select case(TRIM(mrcc_prog))
 case('orca')
  write(6,'(A)') 'Note: 1) this is actually an approximate FIC-MRCCSD method &
                 &since the'
  write(6,'(A)') '         H_bar operator is truncated after the quadratic terms.'
  write(6,'(A)') '      2) this calculation will output the FIC-NEVPT2 energy&
                 & as a byproduct.'
  write(6,'(A)') '      3) FIC-MRCC is supported since ORCA 5.0. Do not use&
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
  call submit_orca_job(orca_path, inpname)
 case('nwchem')
  ! call prt_mrcc_nwchem_inp(inpname)
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_mrcc: MRCC_prog='//TRIM(mrcc_prog)
  write(6,'(A)') 'Currently only MRCC_prog=ORCA/NWChem is supported.'
  stop
 end select

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_mrcc: MRCC job failed.'
  write(6,'(A)') 'Please find the error information in file '//TRIM(outname)
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
   write(6,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e,    'a.u.'
   write(6,'(A,F18.8,1X,A4)')  'E(corr_PT2)  = ', corr_e(1),'a.u.'
   write(6,'(A,F18.8,1X,A4)')  'E(FIC-NEVPT2)= ', nevpt2_e, 'a.u.'
   write(6,'(A,F18.8,1X,A4)')  'E(corr_CCSD) = ', corr_e(2),'a.u.'
   write(6,'(A,F18.8,1X,A4)')  'E(FIC-MRCCSD)= ', mrcc_e,   'a.u.'
  case(2) ! Mk-MRCCSD
   write(6,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e,    'a.u.'
   write(6,'(A,F18.8,1X,A4)')  'E(corr_CCSD) = ', corr_e(2),'a.u.'
   write(6,'(A,F18.8,1X,A4)')  'E(Mk-MRCCSD) = ', mrcc_e,   'a.u.'
  case(4) ! BW-MRCCSD
   write(6,'(/,A,F18.8,1X,A4)')'E(ref)       = ', ref_e,    'a.u.'
   write(6,'(A,F18.8,1X,A4)')  'E(corr_CCSD) = ', corr_e(2),'a.u.'
   write(6,'(A,F18.8,1X,A4)')  'E(BW-MRCCSD) = ', mrcc_e,   'a.u.'
  case default
   write(6,'(A,I0)') 'ERROR in subroutine do_mrcc: invalid mrcc_type=',mrcc_type
   stop
  end select
 case('nwchem')
  select case(mrcc_type)
  case(2) ! Mk-MRCCSD
  case(3) ! Mk-MRCCSD(T)
  case(4) ! BW-MRCCSD
  case(5) ! BW-MRCCSD(T)
  case default
   write(6,'(A,I0)') 'ERROR in subroutine do_mrcc: invalid mrcc_type=',mrcc_type
   stop
  end select
 end select

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mrcc at '//TRIM(data_string)
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
end subroutine prt_mrcc_orca_inp

subroutine read_mrcc_energy_from_output(mrcc_prog, mrcc_type, outname, ref_e, &
                                        corr_e)
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
   write(6,'(/,A)') "ERROR in subroutine read_mrcc_energy_from_output: ' Summary&
                   & of multi' not found"
   write(6,'(A)') 'in file '//TRIM(outname)
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
   write(6,'(/,A)') "ERROR in subroutine read_mrcc_energy_from_output: 'MULT='&
                   & not found in file "//TRIM(outname)
   stop
  end if

  i = index(buf,'EC=')
  read(buf(i+3:),*) corr_e(1)

 case('gvb_bcci2b','gvb_bccc2b','gvb_bccc3b')
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:5) == 'GVB e') exit
  end do ! for while

  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_mrcc_energy_from_output: no &
                   &'GVB e' found in"
   write(6,'(A)') 'file '//TRIM(outname)
   stop
  end if

  i = index(buf, '=')
  read(buf(i+1:),*) ref_e

  ! find if any 'ERROR: amplitude iterations fail'
  open(newunit=fid,file=TRIM(outname),status='old',position='append')
  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:32) == 'ERROR: amplitude iterations fail') then
    close(fid)
    write(6,'(A)') 'ERROR in subroutine read_mrcc_energy_from_output: amplitud&
                   &e iterations not converge.'
    stop
   end if
   if(buf(1:3) == 'ITN') exit
  end do ! for while

  close(fid)
  open(newunit=fid,file=TRIM(outname),status='old',position='append')

  do while(.true.)
   BACKSPACE(fid)
   BACKSPACE(fid)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(index(buf,'=') > 0) exit
  end do ! for while

  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  i = index(buf, '=')
  read(buf(i+1:),*) corr_e(1) ! GVB-BCCCnb correlation energy

  read(fid,'(A)') buf
  i = index(buf, '=')
  read(buf(i+1:),*) corr_e(2) ! GVB-BCCCnb total energy
  close(fid)

 !case('nwchem')
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_mrcc_energy_from_output: invalid&
                  & mrcc_prog='//TRIM(mrcc_prog)
  stop
 end select
end subroutine read_mrcc_energy_from_output

! perform GVB-BCCC2b/3b calculations
subroutine do_gvb_bccc()
 use mol, only: mult, ndb, npair, nopen, mrcc_e
 use mr_keyword, only: mem, nproc, mrcc_type, mrcc_prog, datname
 implicit none
 integer :: i, j, npair1, system, RENAME
 real(kind=8) :: ref_e, corr_e(2)
 real(kind=8), parameter :: e_H = -0.454397401648d0 ! ROHF/STO-2G
 character(len=24) :: data_string
 character(len=240) :: fchname1, fchname2, datname2, pyname, outname
 character(len=240) :: fcidump, ampname1, ampname2, inpname

 write(6,'(//,A)') 'Enter subroutine do_gvb_bccc...'
 write(6,'(A)',advance='no') 'GVB-BCCC'

 select case(mrcc_type)
 case(6)
  write(6,'(A)',advance='no') '2'
  mrcc_prog = 'gvb_bccc2b'
 case(7)
  write(6,'(A)',advance='no') '3'
  mrcc_prog = 'gvb_bccc3b'
 case default
  write(6,'(A,I0)') 'ERROR in subroutine do_gvb_bccc: invalid mrcc_type=',&
                     mrcc_type
  stop
 end select

 write(6,'(A)') 'b based on GVB orbitals'
 write(6,'(A)') 'Frozen_core = T. Frozen_vir = T. Using program gvb_bccc'

 ! print warnings for non-singlet calculations
 if(mult > 1) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Remark: non-singlet dectected. The molecule will be augmented&
                 & by proper number'
  write(6,'(A)') 'of H atoms far away. Then it becomes a singlet complex, and s&
                 &inglet GVB-BCCC w-'
  write(6,'(A)') 'ill be invoked. ROHF energies of H atoms using STO-2G have be&
                 &en subtracted in'
  write(6,'(A)') 'the following energies.'
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning: this algorithm cannot deal with a molecule whose gro&
                 &und state is sing-'
  write(6,'(A)') 'let. For example, if you use this algorithm to calculate the &
                 &triplet of NH3 mo-'
  write(6,'(A)') 'lecule, the GVB-BCCC wave function will collapse to singlet.'
  write(6,'(A)') REPEAT('-',79)
 end if

 i = index(datname, '_s.dat', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine do_gvb_bccc: '_s.dat' suffix not found&
                   & in filname "//TRIM(datname)
  stop
 end if
 fchname1 = datname(1:i+1)//'.fch'

 if(mult == 1) then
  ! singlet: generate FCIDUMP and perform linearized BCCC2b
  fchname2 = datname(1:i-1)//'.fch'
  datname2 = datname(1:i-1)//'.dat'
  pyname = datname(1:i-1)//'.py'
  outname = datname(1:i-1)//'.out'
  fcidump = datname(1:i-1)//'.FCIDUMP'
  call copy_file(fchname1, fchname2, .false.)
  i = system('dat2fch '//TRIM(datname2)//' '//TRIM(fchname2))
  i = 2*npair + nopen
  call prt_py_script_to_gen_fcidump(fchname2, i, i, mem, nproc)
  call submit_pyscf_job(pyname)
  call delete_files(3, [pyname, outname, fchname2])
  npair1 = npair

 else ! non-singlet: localization, add H atoms, and generate FCIDUMP
  pyname = datname(1:i+1)//'.py'
  outname = datname(1:i+1)//'.out'
  datname2 = datname(1:i-1)//'aH.dat'
  fcidump = datname(1:i-1)//'aH.FCIDUMP'
  call prt_py_script_loc_add_gen(fchname1, ndb, npair, nopen)
  call submit_pyscf_job(pyname)
  call delete_file(outname)
  npair1 = npair + nopen
 end if

 i = index(datname, '_s.dat', back=.true.)
 inpname = datname(1:i-1)//'_bccc.input'
 outname = datname(1:i-1)//'_linbccc2b.out'
 ampname1 = datname(1:i-1)//'_linbccc2b.amp'
 ampname2 = 'aaa'   ! no initial guess, nonexistent filename
 call prt_bccc_inp(1, npair1, inpname, fcidump, datname2, ampname2)
 call submit_gvb_bcci_job(nproc, 2, inpname, outname)
 ampname2 = datname(1:i-1)//'_bccc.amp'
 j = RENAME(TRIM(ampname2), TRIM(ampname1))

 ! perform GVB-BCCC2b
 outname = datname(1:i-1)//'_bccc2b.out'
 call prt_bccc_inp(1, npair1, inpname, fcidump, datname2, ampname1)
 call submit_gvb_bccc_job(1, nproc, 2, inpname, outname)
 ampname1 = datname(1:i-1)//'_bccc2b.amp'
 j = RENAME(TRIM(ampname2), TRIM(ampname1))
 call read_mrcc_energy_from_output(mrcc_prog, 4, outname, ref_e, corr_e)
 mrcc_e = corr_e(2)
 if(mult /= 1) then
  ref_e = ref_e - DBLE(nopen)*e_H
  corr_e(2) = corr_e(2) - DBLE(nopen)*e_H
 end if
 write(6,'(/,A,F18.8,1X,A4)')'E(GVB)        = ', ref_e    , 'a.u.'
 write(6,'(A,F18.8,1X,A4)')  'E(E_2bcorr)   = ', corr_e(1), 'a.u.'
 write(6,'(A,F18.8,1X,A4)')  'E(GVB-BCCC2b) = ', corr_e(2), 'a.u.'

 if(mrcc_type > 6) then   ! perform GVB-BCCC3b
  outname = datname(1:i-1)//'_bccc3b.out'
  call prt_bccc_inp(1, npair1, inpname, fcidump, datname2, ampname1)
  call submit_gvb_bccc_job(1, nproc, 3, inpname, outname)
  ampname2 = datname(1:i-1)//'_bccc.amp'
  ampname1 = datname(1:i-1)//'_bccc3b.amp'
  j = RENAME(TRIM(ampname2), TRIM(ampname1))
  call read_mrcc_energy_from_output(mrcc_prog, 5, outname, ref_e, corr_e)
  mrcc_e = corr_e(2)
  if(mult /= 1) then
   ref_e = ref_e - DBLE(nopen)*e_H
   corr_e(2) = corr_e(2) - DBLE(nopen)*e_H
  end if
  write(6,'(A,F18.8,1X,A4)')  'E(E_3bcorr)   = ', corr_e(1), 'a.u.'
  write(6,'(A,F18.8,1X,A4)')  'E(GVB-BCCC3b) = ', corr_e(2), 'a.u.'
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_gvb_bccc at '//TRIM(data_string)
end subroutine do_gvb_bccc

! print Python script to generate FCIDUMP
subroutine prt_py_script_to_gen_fcidump(fchname, nacto, nacte, mem, nproc)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nacto, nacte, mem, nproc
! nacto: number of active orbitals
! nacte: number of active electrons
 character(len=240) :: pyname
 character(len=240), intent(in) :: fchname

 i = index(fchname, '.fch', back=.true.)
 pyname = fchname(1:i-1)//'.py'

 open(newunit=fid,file=TRIM(pyname),status='replace')
 write(fid,'(A)') 'from mokit.lib.gaussian import gen_fcidump'
 write(fid,'(4(A,I0),A)') "gen_fcidump('"//TRIM(fchname)//"',", nacto, ',', &
                           nacte, ',', mem*1000, ',', nproc, ')'
 close(fid)
end subroutine prt_py_script_to_gen_fcidump

! print Python script to perform orbital localization for singly occupied
! open-shell orbitals, add H atoms and generate FCIDUMP for non-singlet
subroutine prt_py_script_loc_add_gen(fchname, ndb, npair, nopen)
 use mr_keyword, only: mem, nproc
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: ndb, npair, nopen
 character(len=240) :: pyname, lmofch, addH_fch, inpname, addH_dat, datname,&
                       addH_fch2
 character(len=240), intent(in) :: fchname

 k = ndb + npair
 i = index(fchname, '_s.fch', back=.true.)
 if(i == 0) then
  write(6,'(/,A)') "ERROR in subroutine prt_py_script_loc_add_gen: no '_s.fch'&
                   & suffix found in"
  write(6,'(A)') 'filename '//TRIM(fchname)
  stop
 end if
 pyname   = fchname(1:i+1)//'.py'
 lmofch   = fchname(1:i+1)//'_LMO.fch'
 addH_fch = fchname(1:i+1)//'_LMOaH.fch'
 inpname  = fchname(1:i+1)//'_LMOaH.inp'
 addH_dat = fchname(1:i-1)//'aH.dat'
 addH_fch2= fchname(1:i-1)//'aH.fch'
 datname  = fchname(1:i+1)//'.dat'

 open(newunit=fid,file=TRIM(pyname),status='replace')
 write(fid,'(A)') 'from mokit.lib.gaussian import loc, gen_fcidump'
 write(fid,'(A)') 'from mokit.lib.rwgeom import copy_and_add_pair_coeff'
 write(fid,'(A)') 'from os import system, rename'
 write(fid,'(A)') 'from shutil import copyfile'
 write(fid,'(/,2(A,I0),A)') "loc(fchname='"//TRIM(fchname)//"',idx=range(",&
                             k, ',', k+nopen, '))'
 write(fid,'(A)') "system('addH2singlet "//TRIM(lmofch)//" -gvb')"
 write(fid,'(A,I0,A)') "system('fch2inp "//TRIM(addH_fch)//" -gvb ",npair+nopen,"')"
 write(fid,'(A)') "rename('"//TRIM(inpname)//"','"//TRIM(addH_dat)//"')"
 write(fid,'(A,I0,A)') "copy_and_add_pair_coeff('"//TRIM(addH_dat)//"','"//&
                        TRIM(datname)//"',",nopen,')'
 write(fid,'(A)') "copyfile('"//TRIM(addH_fch)//"','"//TRIM(addH_fch2)//"')"
 write(fid,'(A)') "system('dat2fch "//TRIM(addH_dat)//' '//TRIM(addH_fch2)//"')"
 i = 2*(npair + nopen)
 write(fid,'(4(A,I0),A)') "gen_fcidump('"//TRIM(addH_fch2)//"',",i,',',i, ',',&
                           mem*1000, ',', nproc, ')'
 close(fid)
end subroutine prt_py_script_loc_add_gen

! print GVB-BCCC input file
subroutine prt_bccc_inp(mult, npair, inpname, fcidump, datname, ampname)
 implicit none
 integer :: fid
 integer, intent(in) :: mult, npair
 character(len=240), intent(in) :: inpname, fcidump, datname, ampname

 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A)') 'fcidump = '//TRIM(fcidump)
 write(fid,'(A)') 'datname = '//TRIM(datname)
 write(fid,'(A)') 'ampname = '//TRIM(ampname)

 select case(mult)
 case(1)
  write(fid,'(A,I0)') 'npair = ', npair
 case(3)
  write(fid,'(A,I0)') 'npair = ', npair+1
  write(fid,'(A,I0)') 'nif = ', 2*(npair+1)
 case default
  write(fid,'(A,I0)') 'ERROR in subroutine prt_bccc_inp: invalid mult=',mult
  stop
 end select

 close(fid)
end subroutine prt_bccc_inp

