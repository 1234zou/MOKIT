! written by jxzou at 20210619: move MRCI related subroutines here from automr.f90

! do uncontracted/ic-/FIC- MRCISD(+Q) for npair<=7, or <=CAS(14,14)
subroutine do_mrcisd()
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, casci, casscf, CIonly, ist, hf_fch, mrcisd,&
  mrcisd_prog, CtrType, casnofch, molcas_path, orca_path, gau_path, molpro_path,&
  psi4_path, bgchg, casci_prog, casscf_prog, chgname, F12, RI
 use mol, only: nbf, nif, npair, nopen, npair0, ndb, casci_e, casscf_e, davidson_e,&
  mrcisd_e, ptchg_e, nuc_pt_e
 use util_wrapper, only: unfchk, mkl2gbw
 integer :: i, system
 real(kind=8) :: e
 character(len=24) :: data_string
 character(len=240) :: string, chkname, inpname, outname, mklname
 character(len=47) :: error_warn='ERROR in subroutine do_mrcisd: invalid CtrType='

 if(.not. mrcisd) return
 write(iout,'(//,A)') 'Enter subroutine do_mrcisd...'

 if(npair0 > 7) then
  write(iout,'(A)') 'ERROR in subroutine do_mrcisd: reference wfn larger than (14,14).'
  write(iout,'(A)') 'DMRG-MRCISD is not supported currently.'
  stop
 end if

 if(.not. CIonly) then
  if(TRIM(casscf_prog) == 'orca') then
   write(iout,'(A)') 'Warning: ORCA is used as the CASSCF solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(iout,'(A)') 'This will affect the CI energy up to 10^-5 a.u.&
                    & Such small error is usually not important.'
   write(iout,'(A)') 'If you care about the accuracy, please use another CASSCF solver.'
  end if

  string = 'MRCISD based on optimized CASSCF orbitals.'
 else ! CIonly = .True.

  if(TRIM(casci_prog) == 'orca') then
   write(iout,'(A)') 'Warning: ORCA is used as the CASCI solver,&
                    & the NO coefficients in .mkl file are only 7-digits.'
   write(iout,'(A)') 'This will affect the CI energy up to 10^-5 a.u.&
                    & Such small error is usually not important.'
   write(iout,'(A)') 'If you care about the accuracy, please use another CASCI solver.'
  end if

  string = 'MRCISD based on CASCI orbitals.'
  write(iout,'(A)') 'Warning: the CASSCF orbital optimization is strongly recommended'
  write(iout,'(A)') 'to be performed before MRCI, unless it is too time-consuming.'
 end if
 write(iout,'(A)') TRIM(string)

 write(iout,'(A)') 'Frozen_Core = F, MRCISD computation using program '//TRIM(mrcisd_prog)

 select case(TRIM(mrcisd_prog))
 case('openmolcas')
  call check_exe_exist(molcas_path)
  i = system('fch2inporb '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  chkname = casnofch(1:i-1)//'.INPORB'
  string  = casnofch(1:i-1)//'.input'
  i = index(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'MRCISD.INPORB'
  inpname = casnofch(1:i)//'MRCISD.input'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(chkname), TRIM(mklname))
  i = RENAME(TRIM(string), TRIM(inpname))
  chkname = ' '
  call prt_mrcisd_molcas_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  i = system(TRIM(molcas_path)//' '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1")

 case('orca')
  call check_exe_exist(orca_path)
  i = system('fch2mkl '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  chkname = casnofch(1:i-1)//'_o.mkl'
  string  = casnofch(1:i-1)//'_o.inp'
  i = index(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'MRCISD.mkl'
  inpname = casnofch(1:i)//'MRCISD.inp'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(chkname), TRIM(mklname))
  i = RENAME(TRIM(string), TRIM(inpname))
  chkname = ' '
  call prt_mrcisd_orca_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(mklname)
  i = system(TRIM(orca_path)//' '//TRIM(inpname)//' >'//TRIM(outname)//" 2>&1")

 case('gaussian')
  call check_exe_exist(gau_path)
  i = index(casnofch, '_NO', back=.true.)
  chkname = casnofch(1:i)//'MRCISD.chk'
  inpname = casnofch(1:i)//'MRCISD.gjf'
  outname = casnofch(1:i)//'MRCISD.log'
  call unfchk(casnofch, chkname)
  call prt_mrcisd_gau_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  i = system(TRIM(gau_path)//' '//TRIM(inpname))

 case('molpro')
  call check_exe_exist(molpro_path)

  i = system('fch2com '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'.com'
  chkname = casnofch
  call convert2molpro_fname(chkname, '.a')

  i = index(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.com'
  outname = casnofch(1:i)//'MRCISD.out'
  mklname = inpname
  call convert2molpro_fname(mklname, '.a')
  i = RENAME(TRIM(string), TRIM(inpname))
  i = RENAME(TRIM(chkname), TRIM(mklname))

  call prt_mrcisd_molpro_inp(inpname)
  i = CEILING(DBLE(mem*125)/DBLE(nproc))
  write(string,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -m ', i,'m '//TRIM(inpname)
  write(iout,'(A)') '$'//TRIM(string)
  i = system(TRIM(string))

 case('psi4')
  call check_exe_exist(psi4_path)

  i = system('fch2psi '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'_psi.inp'

  i = index(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.inp'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(string), TRIM(inpname))

  call prt_mrcisd_psi4_inp(inpname)
  write(string,'(A,I0)') 'psi4 '//TRIM(inpname)//' '//TRIM(outname)&
                         //' -n ',nproc
  write(iout,'(A)') '$'//TRIM(string)
  i = system(TRIM(string))

 case('dalton')
  i = system('fch2dal '//TRIM(casnofch))
  i = index(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'.dal'
  chkname = casnofch(1:i-1)//'.mol'
  i = index(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.dal'
  mklname = casnofch(1:i)//'MRCISD.mol'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(string), TRIM(inpname))
  i = RENAME(TRIM(chkname), TRIM(mklname))
  call prt_mrcisd_dalton_inp(inpname)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  i = index(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'MRCISD'
  chkname = casnofch(1:i)//'MRCISD.sout'
  write(string,'(2(A,I0),A)') 'dalton -gb ',mem,' -omp ',nproc,' -ow '//TRIM(mklname)&
                              //' >'//TRIM(chkname)//" 2>&1"
  write(iout,'(A)') '$'//TRIM(string)
  i = system(TRIM(string))

 case default
  write(iout,'(A)') 'ERROR in subroutine do_mrcisd: invalid program='//TRIM(mrcisd_prog)
  stop
 end select

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine do_mrcisd: error termination.'
  write(iout,'(A)') 'Filename='//TRIM(outname)
  stop
 end if

 ! read Davidson correction and MRCISD energy from OpenMolcas/ORCA/Gaussian output file
 call read_mrcisd_energy_from_output(CtrType, mrcisd_prog, outname, ptchg_e,&
      nuc_pt_e, davidson_e, e)
 mrcisd_e = e + davidson_e ! E(MRCISD+Q)
 if(CIonly) then   ! E(MRCISD) - (E(CASCI) or E(CASSCF))
  e = e - casci_e
 else
  e = e - casscf_e
 end if

 select case(TRIM(mrcisd_prog))
 case('openmolcas')
  select case(CtrType)
  case(1) ! uncontracted MRCISD
   write(iout,'(/,A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E(MRCISD)      =', mrcisd_e, 'a.u.'
  case(2) ! ic-MRCISD
   write(iout,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E_corr(icMRCISD) =', e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E(icMRCISD+Q)    =', mrcisd_e, 'a.u.'
  case default
   write(iout,'(A,I0)') error_warn, CtrType
   stop
  end select
 case('orca')
  write(iout,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
  select case(CtrType)
  case(1) ! uncontracted MRCISD
   write(iout,'(A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E(MRCISD+Q)    =', mrcisd_e, 'a.u.'
  case(3) ! FIC-MRCISD
   write(iout,'(A,F18.8,1X,A4)') 'E_corr(FIC-MRCISD) =', e, 'a.u.'
   write(iout,'(A,F18.8,1X,A4)') 'E(FIC-MRCISD+Q)    =', mrcisd_e, 'a.u.'
  case default
   write(iout,'(A,I0)') error_warn, CtrType
   stop
  end select
 case('gaussian','psi4','dalton') ! only uncontracted MRCISD
  if(CtrType /= 1) then
   write(iout,'(A,I0)') error_warn, CtrType
   stop
  end if
  write(iout,'(/,A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(MRCISD)      =', mrcisd_e, 'a.u.'
 case('molpro')
  write(iout,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
  write(iout,'(A,F18.8,1X,A4)') 'E(MRCISD)      =', mrcisd_e, 'a.u.'
 end select

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_mrcisd at '//TRIM(data_string)
 return
end subroutine do_mrcisd

! print MRCISD keywords into OpenMolcas .input file
subroutine prt_mrcisd_molcas_inp(inpname)
 use print_id, only: iout
 use mol, only: nif, ndb, nopen, nacta, nactb, npair, npair0, charge, mult
 use mr_keyword, only: CtrType
 implicit none
 integer :: i, idx, nvir, ne, fid
 character(len=240), intent(in) :: inpname
 character(len=240) :: buf, inporb

 i = index(inpname, '.input', back=.true.)
 inporb = inpname(1:i-1)//'.INPORB'

 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(1:4)=="&SCF" .or. buf(1:7)=="&RASSCF") exit
 end do ! for while

 BACKSPACE(fid)
 write(fid,'(/,A)') "&RASSCF"
 write(fid,'(A,I0)') 'Charge = ', charge
 write(fid,'(A,I0)') 'Spin = ', mult
 write(fid,'(A)') 'FILEORB = '//TRIM(inporb)

 idx = ndb + npair - npair0 ! doubly occupied orbitals in CAS
 nvir = nif - idx - 2*npair0 - nopen ! virtual orbitals in CAS
 ne = 2*idx + nacta + nactb ! all electrons

 if(CtrType == 1) then ! uncontracted MRCISD
  write(fid,'(A,I0,A)') 'nActEl = ', ne, ' 2 2'
  write(fid,'(A,I0)') 'RAS1 = ', idx
  write(fid,'(A,I0)') 'RAS2 = ', 2*npair0+nopen
  write(fid,'(A,I0)') 'RAS3 = ', nvir
  write(fid,'(A,/)') 'CIonly'
 else if(CtrType == 2) then ! ic-MRCISD
  write(fid,'(A,I0,A)') 'nActEl = ', nacta+nactb, ' 0 0'
  write(fid,'(A,I0)') 'RAS2 = ', 2*npair0+nopen
  write(fid,'(A)') 'CIonly'

  write(fid,'(/,A)') "&MOTRA"
  write(fid,'(A)') 'LUMORB'
  write(fid,'(A)') 'Frozen = 0'

  write(fid,'(/,A)') "&GUGA"
  write(fid,'(A,I0)') 'Spin = ', mult
  write(fid,'(A,I0)') 'Inactive = ', idx
  write(fid,'(A,I0)') 'Active = ', 2*npair0+nopen
  write(fid,'(A,I0)') 'NACTel = ', nacta+nactb
  write(fid,'(A,I0)') 'CIAll = 1'

  write(fid,'(/,A)') "&MRCI"
  write(fid,'(A,/)') 'MAXIterations = 49'
 end if

 close(fid)
 return
end subroutine prt_mrcisd_molcas_inp

! print MRCISD keywords into ORCA .inp file
subroutine prt_mrcisd_orca_inp(inpname1)
 use print_id, only: iout
 use mol, only: nopen, nacta, nactb, npair0, mult
 use mr_keyword, only: mem, nproc, CtrType, DKH2
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

 if(CtrType == 1) then ! uncontracted MRCISD
  write(fid2,'(A)') '%mrci'
  write(fid2,'(2(A,I0),A)') ' NewBlock 1 * nroots 1 refs cas(',nacta+nactb,',',&
                             2*npair0+nopen,') end end'
  write(fid2,'(A)') ' tsel 0.0'
  write(fid2,'(A)') ' tpre 0.0'
  write(fid2,'(A)') ' Etol 1e-7'
  write(fid2,'(A)') ' Rtol 1e-7'
 else if(CtrType == 3) then ! FIC-MRCISD
  write(fid2,'(A)') '%autoci'
  write(fid2,'(A)') ' CItype FICMRCI'
  write(fid2,'(A,I0)') ' nel ',  nacta+nactb
  write(fid2,'(A,I0)') ' norb ', 2*npair0+nopen
  write(fid2,'(A,I0)') ' mult ', mult
  write(fid2,'(A)') ' nroots 1'
  write(fid2,'(A)') ' DavidsonOpt 1'
 end if

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
end subroutine prt_mrcisd_orca_inp

! print MRCISD keywords into Gaussian .gjf file
subroutine prt_mrcisd_gau_inp(gjfname)
 use mol, only: nif, ndb, nopen, nacta, nactb, npair, npair0
 use mr_keyword, only: mem, nproc, DKH2
 implicit none
 integer :: i, ne, nvir, idx, fid
 character(len=240), intent(in) :: gjfname

 idx = ndb + npair - npair0 ! doubly occupied orbitals in CAS
 nvir = nif - idx - 2*npair0 - nopen ! virtual orbitals in CAS
 ne = 2*idx + nacta + nactb ! all electrons

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 i = index(gjfname, '.gjf', back=.true.)
 write(fid,'(A,I0)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A5,I0,A2)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(4(A,I0),A)',advance='no') '#p CAS(',ne,',',nif,',ras(2,',idx,',2,',&
  nvir,'))/chkbasis nosymm guess=read geom=allcheck scf(maxcycle=-1)'

 if(DKH2) then
  write(fid,'(A,/)') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A,/)') ' int=nobasistransform'
 end if

 close(fid)
 return
end subroutine prt_mrcisd_gau_inp

! print MRCISD keywords into Molpro input file
subroutine prt_mrcisd_molpro_inp(inpname)
 implicit none
 integer :: i, fid
 character(len=240), intent(in) :: inpname

 call prt_cas_molpro_inp(inpname, .false., .false.)
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 BACKSPACE(fid)
 write(fid,'(A)') '{MRCIC;CORE}'
 close(fid)

 return
end subroutine prt_mrcisd_molpro_inp

! print MRCISD keywords into a PSI4 input file
subroutine prt_mrcisd_psi4_inp(inpname)
 use mol, only: nif, ndb, nopen, nacto, npair, npair0
 use mr_keyword, only: iout, mem
 implicit none
 integer :: i, fid, idx, nvir
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 idx = ndb + npair - npair0 ! doubly occupied orbitals in CAS
 nvir = nif - idx - 2*npair0 - nopen ! virtual orbitals in CAS

 call modify_memory_in_psi4_inp(inpname, mem)
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(index(buf,'d_convergence') > 0) exit

  if(buf(1:6) == 'memory') then
   write(iout,'(A)') 'ERROR in subroutine prt_mrcisd_psi4_inp: incomplete&
                    & file '//TRIM(inpname)
   close(fid)
   stop
  end if
 end do ! for while

 write(fid,'(A)') ' ci_maxiter 100'
 write(fid,'(A,I0,A)') ' ras1 [',idx  ,']'
 write(fid,'(A,I0,A)') ' ras2 [',nacto,']'
 write(fid,'(A,I0,A)') ' ras3 [',nvir ,']'
 write(fid,'(A)') '}'

 write(fid,'(/,A)') "mrcisd_energy = energy('detci',ref_wfn=scf_wfn)"
 close(fid)
 return
end subroutine prt_mrcisd_psi4_inp

! print MRCISD keywords into a Dalton input file
subroutine prt_mrcisd_dalton_inp(inpname)
 use mol, only: mult, nif, ndb, nopen, nacto, nacte, npair, npair0
 use mr_keyword, only: iout
 implicit none
 integer :: i, fid, fid1, idx, nvir, RENAME
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 idx = ndb + npair - npair0 ! doubly occupied orbitals in CAS
 nvir = nif - idx - 2*npair0 - nopen ! virtual orbitals in CAS

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname1),status='replace')
 write(fid1,'(A,/,A)') '**DALTON INPUT','.RUN WAVE FUNCTIONS'
 write(fid1,'(A)') '**WAVE FUNCTIONS'
 write(fid1,'(A)') '.CI','*CONFIGURATION INPUT'
 write(fid1,'(A,/,I0)') '.SYMMETRY', 1
 write(fid1,'(A,/,I0)') '.SPIN MULTIPLICITY',mult
 write(fid1,'(A,/,I0)') '.ELECTRONS', nacte+2*idx
 write(fid1,'(A,/,I0)') '.INACTIVE', 0
 write(fid1,'(A,/,I0)') '.RAS1 SPACE', idx
 write(fid1,'(A,/,A)') '.RAS1 HOLE2', '0 2'
 write(fid1,'(A,/,I0)') '.RAS2 SPACE', nacto
 write(fid1,'(A,/,I0)') '.RAS3 SPACE', nvir
 write(fid1,'(A,/,A)') '.RAS3 ELECTRONS', '0 2'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12) == '*ORBITAL INP') exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(A)') "ERROR in subroutine prt_mrcisd_dalton_inp: no '*ORBITAL&
                   & INP' found in file "//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(A)') TRIM(buf)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 return
end subroutine prt_mrcisd_dalton_inp

