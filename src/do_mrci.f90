! written by jxzou at 20210619: move MRCI related subroutines here from automr.f90

! do uncontracted/ic-/FIC- MRCISD(+Q) for npair<=7, or <=CAS(14,14)
subroutine do_mrcisd()
 use mr_keyword, only: mem, nproc, dmrgci, dmrgscf, dmrg_no, CIonly, eist, mrcisd,&
  mrcisd_prog, CtrType, casnofch, molcas_omp, molcas_path, orca_path, gau_path,&
  gms_path, gms_scr_path, molpro_path, psi4_path, dalton_mpi, bgchg, casci_prog,&
  casscf_prog, chgname
 use mol, only: casci_e, casscf_e, davidson_e, mrcisd_e, ptchg_e, nuc_pt_e
 use util_wrapper, only: bas_fch2py_wrap, unfchk, mkl2gbw, fch2inporb_wrap
 implicit none
 integer :: i, system, RENAME
 real(kind=8) :: e
 character(len=24) :: data_string
 character(len=240) :: string, cmofch, pyname, chkname, inpname, outname, mklname
 character(len=480) :: datpath
 character(len=47) :: error_warn='ERROR in subroutine do_mrcisd: invalid CtrType='

 if(eist == 1) return ! excited state calculation
 if(.not. mrcisd) return
 write(6,'(//,A)') 'Enter subroutine do_mrcisd...'

 call prt_mrci_orb_type(2, CIonly)

 if(.not. CIonly) then
  if(TRIM(casscf_prog) == 'orca') then
   write(6,'(A)') 'Warning: ORCA is used as the CASSCF solver, the NO coefficie&
                  &nts in .mkl file are only 7-digits.'
   write(6,'(A)') 'This will affect the CI energy up to 10^-5 a.u. Such small e&
                  &rror is usually not important.'
   write(6,'(A)') 'If you care about the accuracy, please use another CASSCF so&
                  &lver.'
  end if
 else ! CIonly = .True.
  if(TRIM(casci_prog) == 'orca') then
   write(6,'(A)') 'Warning: ORCA is used as the CASCI solver, the NO coefficien&
                  &ts in .mkl file are only 7-digits.'
   write(6,'(A)') 'This will affect the CI energy up to 10^-5 a.u. Such small e&
                  &rror is usually not important.'
   write(6,'(A)') 'If you care about the accuracy, please use another CASCI sol&
                  &ver.'
  end if
 end if

 write(6,'(A)') 'Frozen_Core = F, MRCISD computation using program '//&
                 TRIM(mrcisd_prog)

 select case(TRIM(mrcisd_prog))
 case('pyscf')
  if(.not. (dmrgci .or. dmrgscf)) then
   write(6,'(/,A)') 'ERROR in subroutine do_mrcisd: for MRCISD_prog=PySCF, only&
                    & DMRG-FIC-MRCISD is'
   write(6,'(A)') 'supported.'
   stop
  end if
  ! For DMRG-NEVPT2, use CMOs rather than NOs
  if(dmrgci .or. dmrgscf) then ! DMRG-CASCI/CASSCF
   if(dmrg_no) then ! DMRG NOs calculated previously
    i = INDEX(casnofch, '_NO', back=.true.)
    cmofch = casnofch(1:i)//'CMO.fch'
    casnofch = cmofch
   end if
   ! if DMRG NOs not calculated previously, casnofch is set to 'xxx_CMO.fch'
   ! in subroutine do_cas
   i = INDEX(casnofch, '_CMO', back=.true.)
  else                         ! CASCI/CASSCF
   i = INDEX(casnofch, '_NO', back=.true.)
  end if
  pyname = casnofch(1:i)//'MRCISD.py'
  outname = casnofch(1:i)//'MRCISD.out'
  call bas_fch2py_wrap(casnofch, .false., pyname)
  call prt_mrci_script_into_py(pyname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(pyname))
  call submit_pyscf_job(pyname, .true.)

 case('openmolcas')
  call check_exe_exist(molcas_path)
  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.input'
  outname = casnofch(1:i)//'MRCISD.out'
  call fch2inporb_wrap(casnofch, .false., inpname)
  call prt_mrci_molcas_inp(2, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)

 case('orca')
  call check_exe_exist(orca_path)
  i = SYSTEM('fch2mkl '//TRIM(casnofch))
  i = INDEX(casnofch, '.fch', back=.true.)
  chkname = casnofch(1:i-1)//'_o.mkl'
  string  = casnofch(1:i-1)//'_o.inp'
  i = INDEX(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'MRCISD.mkl'
  inpname = casnofch(1:i)//'MRCISD.inp'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(chkname), TRIM(mklname))
  i = RENAME(TRIM(string), TRIM(inpname))
  chkname = ' '
  call prt_mrcisd_orca_inp(inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(mklname)
  call submit_orca_job(orca_path, inpname, .true., .false., .false.)

 case('gaussian')
  call check_exe_exist(gau_path)
  i = INDEX(casnofch, '_NO', back=.true.)
  chkname = casnofch(1:i)//'MRCISD.chk'
  inpname = casnofch(1:i)//'MRCISD.gjf'
  outname = casnofch(1:i)//'MRCISD.log'
  call unfchk(casnofch, chkname)
  call prt_mrcisd_gau_inp(2, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_gau_job(gau_path, inpname, .true.)

 case('molpro')
  call check_exe_exist(molpro_path)

  i = SYSTEM('fch2com '//TRIM(casnofch))
  i = INDEX(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'.com'
  chkname = casnofch
  call convert2molpro_fname(chkname, '.a')

  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.com'
  outname = casnofch(1:i)//'MRCISD.out'
  mklname = inpname
  call convert2molpro_fname(mklname, '.a')
  i = RENAME(TRIM(string), TRIM(inpname))
  i = RENAME(TRIM(chkname), TRIM(mklname))

  call prt_mrcisd_molpro_inp(inpname)
  i = CEILING(DBLE(mem*125)/DBLE(nproc))
  write(datpath,'(2(A,I0),A)') TRIM(molpro_path)//' -n ',nproc,' -t 1 -m ',i,&
                               'm '//TRIM(inpname)

 case('psi4')
  call check_exe_exist(psi4_path)

  i = SYSTEM('fch2psi '//TRIM(casnofch))
  i = INDEX(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'_psi.inp'

  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.inp'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(string), TRIM(inpname))

  call prt_mrci_psi4_inp(2, inpname)
  call submit_psi4_job(psi4_path, inpname, nproc)

 case('dalton')
  i = SYSTEM('fch2dal '//TRIM(casnofch))
  i = INDEX(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'.dal'
  chkname = casnofch(1:i-1)//'.mol'
  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.dal'
  mklname = casnofch(1:i)//'MRCISD.mol'
  outname = casnofch(1:i)//'MRCISD.out'
  i = RENAME(TRIM(string), TRIM(inpname))
  i = RENAME(TRIM(chkname), TRIM(mklname))
  call prt_mrci_dalton_inp(2, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  i = INDEX(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'MRCISD'
  chkname = casnofch(1:i)//'MRCISD.sout'
  call submit_dalton_job(mklname,mem,nproc,dalton_mpi,.false.,.false.,.false.)

 case('gamess')
  i = SYSTEM('fch2inp '//TRIM(casnofch))
  i = INDEX(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'.inp'
  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISD.inp'
  mklname = casnofch(1:i)//'MRCISD.dat'
  outname = casnofch(1:i)//'MRCISD.gms'
  i = RENAME(TRIM(string), TRIM(inpname))
  call prt_mrci_gms_inp(2, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_gms_job(gms_path, gms_scr_path, inpname, nproc)

 case default
  write(6,'(/,A)') 'ERROR in subroutine do_mrcisd: invalid program='//&
                   TRIM(mrcisd_prog)
  stop
 end select

 if(TRIM(mrcisd_prog) == 'molpro') then
  write(6,'(A)') '$'//TRIM(datpath)
  i = SYSTEM(TRIM(datpath))
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine do_mrcisd: Molpro MRCISD job failed.'
   write(6,'(A)') 'Please open file '//TRIM(outname)//' and check why.'
   stop
  end if
 end if

 ! read Davidson correction and MRCISD energy from OpenMolcas/ORCA/Gaussian/
 ! Dalton/GAMESS output file
 call read_mrci_energy_from_output(CtrType, mrcisd_prog, outname, ptchg_e,&
      nuc_pt_e, davidson_e, e)

 mrcisd_e = e + davidson_e ! E(MRCISD+Q)
 if(CIonly) then   ! E(MRCISD) - (E(CASCI) or E(CASSCF))
  e = e - casci_e
 else
  e = e - casscf_e
 end if

 select case(TRIM(mrcisd_prog))
 case('openmolcas', 'psi4')
  if(CtrType == 1) then
   call calc_davidson_corr_from_out(mrcisd_prog, outname, e, davidson_e)
   mrcisd_e = mrcisd_e + davidson_e
  end if
 end select

 select case(TRIM(mrcisd_prog))
 case('pyscf')
  write(6,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E_corr(FIC-MRCISD) =', e, 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(FIC-MRCISD+Q)    =', mrcisd_e, 'a.u.'
 case('openmolcas')
  write(6,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
  select case(CtrType)
  case(1) ! uncontracted MRCISD
   write(6,'(A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
   write(6,'(A,F18.8,1X,A4)') 'E(MRCISD+Q)    =', mrcisd_e, 'a.u.'
  case(2) ! ic-MRCISD
   write(6,'(A,F18.8,1X,A4)') 'E_corr(icMRCISD) =', e, 'a.u.'
   write(6,'(A,F18.8,1X,A4)') 'E(icMRCISD+Q)    =', mrcisd_e, 'a.u.'
  case default
   write(6,'(A,I0)') error_warn, CtrType
   stop
  end select
 case('orca')
  write(6,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
  select case(CtrType)
  case(1) ! uncontracted MRCISD
   write(6,'(A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
   write(6,'(A,F18.8,1X,A4)') 'E(MRCISD+Q)    =', mrcisd_e, 'a.u.'
  case(3) ! FIC-MRCISD
   write(6,'(A,F18.8,1X,A4)') 'E_corr(FIC-MRCISD) =', e, 'a.u.'
   write(6,'(A,F18.8,1X,A4)') 'E(FIC-MRCISD+Q)    =', mrcisd_e, 'a.u.'
  case default
   write(6,'(A,I0)') error_warn, CtrType
   stop
  end select
 case('gaussian','dalton') ! only uncontracted MRCISD
  if(CtrType /= 1) then
   write(6,'(A,I0)') error_warn, CtrType
   stop
  end if
  write(6,'(/,A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(MRCISD)      =', mrcisd_e, 'a.u.'
 case('molpro','gamess','psi4')
  write(6,'(/,A,F18.8,1X,A4)') 'Davidson correction=', davidson_e, 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E_corr(MRCISD) =', e, 'a.u.'
  write(6,'(A,F18.8,1X,A4)') 'E(MRCISD+Q)    =', mrcisd_e, 'a.u.'
 end select

 if(TRIM(mrcisd_prog) == 'gamess') then
  write(6,'(/,A)') 'You may notice that the E(MRCISD+Q) above is slightly&
                  & different with that in'
  write(6,'(A)') '.gms file. This is because GAMESS uses renormalized Davidson&
                 & size extensivity'
  write(6,'(A)') 'correction. While MOKIT adopts the simple Davidson size exte&
                 &nsivity correction'
  write(6,'(A)') 'E_corr*(1-c^2).'
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mrcisd at '//TRIM(data_string)
end subroutine do_mrcisd

! print DMRG-FIC-MRCISD script into a given .py file
subroutine prt_mrci_script_into_py(pyname)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, maxM, X2C, RI, RIJK_bas, iroot, target_root,&
  block_mpi
 implicit none
 integer :: i, nroots, fid1, fid2, RENAME
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, pyname1
 character(len=240), intent(in) :: pyname

 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
 pyname1 = TRIM(pyname)//'.t'

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:17) == 'from pyscf import') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 buf = TRIM(buf)//', mcscf, dmrgscf, lib'
 write(fid2,'(A)') TRIM(buf)
 write(fid2,'(A)') 'from pyblock2.icmr.icmrcisd_full import WickICMRCISD'

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 write(fid2,'(/,A)',advance='no') "dmrgscf.settings.MPIPREFIX = '"
 if(block_mpi) write(fid2,'(A,I0)',advance='no') 'mpirun -n ', nproc
 write(fid2,'(A)') "'"
 write(fid2,'(A,I0)') 'nproc = ', nproc
 write(fid2,'(A)') 'lib.num_threads(nproc)'

 do while(.true.)
  read(fid1,'(A)') buf
  write(fid2,'(A)') TRIM(buf)
  if(buf(1:12) == 'mf.max_cycle') exit
 end do ! for while
 write(fid2,'(A,I0,A)') 'mf.max_memory = ', mem*1000, ' # MB'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')

 write(fid2,'(A)') '# generate CASCI wfn'
 if(X2C) then
  write(fid2,'(A)',advance='no') 'mc = mcscf.CASCI(mf.x2c1e(),'
 else
  write(fid2,'(A)',advance='no') 'mc = mcscf.CASCI(mf,'
 end if
 write(fid2,'(3(I0,A))',advance='no') nacto,',(',nacta,',',nactb,')'
 if(RI) then
  write(fid2,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
 else
  write(fid2,'(A)') ')'
 end if

 write(fid2,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=', maxM, ')'
 if(block_mpi) then
  i = CEILING(0.5*REAL(mem)/REAL(nproc))
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', (mem-nproc*i)*1000, ' # MB'
 else
  i = CEILING(0.4*REAL(mem))
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', (mem-i)*1000, ' # MB'
  write(fid2,'(A)') 'mc.fcisolver.threads = nproc'
 end if
 write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ', i, ' # GB'

 if(iroot > 0) then ! FIC-MRCISD based on SS-CASSCF
  call read_ss_root_from_txt(nroots, target_root)
  write(fid2,'(A,I0)') 'mc.fcisolver.nroots = ', nroots
 end if

 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.kernel()'

 write(fid2,'(/,A)') 'mol.verbose = 5'
 write(fid2,'(A)') 'mm = WickICMRCISD(mc)'
 write(fid2,'(A)') 'mm.kernel()'
 close(fid2)
 i = RENAME(pyname1, pyname)
end subroutine prt_mrci_script_into_py

! print MRCISD/MRCISDT keywords into OpenMolcas .input file
subroutine prt_mrci_molcas_inp(order, inpname)
 use mol, only: nif, ndb, nopen, nacta, nactb, npair, npair0, charge, mult
 use mr_keyword, only: CtrType
 implicit none
 integer :: i, idx, nvir, ne, fid
 integer, intent(in) :: order ! 2/3 for MRCISD/MRCISDT
 character(len=240), intent(in) :: inpname
 character(len=240) :: buf, inporb

 i = INDEX(inpname, '.input', back=.true.)
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
  write(fid,'(A,I3,2I2)') 'nActEl = ', ne, order, order
  write(fid,'(A,I0)') 'RAS1 = ', idx
  write(fid,'(A,I0)') 'RAS2 = ', 2*npair0+nopen
  write(fid,'(A,I0)') 'RAS3 = ', nvir
  if(order == 2) write(fid,'(A)') 'PRWF = 1d-5'
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
end subroutine prt_mrci_molcas_inp

! print MRCISD keywords into ORCA .inp file
subroutine prt_mrcisd_orca_inp(inpname1)
 use mol, only: nopen, nacta, nactb, npair0, mult
 use mr_keyword, only: mem, nproc, CtrType, DKH2, hardwfn, crazywfn
 implicit none
 integer :: i, fid1, fid2
 integer :: RENAME
 character(len=240), intent(in) :: inpname1
 character(len=240) :: buf, inpname2

 inpname2 = TRIM(inpname1)//'.t'
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
  write(fid2,'(A)') ' CItype MRCI'
  write(fid2,'(2(A,I0),A)') ' NewBlock 1 * nroots 1 refs cas(',nacta+nactb,',',&
                             2*npair0+nopen,') end end'
  write(fid2,'(A)') ' tsel 0.0'
  write(fid2,'(A)') ' tpre 0.0'
  write(fid2,'(A)') ' Etol 1e-7'
  write(fid2,'(A)') ' Rtol 1e-7'
  if(hardwfn) then
   write(fid2,'(A)') ' NGuessMat 1700'
   write(fid2,'(A)') ' NGuessMatRefCI 1700'
  else if(crazywfn) then
   write(fid2,'(A)') ' NGuessMat 2500'
   write(fid2,'(A)') ' NGuessMatRefCI 2500'
  else
   write(fid2,'(A)') ' NGuessMat 1000'
   write(fid2,'(A)') ' NGuessMatRefCI 1000'
  end if
 else if(CtrType == 3) then ! FIC-MRCISD
  write(fid2,'(A)') '%autoci'
  write(fid2,'(A)') ' CItype FICMRCI'
  write(fid2,'(A,I0)') ' nel ',  nacta+nactb
  write(fid2,'(A,I0)') ' norb ', 2*npair0+nopen
  write(fid2,'(A,I0)') ' mult ', mult
  write(fid2,'(A)') ' nroots 1'
  write(fid2,'(A)') ' DavidsonOpt 1'
  ! NGuessMat can only be used in ORCA 6
 end if

 if(hardwfn) then
  write(fid2,'(A)') ' MaxIter 400'
 else if(crazywfn) then
  write(fid2,'(A)') ' MaxIter 600'
 else
  write(fid2,'(A)') ' MaxIter 200'
 end if
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
end subroutine prt_mrcisd_orca_inp

! print MRCISD/MRCISDT keywords into Gaussian .gjf file
subroutine prt_mrcisd_gau_inp(order, gjfname)
 use mol, only: nif, ndb, nopen, nacta, nactb, npair, npair0
 use mr_keyword, only: mem, nproc, DKH2
 implicit none
 integer :: i, ne, nvir, idx, fid
 integer, intent(in) :: order
 character(len=240), intent(in) :: gjfname

 idx = ndb + npair - npair0 ! doubly occupied orbitals in CAS
 nvir = nif - idx - 2*npair0 - nopen ! virtual orbitals in CAS
 ne = 2*idx + nacta + nactb ! all electrons

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 i = INDEX(gjfname, '.gjf', back=.true.)
 write(fid,'(A,I0)') '%chk='//gjfname(1:i-1)//'.chk'
 write(fid,'(A5,I0,A2)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(6(A,I0),A)',advance='no') '#p CASSCF(',ne,',',nif,',ras(',order,',',&
  idx,',',order,',',nvir,')) chkbasis nosymm guess=read geom=allcheck scf(maxcycle=-1)'

 if(DKH2) then
  write(fid,'(A,/)') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A,/)') ' int=nobasistransform'
 end if

 close(fid)
end subroutine prt_mrcisd_gau_inp

! print MRCISD keywords into Molpro input file
subroutine prt_mrcisd_molpro_inp(inpname)
 implicit none
 integer :: fid
 character(len=240), intent(in) :: inpname

 call prt_cas_molpro_inp(inpname, .false., .false.)
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 BACKSPACE(fid)
 write(fid,'(A)') '{MRCIC;CORE}'
 close(fid)
end subroutine prt_mrcisd_molpro_inp

! print MRCISD/MRCISDT keywords into a PSI4 input file
subroutine prt_mrci_psi4_inp(order, inpname)
 use mol, only: nif, ndb, nopen, nacto, npair, npair0
 use mr_keyword, only: mem
 implicit none
 integer :: fid, idx, nvir
 integer, intent(in) :: order
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
   write(6,'(A)') 'ERROR in subroutine prt_mrcisd_psi4_inp: incomplete&
                  & file '//TRIM(inpname)
   close(fid)
   stop
  end if
 end do ! for while

 write(fid,'(A)') ' ci_maxiter 100'
 write(fid,'(A,I0,A)') ' ras1 [',idx  ,']'
 write(fid,'(A,I0,A)') ' ras2 [',nacto,']'
 write(fid,'(A,I0,A)') ' ras3 [',nvir ,']'
 if(order == 2) write(fid,'(A)') ' num_dets_print 9999'
 write(fid,'(A,I0)') ' ex_level ', order
 write(fid,'(A)') '}'

 write(fid,'(/,A)') "mrci_energy = energy('detci',ref_wfn=scf_wfn)"
 close(fid)
end subroutine prt_mrci_psi4_inp

! print MRCISD/MRCISDT keywords into a Dalton input file
subroutine prt_mrci_dalton_inp(order, inpname)
 use mol, only: mult, nif, ndb, nopen, nacto, nacte, npair, npair0
 use mr_keyword, only: DKH2
 implicit none
 integer :: i, fid, fid1, idx, nvir, RENAME
 integer, intent(in) :: order
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 idx = ndb + npair - npair0 ! doubly occupied orbitals in CAS
 nvir = nif - idx - 2*npair0 - nopen ! virtual orbitals in CAS

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname1),status='replace')
 write(fid1,'(A)') '**DALTON INPUT'
 if(DKH2) write(fid1,'(A)') '.DOUGLAS-KROLL'
 write(fid1,'(A)') '.RUN WAVE FUNCTIONS'
 write(fid1,'(A)') '**WAVE FUNCTIONS'
 write(fid1,'(A)') '.CI','*CONFIGURATION INPUT'
 write(fid1,'(A,/,I0)') '.SYMMETRY', 1
 write(fid1,'(A,/,I0)') '.SPIN MULTIPLICITY',mult
 write(fid1,'(A,/,I0)') '.ELECTRONS', nacte+2*idx
 write(fid1,'(A,/,I0)') '.INACTIVE', 0
 write(fid1,'(A,/,I0)') '.RAS1 SPACE', idx
 write(fid1,'(A)') '.RAS1 HOLE2'
 write(fid1,'(A,I0)') '0 ', order
 write(fid1,'(A,/,I0)') '.RAS2 SPACE', nacto
 write(fid1,'(A,/,I0)') '.RAS3 SPACE', nvir
 write(fid1,'(A)') '.RAS3 ELECTRONS'
 write(fid1,'(A,I0)') '0 ', order

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12) == '*ORBITAL INP') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine prt_mrci_dalton_inp: no '*ORBITAL INP' &
                 &found in file "//TRIM(inpname)
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
end subroutine prt_mrci_dalton_inp

! add MRCISD+Q/MRCISDT keywords into a GAMESS input file
subroutine prt_mrci_gms_inp(order, inpname)
 use mr_keyword, only: mem, nproc, cart
 use mol, only: ndb, nacto, nacte, npair, npair0, nif, charge, mult
 implicit none
 integer :: i, ndb0, ne, fid, fid1, RENAME
 integer, intent(in) :: order
 character(len=240), intent(in) :: inpname
 character(len=240) :: buf, inpname1

 ndb0 = ndb + npair - npair0
 ne = 2*ndb0 + nacte
 inpname1 = TRIM(inpname)//'.t'

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')
 write(fid1,'(2(A,I0),A)') ' $CONTRL SCFTYP=NONE RUNTYP=ENERGY ICHARG=',&
                           charge, ' MULT=', mult, ' NOSYM=1 ICUT=11'
 if(cart) then
  write(fid1,'(A)') '  CITYP=ORMAS $END'
 else
  write(fid1,'(A)') '  ISPHER=1 CITYP=ORMAS $END'
 end if
 write(fid1,'(A,I0,A)') ' $SYSTEM MWORDS=',CEILING(DBLE(mem*125)/DBLE(nproc)),&
                        ' $END'
 write(fid1,'(2(A,I0),A)') ' $CIDET NCORE=0 NACT=', nif, ' NELS=', ne, ' $END'
 write(fid1,'(A,7(I0,A))')' $ORMAS NSPACE=3 MSTART(1)=1,',ndb0+1,',',&
  ndb0+nacto+1,' MINE(1)=',ndb0*2-order,',',nacte-order,',0 MAXE(1)=',ndb0*2,&
  ',',nacte+order,',',order,' $END'

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(2:7) == '$GUESS') exit
 end do

 BACKSPACE(fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_mrci_gms_inp

! calculate the Davidson size-extensivity correction energy for OpenMolcas
! davidson_e = E_corr*(1-c^2), where c is the reference weight
subroutine calc_davidson_corr_from_out(mrcisd_prog, outname, E_corr, davidson_e)
 use mol, only: ndb, nif, npair, npair0, nacto
 implicit none
 integer :: i, j, k, idx1, idx2, fid1, fid2, system
 integer, allocatable :: det_occ(:,:)
 ! occupation number of each det, 2/1/0/-1
 real(kind=8) :: c
 real(kind=8), intent(in) :: E_corr
 real(kind=8), intent(out) :: davidson_e
 real(kind=8), allocatable :: ci(:)
 character(len=:), allocatable :: str
 character(len=10), intent(in) :: mrcisd_prog
 character(len=240) :: buf, h5name, pyname, pyout
 character(len=240), intent(in) :: outname

 c = 0d0; davidson_e = 0d0 ! initialization
 idx1 = ndb + npair - npair0
 idx2 = idx1 + nacto

 if(TRIM(mrcisd_prog) == 'psi4') then
  open(newunit=fid1,file=TRIM(outname),status='old',position='rewind')
  do while(.true.)
   read(fid1,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(13:22) == 'most impor') exit
  end do ! for while
  if(i /= 0) then
   close(fid1)
   write(6,'(A)') "ERROR in subroutine calc_davidson_corr_from_out: no 'most i&
                  &mpor' found in"
   write(6,'(A)') 'file '//TRIM(outname)
   stop
  end if
  read(buf(7:),*) k ! number of important determinants
  allocate(det_occ(nif,k), ci(k))
  read(fid1,'(A)') buf
  do i = 1, k, 1
   read(fid1,'(A)') buf
   read(buf(10:),*) ci(i)
   buf = buf(37:)
   call det_str2occ_psi4(buf, nif, det_occ(:,i))
  end do ! for i
  close(fid1)
  do i = 1, k, 1
   if(ALL(det_occ(1:idx1,i)==2) .and. ALL(det_occ(idx2+1:nif,i)==0)) c=c+ci(i)*ci(i)
  end do ! for i
  deallocate(det_occ,ci)
 else ! openmolcas
  str = REPEAT(' ',nif)
  i = INDEX(outname, '.', back=.true.)
  h5name = outname(1:i-1)//'.rasscf.h5'
  pyname = outname(1:i-1)//'_ci.py'
  pyout = outname(1:i-1)//'_ci.o'
  open(newunit=fid1,file=TRIM(outname),status='old',position='rewind')
  do while(.true.)
   read(fid1,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(7:14) == 'conf/sym') exit
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine calc_davidson_corr: no 'conf/sym' fou&
                    &nd in "//TRIM(outname)
   close(fid1)
   stop
  end if

  open(newunit=fid2,file=TRIM(pyname),status='replace')
  write(fid2,'(A)') 'import h5py'
  write(fid2,'(A)') "f = h5py.File('"//TRIM(h5name)//"','r')"
  write(fid2,'(A)') "CI_coeff = f['CI_VECTORS']"
  write(fid2,'(A)',advance='no') 'idx = ['
  k = 0 
  do while(.true.)
   read(fid1,*,iostat=j) i, str
   if(j /= 0) exit
   if(str(1:idx1)==REPEAT('2',idx1) .and. str(idx2+1:nif)==REPEAT('0',nif-idx2)) then
    k = k + 1
    write(fid2,'(I0,A1)',advance='no') i, ','
    if(MOD(k,20) == 0) write(fid2,'(A1,/)',advance='no') '\'
   end if
  end do ! for while

  close(fid1)
  write(fid2,'(A)') ']'
  write(fid2,'(A)') 'ref_weight = 0e0'
  write(fid2,'(A)') 'for i in range(0,len(idx)):'
  write(fid2,'(A)') '  rtmp = CI_coeff[0][idx[i]-1]'
  write(fid2,'(A)') '  ref_weight = ref_weight + rtmp*rtmp'
  write(fid2,'(A)') 'f.close()'
  write(fid2,'(A)') "print('%.12e' %ref_weight)"
  close(fid2)
  i = SYSTEM('python '//TRIM(pyname)//" >"//TRIM(pyout)//" 2>&1")
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine calc_davidson_corr: failed to run '&
                   //TRIM(pyname)
   stop
  end if
  open(newunit=fid2,file=TRIM(pyout),status='old')
  read(fid2,*) c
  close(fid2, status='delete')
  open(newunit=fid2,file=TRIM(pyname),status='old')
  close(fid2, status='delete')
 end if

 davidson_e = E_corr*(1d0 - c)
end subroutine calc_davidson_corr_from_out

! convert determinant string into occupation number
subroutine det_str2occ_psi4(buf, nif, det_occ)
 implicit none
 integer :: i, k, iorb, ncol
 integer, intent(in) :: nif
 integer, intent(out) :: det_occ(nif)
 integer, external :: detect_ncol_in_buf
 character(len=5), allocatable :: str(:)
 character(len=240), intent(in) :: buf

 det_occ = 0
 ncol = detect_ncol_in_buf(buf)
 allocate(str(ncol))
 forall(i = 1:ncol) str(i) = REPEAT(' ',5)
 read(buf,*) (str(i),i=1,ncol)

 do i = 1, ncol, 1
  k = INDEX(str(i),'A')
  read(str(i)(1:k-1),*) iorb
  select case(str(i)(k+1:k+1))
  case('A')
   det_occ(iorb) = 1
  case('B')
   det_occ(iorb) = -1
  case('X')
   det_occ(iorb) = 2
  case default
   write(6,'(A)') 'ERROR in subroutine det_str2occ_psi4: invalid occupation&
                    & type='//str(i)(k+1:k+1)
   write(6,'(A,I0)') 'i=', i
   stop
  end select
 end do ! for i
end subroutine det_str2occ_psi4

subroutine prt_mrci_orb_type(order, CIonly)
 implicit none
 integer, intent(in) :: order
 character(len=8) :: mrci
 logical, intent(in) :: CIonly

 select case(order)
 case(2)
  mrci = 'MRCISD'
 case(3)
  mrci = 'MRCISDT'
 case(4)
  mrci = 'MRCISDTQ'
 case default
  write(6,'(/,A)') 'ERROR in subroutine prt_mrci_orb_type: order out of range.'
  write(6,'(A,I0)') 'order=', order
  stop
 end select

 if(CIonly) then
  write(6,'(A)') 'Warning: the CASSCF orbital optimization is strongly recommended'
  write(6,'(A)') 'to be performed before MRCI, unless it is too time-consuming.'
  write(6,'(A)') TRIM(mrci)//' based on CASCI orbitals.'
 else
  write(6,'(A)') TRIM(mrci)//' based on CASSCF orbitals.'
 end if
end subroutine prt_mrci_orb_type

! do uncontacted MRCISDT using one of Gaussian/OpenMolcas/Dalton
subroutine do_mrcisdt()
 use mr_keyword, only: mem, nproc, eist, mrcisdt, mrcisdt_prog, molcas_omp, &
  molcas_path, gau_path, psi4_path, dalton_mpi, gms_path, gms_scr_path, casnofch,&
  chgname, CIonly, bgchg, CtrType
 use mol , only: npair0, casci_e, casscf_e, davidson_e, mrcisd_e, ptchg_e, nuc_pt_e
 use util_wrapper, only: unfchk, fch2inporb_wrap
 implicit none
 integer :: i, system, RENAME
 real(kind=8) :: e
 character(len=24) :: data_string
 character(len=240) :: string, chkname, inpname, outname, mklname

 if(eist == 1) return ! excited state calculation
 if(.not. mrcisdt) return
 write(6,'(//,A)') 'Enter subroutine do_mrcisdt...'

 if(npair0 > 7) then
  write(6,'(A)') 'ERROR in subroutine do_mrcisdt: reference wfn larger than (14,14).'
  write(6,'(A)') 'DMRG-MRCISDT is not supported.'
  stop
 end if

 call prt_mrci_orb_type(3, CIonly)
 write(6,'(A)') 'Frozen_Core = F, MRCISDT computation using program '//&
                 TRIM(mrcisdt_prog)

 select case(TRIM(mrcisdt_prog))
 case('openmolcas')
  call check_exe_exist(molcas_path)
  call fch2inporb_wrap(casnofch, .false., inpname)
  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISDT.input'
  outname = casnofch(1:i)//'MRCISDT.out'
  call prt_mrci_molcas_inp(3, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)

 case('gaussian')
  call check_exe_exist(gau_path)
  i = INDEX(casnofch, '_NO', back=.true.)
  chkname = casnofch(1:i)//'MRCISDT.chk'
  inpname = casnofch(1:i)//'MRCISDT.gjf'
  outname = casnofch(1:i)//'MRCISDT.log'
  call unfchk(casnofch, chkname)
  call prt_mrcisd_gau_inp(3, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_gau_job(gau_path, inpname, .true.)

 case('psi4')
  call check_exe_exist(psi4_path)
  i = SYSTEM('fch2psi '//TRIM(casnofch))
  i = INDEX(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'_psi.inp'
  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISDT.inp'
  outname = casnofch(1:i)//'MRCISDT.out'
  i = RENAME(TRIM(string), TRIM(inpname))
  call prt_mrci_psi4_inp(3, inpname)
  call submit_psi4_job(psi4_path, inpname, nproc)

 case('dalton')
  i = SYSTEM('fch2dal '//TRIM(casnofch))
  i = INDEX(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'.dal'
  chkname = casnofch(1:i-1)//'.mol'
  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISDT.dal'
  mklname = casnofch(1:i)//'MRCISDT.mol'
  outname = casnofch(1:i)//'MRCISDT.out'
  i = RENAME(TRIM(string), TRIM(inpname))
  i = RENAME(TRIM(chkname), TRIM(mklname))
  call prt_mrci_dalton_inp(3, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  i = INDEX(casnofch, '_NO', back=.true.)
  mklname = casnofch(1:i)//'MRCISDT'
  chkname = casnofch(1:i)//'MRCISDT.sout'
  call submit_dalton_job(mklname,mem,nproc,dalton_mpi,.false.,.false.,.false.)

 case('gamess')
  i = SYSTEM('fch2inp '//TRIM(casnofch))
  i = INDEX(casnofch, '.fch', back=.true.)
  string = casnofch(1:i-1)//'.inp'
  i = INDEX(casnofch, '_NO', back=.true.)
  inpname = casnofch(1:i)//'MRCISDT.inp'
  mklname = casnofch(1:i)//'MRCISDT.dat'
  outname = casnofch(1:i)//'MRCISDT.gms'
  i = RENAME(TRIM(string), TRIM(inpname))
  call prt_mrci_gms_inp(3, inpname)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_gms_job(gms_path, gms_scr_path, inpname, nproc)

 case default
  write(6,'(/,A)') 'ERROR in subroutine do_mrcisdt: invalid program='//&
                   TRIM(mrcisdt_prog)
  stop
 end select

 ! read MRCISDT energy from OpenMolcas/Gaussian/Dalton output file
 call read_mrci_energy_from_output(CtrType, mrcisdt_prog, outname, ptchg_e,&
      nuc_pt_e, davidson_e, e)

 mrcisd_e = e ! E(MRCISDT)
 if(CIonly) then   ! E(MRCISDT) - (E(CASCI) or E(CASSCF))
  e = e - casci_e
 else
  e = e - casscf_e
 end if

 write(6,'(/,A,F18.8,1X,A4)') 'E_corr(MRCISDT) =', e, 'a.u.'
 write(6,'(A,F18.8,1X,A4)') 'E(MRCISDT)      =', mrcisd_e, 'a.u.'

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_mrcisdt at '//TRIM(data_string)
end subroutine do_mrcisdt

