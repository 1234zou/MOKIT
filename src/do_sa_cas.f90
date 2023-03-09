! written by jxzou at 20220401: generate the SA-CASSCF input file

! perform SA-CASSCF calculations
subroutine do_sa_cas()
 use mol, only: nif, nbf, ndb, nopen, nacta, nactb, nacto, nacte, npair, &
  npair0, sa_cas_e, ci_mult
 use mr_keyword, only: ist, nacto_wish, nacte_wish, hf_fch, casscf, bgchg, &
  casscf_prog, dmrgscf_prog, nevpt2_prog, chgname, excited, nstate, nevpt2,&
  caspt2, noQD, on_thres, orca_path
 use phys_cons, only: au2ev
 implicit none
 integer :: i, system
 real(kind=8) :: unpaired_e ! unpaired electrons
 real(kind=8), allocatable :: e_ev(:), nevpt2_e(:)
 character(len=10) :: cas_prog = ' '
 character(len=24) :: data_string = ' '
 character(len=240) :: inpname, outname

 if(.not. excited) return

 if(ist == 5) then
  write(6,'(A)') 'Radical index for input NOs:'
  call calc_unpaired_from_fch(hf_fch, 3, .false., unpaired_e)
  ! read nbf, nif, nopen, nacto, ... variables from NO .fch(k) file
  call read_no_info_from_fch(hf_fch, on_thres, nbf, nif, ndb, nopen, nacta, &
                             nactb, nacto, nacte)
  npair0 = nactb; npair = npair0

  ! if the user has specified the active space size, overwrite
  if(nacto_wish>0 .and. nacte_wish>0) then
   nacto = nacto_wish
   nacte = nacte_wish
   nacta = (nacte + nopen)/2
   nactb = (nacte - nopen)/2
  end if
 end if

 if(nacto*nacte == 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_sa_cas: nacto*nacte is 0.'
  write(6,'(A)') 'Some variables have not been initialized.'
  stop
 end if
 if(nacto /= nacte) then
  write(6,'(A)') 'ERROR in subroutine do_sa_cas: nacto/=nacte not supported.'
  stop
 end if

 casscf = .true.
 if(nacto > 15) casscf = .false.
 if(casscf) then
  data_string = 'SA-CASSCF'
  cas_prog = casscf_prog
 else
  data_string = 'DMRG-SA-CASSCF'
  cas_prog = dmrgscf_prog
 end if
 write(6,'(/,2(A,I0),A)') TRIM(data_string)//'(', nacte, 'e,', nacto,&
                          'o) using program '//TRIM(cas_prog)

 i = index(hf_fch, '.fch', back=.true.)
 select case(TRIM(cas_prog))
 case('pyscf')
  inpname = hf_fch(1:i-1)//'_SA.py'
  outname = hf_fch(1:i-1)//'_SA.out'
  call prt_sacas_script_into_py(inpname, hf_fch)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_pyscf_job(inpname)
 case('gaussian')
  inpname = hf_fch(1:i-1)//'_SA.gjf'
  outname = hf_fch(1:i-1)//'_SA.log'
  call prt_sacas_gjf(inpname, hf_fch)
 case('orca')
  inpname = hf_fch(1:i-1)//'_SA.inp'
  outname = hf_fch(1:i-1)//'_SA.out'
  call prt_sacas_orca_inp(inpname, hf_fch)
  call submit_orca_job(orca_path, inpname)
 case('gamess')
  inpname = hf_fch(1:i-1)//'_SA.inp'
  outname = hf_fch(1:i-1)//'_SA.gms'
  call prt_sacas_gms_inp(inpname, hf_fch)
 case default
  write(6,'(A)') 'ERROR in subroutine do_sa_cas: CASSCF_prog='//TRIM(cas_prog)&
                 //' unrecognized or unsupported.'
  stop
 end select

 allocate(sa_cas_e(0:nstate), ci_mult(0:nstate)) ! 0 means ground state
 call read_sa_cas_energies_from_output(cas_prog,outname,nstate,sa_cas_e,ci_mult)

 allocate(e_ev(nstate), source=0d0)
 forall(i = 1:nstate) e_ev(i) = (sa_cas_e(i) - sa_cas_e(0))*au2ev
 write(6,'(A)') 'Multi-root energies in SA-CASSCF(0 for ground state):'
 write(6,'(A,A)') REPEAT(' ',52),'E_ex/eV'
 write(6,'(A,I3,A,F16.8,A,F6.3)') 'State ',0,', E =',sa_cas_e(0),' a.u., <S**2> =',&
                                   ci_mult(0)
 do i = 1, nstate, 1
  write(6,'(A,I3,A,F16.8,A,F6.3,2X,F6.2)') 'State ',i,', E =',sa_cas_e(i),&
                                     ' a.u., <S**2> =', ci_mult(i), e_ev(i)
 end do ! for i
 deallocate(e_ev)

 if(nevpt2) then
  if(TRIM(nevpt2_prog) /= TRIM(cas_prog)) then
   write(6,'(A)') 'ERROR in subroutine do_sa_cas: currently only NEVPT2_prog=&
                  &CASSCF_prog is allowed.'
   write(6,'(A)') 'But NEVPT2_prog='//TRIM(nevpt2_prog)//' and CASSCF_prog='//&
                   TRIM(cas_prog)
   stop
  end if

  allocate(nevpt2_e(0:nstate))
  select case(TRIM(nevpt2_prog))
  case('pyscf')
   call read_multiroot_nevpt2_from_pyscf(outname, nstate, nevpt2_e)
  case('orca')
   call read_multiroot_nevpt2_from_orca(outname, nstate, nevpt2_e)
  case default
   stop
  end select
  deallocate(nevpt2_e)
 end if
end subroutine do_sa_cas

! print (DMRG-)SA-CASSCF script into a given .py file
subroutine prt_sacas_script_into_py(pyname, gvb_fch)
 use mol, only: nacto, nacta, nactb
 use mr_keyword, only: mem, nproc, casscf, dmrgscf, maxM, hardwfn, crazywfn, &
  dkh2_or_x2c, RI, RIJK_bas, hf_fch, mixed_spin, nstate, nevpt2
 implicit none
 integer :: i, fid1, fid2, system, RENAME
 real(kind=8) :: ss ! spin square S(S+1)
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, pyname1, cmofch
 character(len=240), intent(in) :: pyname, gvb_fch

 i = system('bas_fch2py '//TRIM(gvb_fch))
 i = index(gvb_fch, '.fch', back=.true.)
 pyname1 = gvb_fch(1:i-1)//'.py'
 i = RENAME(TRIM(pyname1), TRIM(pyname))

 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
 pyname1 = TRIM(pyname)//'.t'
 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do

 if(dmrgscf) then
  write(fid2,'(A)') 'from pyscf import mcscf, dmrgscf, lib'
 else
  write(fid2,'(A)') 'from pyscf import mcscf, lib'
 end if
 if(nevpt2) write(fid2,'(A)') 'from pyscf import mrpt'
 write(fid2,'(A)') 'from mokit.lib.py2fch import py2fch'
 write(fid2,'(A)') 'from shutil import copyfile'
 write(fid2,'(A,/)') 'import numpy as np'
 if(dmrgscf) then
  write(fid2,'(A,I0,A)') "dmrgscf.settings.MPIPREFIX ='mpirun -n ",nproc,"'"
 end if
 write(fid2,'(A,I0,A1,/)') 'lib.num_threads(',nproc,')'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:9) == 'mf.kernel') exit
  write(fid2,'(A)') TRIM(buf)
 end do
 write(fid2,'(A,I0,A)') 'mf.max_memory = ', mem*1000, ' # MB'
 write(fid2,'(A)') 'mf.kernel()'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')

 ! mem*500 is in fact mem*1000/2. The mc.max_memory and fcisolver.max_memory seem
 ! not to share common memory, they used two memory, so I have to make them half
 if(casscf) then ! SA-CASSCF
  if(dkh2_or_x2c) then
   write(fid2,'(A)',advance='no') 'mc = mcscf.CASSCF(mf.x2c1e(),'
  else
   write(fid2,'(A)',advance='no') 'mc = mcscf.CASSCF(mf,'
  end if
  write(fid2,'(3(I0,A))',advance='no') nacto,',(',nacta,',',nactb,')'
  if(RI) then
   write(fid2,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
  else
   write(fid2,'(A)') ')'
  end if
  write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ',mem*200,' # MB'
 else ! DMRG-SA-CASSCF
  if(dkh2_or_x2c) then
   write(fid2,'(A)',advance='no') 'mc = dmrgscf.DMRGSCF(mf.x2c1e(),'
  else
   write(fid2,'(A)',advance='no') 'mc = dmrgscf.DMRGSCF(mf,'
  end if
  write(fid2,'(3(I0,A))') nacto,',(',nacta,',',nactb,'))'
  write(fid2,'(A,I0)') 'mc.fcisolver.maxM = ', maxM
  write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ',CEILING(DBLE(mem)/DBLE((5*nproc))),' # GB'
 end if

 write(fid2,'(A)',advance='no') 'mc = mc.state_average_(['
 write(fid2,'(5(A,I0,A))',advance='no') ('1e0/',nstate+1,'e0,',i=1,nstate+1,1)
 if(hardwfn) then
  write(fid2,'(A)') '0.0,0.0,0.0])'
 else if(crazywfn) then
  write(fid2,'(A)') '0.0,0.0,0.0,0.0,0.0,0.0])'
 else
  write(fid2,'(A)') '])'
 end if
 write(fid2,'(A)') 'mc.conv_tol = 1e-9'
 write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*800, ' # MB'
 write(fid2,'(A)') 'mc.max_cycle = 200'
 write(fid2,'(A,I0)') 'mc.fcisolver.spin = ', nacta-nactb

 call prt_hard_or_crazy_casci_pyscf(fid2, nacta-nactb, hardwfn,crazywfn,.false.)
 ss = DBLE(nacta - nactb)*0.5d0
 ss = ss*(ss+1d0)
 if(.not. mixed_spin) write(fid2,'(A,F7.3,A)') 'mc.fix_spin_(ss=',ss,')'
 write(fid2,'(A)') 'mc.verbose = 5'
 write(fid2,'(A)') 'mc.kernel()'
 if(nevpt2) write(fid2,'(A)') 'mo = mc.mo_coeff.copy()'

 i = index(hf_fch, '.fch', back=.true.)
 cmofch = hf_fch(1:i-1)//'_SA.fch'
 write(fid2,'(/,A)') '# save MOs into .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(gvb_fch)//"', '"//TRIM(cmofch)//"')"
 write(fid2,'(A)') 'noon = np.zeros(nif)'
 write(fid2,'(A)') "py2fch('"//TRIM(cmofch)//"',nbf,nif,mc.mo_coeff,'a',noon,False)"
 ! mc.mo_occ only exists for PySCF >= 1.7.4

 if(nevpt2) then
  write(fid2,'(/,A)') '# perform multi-root CASCI'
  if(dkh2_or_x2c) then
   write(fid2,'(A)',advance='no') 'mc = mcscf.CASCI(mf.x2c1e(),'
  else
   write(fid2,'(A)',advance='no') 'mc = mcscf.CASCI(mf,'
  end if
  write(fid2,'(3(I0,A))',advance='no') nacto,',(',nacta,',',nactb,')'

  if(casscf) then ! CASSCF
   if(RI) then
    write(fid2,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
   else
    write(fid2,'(A)') ')'
   end if
   write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ',mem*200,' # MB'
  else            ! DMRG-SA-CASSCF
   write(fid2,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=',maxM,')'
   write(fid2,'(A,I0,A)') 'mc.fcisolver.memory = ',CEILING(DBLE(mem)/DBLE((5*nproc))),' # GB'
  end if

  write(fid2,'(A)',advance='no') 'mc.fcisolver.nroots = '
  if(hardwfn) then
   write(fid2,'(I0)') nstate+4
  else if(crazywfn) then
   write(fid2,'(I0)') nstate+7
  else
   write(fid2,'(I0)') nstate+1
  end if
  call prt_hard_or_crazy_casci_pyscf(fid2, nacta-nactb,hardwfn,crazywfn,.false.)
  if(.not. mixed_spin) write(fid2,'(A,F7.3,A)') 'mc.fix_spin_(ss=',ss,')'
  write(fid2,'(A)') 'mc.verbose = 4'
  write(fid2,'(A)') 'mc.kernel(mo)'

  write(fid2,'(/,A)') '# NEVPT2 based on multi-root CASCI'
  do i = 0, nstate, 1
   write(fid2,'(A,I0,A)') 'mrpt.NEVPT(mc, root=',i,').kernel()'
  end do ! for i
 end if
 close(fid2)

 i = RENAME(TRIM(pyname1), TRIM(pyname))
end subroutine prt_sacas_script_into_py

! print SA-CASSCF input file of Gaussian
subroutine prt_sacas_gjf(gjfname, hf_fch)
 use mol, only: nacto, nacte
 use mr_keyword, only: mem, nproc, mixed_spin, nstate
 use util_wrapper, only: unfchk
 implicit none
 integer :: i, fid
 real(kind=8), allocatable :: weight(:)
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname, hf_fch

 i = index(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'
 call unfchk(hf_fch, chkname)

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(2(A,I0))') '#p CASSCF(',nacte,',',nacto
 if(mixed_spin) write(fid,'(A)',advance='no') ',SlaterDet'
 write(fid,'(A)') ',StateAverage) chkbasis nosymm int=nobasistransform'
 write(fid,'(A,/)') 'scf(maxcycle=500) guess=read geom=allcheck'

 allocate(weight(nstate))
 weight = 1d0/DBLE(nstate)
 do i = 1, nstate, 1
  write(fid,'(F10.8)') weight(i)
 end do ! for i

 deallocate(weight)
 write(fid,'(A)',advance='no')
 close(fid)
end subroutine prt_sacas_gjf

! print SA-CASSCF input file of ORCA
subroutine prt_sacas_orca_inp(inpname, hf_fch)
 use util_wrapper, only: fch2mkl_wrap, mkl2gbw
 use mol, only: nacto, nacte, mult
 use mr_keyword, only: mem, nproc, dkh2_or_x2c, nevpt2, FIC, DLPNO, F12, RI, &
  RIJK_bas, mixed_spin, nstate, hardwfn, crazywfn
 implicit none
 integer :: i, fid, fid1
 character(len=240) :: buf, mklname, gbwname, inpname1
 character(len=240), intent(in) :: inpname, hf_fch

 call fch2mkl_wrap(hf_fch)
 i = index(hf_fch, '.fch', back=.true.)
 inpname1 = hf_fch(1:i-1)//'_o.inp'
 mklname = hf_fch(1:i-1)//'_o.mkl'

 i = index(inpname, '.inp', back=.true.)
 gbwname = inpname(1:i-1)//'.gbw'
 call mkl2gbw(mklname, gbwname)
 call delete_file(mklname)

 open(newunit=fid,file=TRIM(inpname1),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname),status='replace')
 write(fid1,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid1,'(A,I0)') '%maxcore ', CEILING(1000d0*DBLE(mem)/DBLE(nproc))
 write(fid1,'(A)',advance='no') '! CASSCF'
 ! RIJK in CASSCF must be combined with CONVentional
 if(RI) write(fid1,'(A)',advance='no') ' RIJK conv '//TRIM(RIJK_bas)
 write(fid1,'(A)') ' VeryTightSCF'

 if(dkh2_or_x2c) then
  write(fid1,'(A)') '%rel'
  write(fid1,'(A)') ' method DKH'
  write(fid1,'(A)') ' order 2'
  write(fid1,'(A)') 'end'
 end if
 write(fid1,'(A)') '%casscf'
 write(fid1,'(A,I0)') ' nel ', nacte
 write(fid1,'(A,I0)') ' norb ', nacto
 if(mixed_spin) then
  write(fid1,'(A,I0,A1,I0)') ' mult ',mult,',',mult+2
  if(MOD(nstate,2) == 0) then
   write(fid1,'(A,I0,A1,I0)') ' nroots ',nstate/2+1,',',nstate/2
  else
   write(fid1,'(A,I0,A1,I0)') ' nroots ',(nstate+1)/2,',',(nstate+1)/2
  end if
 else
  write(fid1,'(A,I0)') ' mult ', mult
  write(fid1,'(A,I0)') ' nroots ', nstate+1
 end if
 write(fid1,'(A)') ' maxiter 200'
 if(RI) write(fid1,'(A)') ' TrafoStep RI'
 call prt_hard_or_crazy_casci_orca(fid1, hardwfn, crazywfn)

 if(nevpt2) then
  if(FIC) then
   if(DLPNO) then
    write(fid1,'(A)') ' PTMethod DLPNO_NEVPT2'
   else
    write(fid1,'(A)') ' PTMethod FIC_NEVPT2'
   end if
  else
   write(fid1,'(A)') ' PTMethod SC_NEVPT2'
  end if
  if(F12) then
   write(fid1,'(A)') ' PTSettings'
   write(fid1,'(A)') '  F12 true'
   write(fid1,'(A)') ' end'
  end if
 end if

 write(fid1,'(A)') 'end'
 if(nevpt2) then
  write(fid1,'(A)') '%method'
  write(fid1,'(A)') ' FrozenCore FC_NONE'
  write(fid1,'(A)') 'end'
 end if

 do i = 1, 3   ! skip 3 lines
  read(fid,'(A)') buf
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
end subroutine prt_sacas_orca_inp

! print SA-CASSCF input file of GAMESS
subroutine prt_sacas_gms_inp(inpname, hf_fch)
 use util_wrapper, only: fch2inp_wrap
 use mol, only: ndb, npair, npair0, nacto, nacte, charge, mult
 use mr_keyword, only: mem, nproc, hardwfn, crazywfn, dkh2_or_x2c, cart, &
  mixed_spin, nstate
 implicit none
 integer :: i, ncore, fid1, fid2, RENAME
 real(kind=8), allocatable :: weight(:)
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname, hf_fch

 ncore = ndb + npair - npair0
 call fch2inp_wrap(hf_fch, .false., 0, 0)
 i = index(hf_fch, '.fch', back=.true.)
 inpname1 = hf_fch(1:i-1)//'.inp'
 i = RENAME(TRIM(inpname1), TRIM(inpname))
 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 write(fid2,'(A)',advance='no') ' $CONTRL SCFTYP=MCSCF RUNTYP=ENERGY ICHARG='
 write(fid2,'(2(I0,A))') charge, ' MULT=', mult, ' NOSYM=1'
 write(fid2,'(A)',advance='no') '  ICUT=11'
 if(dkh2_or_x2c) write(fid2,'(A)',advance='no') ' RELWFN=DK'
 if(.not. cart) then
  write(fid2,'(A)') ' ISPHER=1 $END'
 else
  write(fid2,'(A)') ' $END'
 end if

 write(fid2,'(A,I0,A)') ' $SYSTEM MWORDS=',CEILING(DBLE(mem*125)/DBLE(nproc)),&
                        ' $END'
 write(fid2,'(A)',advance='no') ' $DET'
 write(fid2,'(3(A,I0),A)',advance='no') ' NCORE=',ncore,' NELS=',nacte,' NACT=',&
                                        nacto,' ITERMX=500'
 if(mixed_spin) write(fid2,'(A)',advance='no') ' PURES=.FALSE.'
 write(fid2,'(A)',advance='no') ' NSTATE='
 if(hardwfn) then
  write(fid2,'(I0)') nstate+6
 else if(crazywfn) then
  write(fid2,'(I0)') nstate+3
 else
  write(fid2,'(I0)') nstate
 end if

 write(fid2,'(A)',advance='no') '  WSTATE(1)='
 allocate(weight(nstate+1))
 weight = 1d0/DBLE(nstate+1)
 do i = 1, nstate, 1
  write(fid2,'(F10.8,A1)',advance='no') weight(i),','
 end do ! for i
 write(fid2,'(F10.8,A1)',advance='no') weight(nstate+1)
 deallocate(weight)

 write(fid2,'(A)') ' $END'
 write(fid2,'(A)') ' $MCSCF MAXIT=200 $END'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(2:7) == '$GUESS') exit
 end do

 BACKSPACE(fid1)
 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_sacas_gms_inp

subroutine submit_pyscf_job(pyname)
 implicit none
 integer :: i, system
 character(len=240) :: outname
 character(len=480) :: buf
 character(len=240), intent(in) :: pyname

 i = index(pyname, '.py', back=.true.)
 outname = pyname(1:i-1)//'.out'

 write(buf,'(A)') 'python '//TRIM(pyname)//' >'//TRIM(outname)//" 2>&1"
 write(6,'(A)') '$'//TRIM(buf)
 i = system(TRIM(buf))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subrouitine submit_pyscf_job: PySCF job failed.'
  write(6,'(A)') 'Please open file '//TRIM(outname)//' and check.'
  stop
 end if
end subroutine submit_pyscf_job

subroutine read_sa_cas_energies_from_output(cas_prog, outname, nstate, &
                                            sa_cas_e, ci_mult)
 use mol, only: mult
 implicit none
 integer, intent(in) :: nstate ! ground state included
 character(len=10), intent(in) :: cas_prog
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: sa_cas_e(0:nstate), ci_mult(0:nstate)

 select case(TRIM(cas_prog))
 case('gaussian')
!  call read_sa_cas_energies_from_gau_log(outname, nstate, sa_cas_e)
  ci_mult = mult
 case('pyscf')
  call read_sa_cas_energies_from_pyscf_out(outname, nstate, sa_cas_e, ci_mult)
 case('orca')
  call read_sa_cas_energies_from_orca_out(outname, nstate, sa_cas_e, ci_mult)
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_sa_cas_energies_from_output: &
                   &cas_prog cannot be recognized.'
  write(6,'(A)') 'cas_prog='//TRIM(cas_prog)
  stop
 end select
end subroutine read_sa_cas_energies_from_output

subroutine read_sa_cas_energies_from_pyscf_out(outname, nstate, sa_cas_e, ci_mult)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: sa_cas_e(0:nstate), ci_mult(0:nstate)

 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:16) == 'CASCI energy for') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_sa_cas_energies_from_pyscf_out:&
                  & no 'CASCI energy for'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  stop
 end if

 do i = 0, nstate, 1
  read(fid,'(A)') buf
  j = index(buf, '=')
  read(buf(j+1:),*) sa_cas_e(i)
  buf(j:j) = ' '
  j = index(buf, '=')
  read(buf(j+1:),*) ci_mult(i)
 end do ! for i

 close(fid)
end subroutine read_sa_cas_energies_from_pyscf_out

subroutine read_sa_cas_energies_from_orca_out(outname, nstate, sa_cas_e, ci_mult)
 implicit none
 integer :: i, j, k, m, mult, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: s
 real(kind=8), intent(out) :: sa_cas_e(0:nstate), ci_mult(0:nstate)

 sa_cas_e = 0d0; ci_mult = 0d0

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:24) == 'CAS-SCF STATES FOR BLOCK') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine read_sa_cas_energies_from_orca_out: no&
                  & 'CAS-SCF STATES FOR BLOCK'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  stop
 end if

 ! CAS-SCF STATES FOR BLOCK  1 MULT= 3 NROOTS= 3
 i = index(buf, '=', back=.true.)
 read(buf(i+1:),*) k

 j = index(buf, '=')
 read(buf(j+1:),*) mult
 s = 0.5d0*DBLE(mult-1)
 ci_mult(0:k-1) = s*(s+1d0)

 m = -1
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == 'ROOT') then
   m = m + 1
   i = index(buf, '=')
   read(buf(i+1:),*) sa_cas_e(m)
   if(m+1 == k) exit
  end if
 end do ! for while

 if(nstate+1 == k) then
  close(fid)
  return
 end if

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:24) == 'CAS-SCF STATES FOR BLOCK') exit
 end do ! for while

 i = index(buf, '=', back=.true.)
 read(buf(i+1:),*) k

 j = index(buf, '=')
 read(buf(j+1:),*) mult
 s = 0.5d0*DBLE(mult-1)
 ci_mult(m+1:) = s*(s+1d0)

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == 'ROOT') then
   m = m + 1
   i = index(buf, '=')
   read(buf(i+1:),*) sa_cas_e(m)
   if(m == nstate) exit
  end if
 end do ! for while

 close(fid)

 ! sort the electronic states
 do i = 0, nstate-1, 1
  do j = i+1, nstate, 1
   if(sa_cas_e(i) > sa_cas_e(j)) then
    s = sa_cas_e(i); sa_cas_e(i) = sa_cas_e(j); sa_cas_e(j) = s
    s = ci_mult(i); ci_mult(i) = ci_mult(j); ci_mult(j) = s
   end if
  end do ! for i
 end do ! for j
end subroutine read_sa_cas_energies_from_orca_out

! read multi-root CASCI-based NEVPT2 energies from a PySCF out file
subroutine read_multiroot_nevpt2_from_pyscf(outname, nstate, nevpt2_e)
 use phys_cons, only: au2ev
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), allocatable :: casci_e(:)
 real(kind=8), intent(out) :: nevpt2_e(0:nstate)
 logical :: no_conv

 no_conv = .false.
 allocate(casci_e(0:nstate), source=0d0)
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:10) == 'CASCI conv') exit
  if(buf(1:14) == 'CASCI not conv') then
   no_conv = .true.
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_multiroot_nevpt2_from_pyscf: no &
                  &'CASCI conv' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  stop
 end if

 if(no_conv) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning in subroutine read_multiroot_nevpt2_from_pyscf: CASCI&
                & not converged.'
  write(6,'(A)') 'You should be cautious about the CASCI and NEVPT2 energies.&
                & You are recommended'
  write(6,'(A)') 'to open file '//TRIM(outname)//' and check.'
  write(6,'(A)') 'The program will not stop but continue...'
  write(6,'(A)') REPEAT('-',79)
 end if

 do i = 0, nstate, 1
  read(fid,'(A)') buf
  j = index(buf,'=')
  read(buf(j+1:),*) casci_e(i)
 end do ! for i

 i = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:8) == 'Nevpt2 E') then
   j = index(buf,'=')
   read(buf(j+1:),*) nevpt2_e(i)
   i = i + 1
   if(i == nstate) exit
  end if
 end do ! for while

 close(fid)
 nevpt2_e = nevpt2_e + casci_e

 forall(i = 1:nstate) casci_e(i) = (nevpt2_e(i)-nevpt2_e(0))*au2ev
 write(6,'(/,A)') 'NEVPT2 energies:                    E_ex/eV'
 write(6,'(A,I3,A,F16.8,A)') 'State ',0,', E =',nevpt2_e(0),' a.u.'
 do i = 1, nstate, 1
  write(6,'(A,I3,A,F16.8,A,2X,F6.2)') 'State ',i,', E =',nevpt2_e(i),' a.u.',&
                                      casci_e(i)
 end do ! for i

 deallocate(casci_e)
end subroutine read_multiroot_nevpt2_from_pyscf

! read multi-root CASCI-based NEVPT2 energies from an ORCA out file
subroutine read_multiroot_nevpt2_from_orca(outname, nstate, nevpt2_e)
 use phys_cons, only: au2ev
 integer :: fid
 integer, intent(in) :: nstate
 !character(len=240) :: buf
 character(len=240), intent(in) :: outname
 !real(kind=8), allocatable :: casci_e(:)
 real(kind=8), intent(out) :: nevpt2_e(0:nstate)

 nevpt2_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 close(fid)
end subroutine read_multiroot_nevpt2_from_orca

