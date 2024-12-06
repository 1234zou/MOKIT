! written by jxzou at 20220401: generate the SA-CASSCF input file

! perform SA-CASSCF calculations
subroutine do_sa_cas()
 use mol, only: nif, nbf, ndb, nopen, nacta, nactb, nacto, nacte, npair, &
  npair0, sa_cas_e, ci_ssquare, fosc
 use mr_keyword, only: mem, nproc, ist, nacto_wish, nacte_wish, hf_fch, casscf,&
  dmrgscf, bgchg, casscf_prog, dmrgscf_prog, nevpt2_prog, chgname, excited, &
  nstate, nevpt2, on_thres, orca_path, molcas_omp
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
  write(6,'(/,A)') 'ERROR in subroutine do_sa_cas: nacto/=nacte not supported.'
  stop
 end if

 casscf = .true.
 if(nacto > 15) then
  casscf = .false.
  dmrgscf = .true.
 end if

 if(casscf) then
  data_string = 'SA-CASSCF'
  cas_prog = casscf_prog
 else
  data_string = 'SA-DMRG-CASSCF'
  cas_prog = dmrgscf_prog
 end if
 write(6,'(/,2(A,I0),A)') TRIM(data_string)//'(', nacte, 'e,', nacto,&
                          'o) using program '//TRIM(cas_prog)

 i = INDEX(hf_fch, '.fch', back=.true.)
 select case(TRIM(cas_prog))
 case('pyscf')
  inpname = hf_fch(1:i-1)//'_SA-CAS.py'
  outname = hf_fch(1:i-1)//'_SA-CAS.out'
  call prt_sacas_script_into_py(inpname, hf_fch)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_pyscf_job(inpname, .true.)
 case('gaussian')
  inpname = hf_fch(1:i-1)//'_SA-CAS.gjf'
  outname = hf_fch(1:i-1)//'_SA-CAS.log'
  call prt_sacas_gjf(inpname, hf_fch)
 case('orca')
  inpname = hf_fch(1:i-1)//'_SA-CAS.inp'
  outname = hf_fch(1:i-1)//'_SA-CAS.out'
  call prt_sacas_orca_inp(inpname, hf_fch)
  call submit_orca_job(orca_path, inpname, .true., .false., .false.)
 case('gamess')
  inpname = hf_fch(1:i-1)//'_SA-CAS.inp'
  outname = hf_fch(1:i-1)//'_SA-CAS.gms'
  call prt_sacas_gms_inp(inpname, hf_fch)
 case('openmolcas')
  inpname = hf_fch(1:i-1)//'_SA-CAS.input'
  outname = hf_fch(1:i-1)//'_SA-CAS.out'
  call prt_sacas_molcas_inp(inpname, hf_fch)
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)
  call copy_nto_from_orb2fch(hf_fch, nstate)
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_sa_cas: CASSCF_prog='//TRIM(cas_prog)&
                 //' unrecognized or unsupported.'
  stop
 end select

 allocate(sa_cas_e(0:nstate), ci_ssquare(0:nstate), fosc(nstate))
 call read_sa_cas_energies_from_output(cas_prog,outname,nstate,sa_cas_e,ci_ssquare)
 if(dmrgscf) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning: DMRG NTOs are not implemented yet. Oscillator stren&
                 &gths f_osc'
  write(6,'(A)') 'below will be set to zero.'
  write(6,'(A)') REPEAT('-',79)
 else
  call read_fosc_from_output(cas_prog, outname, nstate, fosc)
 end if

 allocate(e_ev(nstate), source=0d0)
 forall(i = 1:nstate) e_ev(i) = (sa_cas_e(i) - sa_cas_e(0))*au2ev
 write(6,'(A)') 'CASCI energies after SA-CASSCF(0 for ground state):   E_ex/eV &
                &  fosc'
 write(6,'(A,I3,A,F16.8,A,F6.3)') 'State ',0,', E =',sa_cas_e(0),' a.u. <S**2> =',&
                                   ci_ssquare(0)
 do i = 1, nstate, 1
  write(6,'(A,I3,A,F16.8,A,F6.3,2X,F7.3,2X,F7.4)') 'State ', i, ', E =', &
               sa_cas_e(i), ' a.u. <S**2> =', ci_ssquare(i), e_ev(i), fosc(i)
 end do ! for i

 if(nevpt2) then
  if(TRIM(nevpt2_prog) /= TRIM(cas_prog)) then
   write(6,'(/,A)') 'ERROR in subroutine do_sa_cas: currently only NEVPT2_prog=&
                    &CASSCF_prog is allowed.'
   write(6,'(A)') 'But NEVPT2_prog='//TRIM(nevpt2_prog)//' and CASSCF_prog='//&
                   TRIM(cas_prog)
   stop
  end if

  allocate(nevpt2_e(0:nstate))
  select case(TRIM(nevpt2_prog))
  case('pyscf')
   call read_multiroot_nevpt2_from_pyscf(outname, dmrgscf, nstate, nevpt2_e)
  case('orca')
   call read_multiroot_nevpt2_from_orca(outname, nstate, nevpt2_e)
  case default
   write(6,'(/,A)') 'ERROR in subroutine do_sa_cas: NEVPT2_prog cannot be recog&
                    &nized.'
   write(6,'(A)') 'NEVPT2_prog='//TRIM(nevpt2_prog)
   stop
  end select

  forall(i = 1:nstate) e_ev(i) = (nevpt2_e(i) - nevpt2_e(0))*au2ev
  write(6,'(/,A)') 'NEVPT2 energies:                     E_ex/eV'
  write(6,'(A,I3,A,F16.8,A)') 'State ',0,', E =',nevpt2_e(0),' a.u.'
  do i = 1, nstate, 1
   write(6,'(A,I3,A,F16.8,A,2X,F7.3)') 'State ',i,', E =',nevpt2_e(i),' a.u.',&
                                       e_ev(i)
  end do ! for i
  deallocate(nevpt2_e)
 end if

 deallocate(e_ev)
end subroutine do_sa_cas

! print (DMRG-)SA-CASSCF script into a given .py file
subroutine prt_sacas_script_into_py(pyname, gvb_fch)
 use mol, only: nacto, nacte, nacta, nactb
 use mr_keyword, only: mem, nproc, casscf, dmrgscf, maxM, block_mpi, hardwfn, &
  crazywfn, RI, RIJK_bas, hf_fch, mixed_spin, nstate, nevpt2
 use util_wrapper, only: bas_fch2py_wrap
 implicit none
 integer :: i, nacta1, nactb1, fid1, fid2, RENAME
 real(kind=8) :: ss ! spin square S(S+1)
 character(len=21) :: RIJK_bas1
 character(len=240) :: buf, pyname1, cmofch
 character(len=240), intent(in) :: pyname, gvb_fch

 if(mixed_spin) then
  ! set as the lowest spin in order to obtain different spins in SA-CASSCF
  nactb1 = nacte/2
  nacta1 = nacte - nactb1
 else
  nacta1 = nacta
  nactb1 = nactb
 end if
 call bas_fch2py_wrap(gvb_fch, .false., pyname)
 if(RI) call auxbas_convert(RIJK_bas, RIJK_bas1, 1)
 pyname1 = TRIM(pyname)//'.t'

 open(newunit=fid1,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:31) == 'from pyscf import gto, scf, lib') exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(dmrgscf) then
  if(nevpt2) then
   write(fid2,'(A)') 'from pyscf import gto, scf, mcscf, dmrgscf, mrpt, lib'
  else
   write(fid2,'(A)') 'from pyscf import gto, scf, mcscf, dmrgscf, lib'
  end if
 else
  if(nevpt2) then
   write(fid2,'(A)') 'from pyscf import gto, scf, mcscf, mrpt, lib'
  else
   write(fid2,'(A)') 'from pyscf import gto, scf, mcscf, lib'
  end if
 end if
 write(fid2,'(A)') 'from mokit.lib.py2fch import py2fch'
 write(fid2,'(A)') 'from mokit.lib.excited import gen_nto_and_fosc_from_mo_tdm'

 do while(.true.)
  read(fid1,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 write(fid2,'(A)') 'from shutil import copyfile'
 write(fid2,'(A,/)') 'import numpy as np'

 if(dmrgscf) then
  write(fid2,'(A)',advance='no') "dmrgscf.settings.MPIPREFIX = '"
  if(block_mpi) write(fid2,'(A,I0)',advance='no') 'mpirun -n ', nproc
  write(fid2,'(A)') "'"
 end if
 write(fid2,'(A,I0)') 'nproc = ', nproc
 write(fid2,'(A)') 'lib.num_threads(nproc)'

 do while(.true.)
  read(fid1,'(A)') buf
  if(buf(1:15) == 'lib.num_threads') cycle
  if(buf(1:13) == 'mf.max_memory') cycle
  if(buf(1:9) == 'mf.kernel') exit
  write(fid2,'(A)') TRIM(buf)
 end do
 write(fid2,'(A,I0,A)') 'mf.max_memory = ', mem*1000, ' # MB'
 write(fid2,'(A)') TRIM(buf)

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do
 close(fid1,status='delete')

 ! mem*500 is in fact mem*1000/2. The mc.max_memory and fcisolver.max_memory seem
 ! not to share common memory, they used two memory, so I have to make them half
 if(casscf) then ! SA-CASSCF
  write(fid2,'(3(A,I0),A)',advance='no') 'mc = mcscf.CASSCF(mf,', nacto, ',(', &
                                          nacta1, ',', nactb1, ')'
  if(RI) then
   write(fid2,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
  else
   write(fid2,'(A)') ')'
  end if
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*700, ' # MB'
  write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ',mem*300,' # MB'
 else ! DMRG-SA-CASSCF
  write(fid2,'(3(A,I0),A)') 'mc = dmrgscf.DMRGSCF(mf,', nacto, ',(', nacta1, &
                            ',', nactb1, '))'
  write(fid2,'(A,I0)') 'mc.fcisolver.maxM = ', maxM
  call prt_block_mem(fid2, mem, nproc, block_mpi)
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
 write(fid2,'(A)') 'mc.max_cycle = 200'
 if(.not. (mixed_spin .or. dmrgscf)) then
  write(fid2,'(A,I0)') 'mc.fcisolver.spin = ', nacta-nactb
 end if

 if(.not. dmrgscf) then
  call prt_hard_or_crazy_casci_pyscf(fid2,nacta-nactb,hardwfn,crazywfn,.false.)
  ss = DBLE(nacta - nactb)*0.5d0
  ss = ss*(ss + 1d0)
  if(.not. mixed_spin) write(fid2,'(A,F7.3,A)') 'mc.fix_spin_(ss=', ss, ')'
 end if
 write(fid2,'(A)') 'mc.verbose = 4'
 write(fid2,'(A)') 'mc.kernel()'
 write(fid2,'(A)') 'mo = mc.mo_coeff.copy() # backup'

 i = INDEX(hf_fch, '.fch', back=.true.)
 cmofch = hf_fch(1:i-1)//'_SA-CAS.fch'
 write(fid2,'(/,A)') '# save MOs into .fch file'
 write(fid2,'(A)') "copyfile('"//TRIM(gvb_fch)//"','"//TRIM(cmofch)//"')"
 write(fid2,'(A)') 'noon = np.zeros(nif)'
 write(fid2,'(A)') "py2fch('"//TRIM(cmofch)//"',nbf,nif,mo,'a',noon,False,False)"
 ! Note: mc.mo_occ only exists for PySCF >= 1.7.4

 write(fid2,'(/,A)') '# perform multi-root CASCI'
 write(fid2,'(3(A,I0),A)',advance='no') 'mc = mcscf.CASCI(mf,', nacto, ',(', &
                                         nacta1, ',', nactb1, ')'

 if(casscf) then ! multi-root CASCI
  if(RI) then
   write(fid2,'(A)') ").density_fit(auxbasis='"//TRIM(RIJK_bas1)//"')"
  else
   write(fid2,'(A)') ')'
  end if
  write(fid2,'(A,I0,A)') 'mc.max_memory = ', mem*700, ' # MB'
  write(fid2,'(A,I0,A)') 'mc.fcisolver.max_memory = ',mem*300,' # MB'
 else            ! multi-root DMRG-CASCI
  write(fid2,'(A)') ')'
  write(fid2,'(A,I0,A)') 'mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=',maxM,')'
  call prt_block_mem(fid2, mem, nproc, block_mpi)
 end if

 write(fid2,'(A)',advance='no') 'mc.fcisolver.nroots = '
 if(hardwfn) then
  write(fid2,'(I0)') nstate+4
 else if(crazywfn) then
  write(fid2,'(I0)') nstate+7
 else
  write(fid2,'(I0)') nstate+1
 end if

 if(.not. dmrgscf) then
  call prt_hard_or_crazy_casci_pyscf(fid2,nacta-nactb,hardwfn,crazywfn,.false.)
  if(.not. mixed_spin) write(fid2,'(A,F7.3,A)') 'mc.fix_spin_(ss=', ss, ')'
 end if
 write(fid2,'(A)') 'mc.verbose = 4'
 write(fid2,'(A)') 'mc.kernel(mo)'

 ! modified from pyscf-xxx/examples/mcscf/15-transition_dm.py
 if(dmrgscf) then
  write(fid2,'(/,A)') '# NTOs of DMRG-CASCI will be skipped since they are not &
                      &implemented yet.'
 else
  write(fid2,'(/,A)') '# calculate oscillator strengths and NTOs of transition &
                      &|0> -> |i>'
  write(fid2,'(A)') 'charges = mol.atom_charges()'
  write(fid2,'(A)') 'coords = mol.atom_coords()'
  write(fid2,'(A)') "charge_center = np.einsum('z,zx->x', charges, coords)/char&
                    &ges.sum()"
  write(fid2,'(A)') 'with mol.with_common_origin(charge_center):'
  write(fid2,'(A)') "  dip_int = mol.intor('int1e_r')"
  write(fid2,'(A)') 'nroots = mc.fcisolver.nroots'
  write(fid2,'(A)') 'nacto = mc.ncas'
  write(fid2,'(A)') 'idx1 = mc.ncore + 1'
  write(fid2,'(A)') 'idx2 = mc.ncore + nacto'
  write(fid2,'(A)') 'mo_cas = mo[:,mc.ncore:idx2].copy()'
  write(fid2,'(A)') 'for i in range(1,nroots):'
  i = INDEX(cmofch, '.fch', back=.true.)
  write(fid2,'(A)') "  part_fch = '"//cmofch(1:i-1)//"_NTO_P0'+str(i)+'.fch'"
  write(fid2,'(A)') "  hole_fch = '"//cmofch(1:i-1)//"_NTO_H0'+str(i)+'.fch'"
  write(fid2,'(A)') "  copyfile('"//TRIM(cmofch)//"', part_fch)"
  write(fid2,'(A)') "  copyfile('"//TRIM(cmofch)//"', hole_fch)"
  write(fid2,'(A)') '  tdm = mc.fcisolver.trans_rdm1(mc.ci[0],mc.ci[i], nacto, &
                    &mc.nelecas)'
  write(fid2,'(A)') '  ev, part_mo, hole_mo, fosc = gen_nto_and_fosc_from_mo_td&
                    &m(nbf, nacto, mo_cas, \'
  write(fid2,'(A)') '                               tdm, dip_int, mc.e_tot[i]-m&
                    &c.e_tot[0])'
  write(fid2,'(A)') '  noon[mc.ncore:idx2] = ev.copy()'
  write(fid2,'(A)') '  mo[:,mc.ncore:idx2] = part_mo.copy()'
  write(fid2,'(A)') "  py2fch(part_fch,nbf,nif,mo,'a',noon,False,False)"
  write(fid2,'(A)') '  mo[:,mc.ncore:idx2] = hole_mo.copy()'
  write(fid2,'(A)') "  py2fch(hole_fch,nbf,nif,mo,'a',noon,False,False)"
  write(fid2,'(A)') "  print('|0> -> |%d>'%(i),', fosc =',fosc)"
 end if

 if(nevpt2) then
  if(dmrgscf) then
   write(fid2,'(/,A)') '# State-specific DMRG-NEVPT2 based on multi-root DMRG-C&
                       &ASCI'
   call prt_dmrg_nevpt2_setting(fid2)
  else
   write(fid2,'(/,A)') '# State-specific NEVPT2 based on multi-root CASCI'
   write(fid2,'(A,I0)') 'nstate = ', nstate
   write(fid2,'(A)') 'for i in range(nstate+1):'
   write(fid2,'(A)') '  mrpt.NEVPT(mc, root=i).kernel()'
  end if
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

 i = INDEX(gjfname, '.gjf', back=.true.)
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
 use mr_keyword, only: mem, nproc, nevpt2, QD, FIC, DLPNO, F12, RI, RIJK_bas, &
  mixed_spin, nstate, hardwfn, crazywfn
 implicit none
 integer :: i, fid, fid1
 character(len=240) :: buf, mklname, gbwname, inpname1
 character(len=240), intent(in) :: inpname, hf_fch

 call fch2mkl_wrap(hf_fch)
 i = INDEX(hf_fch, '.fch', back=.true.)
 inpname1 = hf_fch(1:i-1)//'_o.inp'
 mklname = hf_fch(1:i-1)//'_o.mkl'

 i = INDEX(inpname, '.inp', back=.true.)
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
  if(QD .or. F12) then
   write(fid1,'(A)') ' PTSettings'
   if(QD) write(fid1,'(A)') '  QDType QD_VanVleck'
   if(F12) write(fid1,'(A)') '  F12 true'
   write(fid1,'(A)') ' end'
  end if
 end if

 write(fid1,'(A)') 'end'
 if(nevpt2) then
  write(fid1,'(A)') '%method'
  write(fid1,'(A)') ' FrozenCore FC_NONE'
  write(fid1,'(A)') 'end'
 end if

 !write(fid1,'(A)') '%mrci'
 !write(fid1,'(A)') ' tsel 0.0'
 !write(fid1,'(A)') ' tpre 0.0'
 !write(fid1,'(A)') ' Etol 1e-7'
 !write(fid1,'(A)') ' Rtol 1e-7'
 !write(fid1,'(A)') ' MaxIter 100'
 !write(fid1,'(A)') ' densities 0,1'
 !write(fid1,'(3(A,I0),A)') ' NewBlock 1 * nroots ',nstate+1,' excitations none &
 !                          &refs cas(',nacte,',',nacto,') end end'
 !write(fid1,'(A)') 'end'

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
 i = INDEX(hf_fch, '.fch', back=.true.)
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

! print SA-CASSCF input file of OpenMolcas
subroutine prt_sacas_molcas_inp(inpname, hf_fch)
 use mr_keyword, only: molcas_path
 use util_wrapper, only: fch2inporb_wrap
 implicit none
 character(len=240), intent(in) :: inpname, hf_fch

 call check_exe_exist(molcas_path)
 call fch2inporb_wrap(hf_fch, .false., inpname)
 call prt_cas_molcas_inp(inpname, .true.)
end subroutine prt_sacas_molcas_inp

! read SA-CASSCF energies from various output files
subroutine read_sa_cas_energies_from_output(cas_prog, outname, nstate, &
                                            sa_cas_e, ci_ssquare)
 implicit none
 integer, intent(in) :: nstate ! ground state included
 character(len=10), intent(in) :: cas_prog
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: sa_cas_e(0:nstate), ci_ssquare(0:nstate)

 select case(TRIM(cas_prog))
 case('gaussian')
!  call read_sa_cas_energies_from_gau_log(outname, nstate, sa_cas_e)
  ci_ssquare = 0d0
 case('pyscf')
  call read_sa_cas_e_from_pyscf_out(outname, nstate, sa_cas_e, ci_ssquare)
 case('orca')
  call read_sa_cas_e_from_orca_out(outname, nstate, sa_cas_e, ci_ssquare)
 case('openmolcas')
  call read_sa_cas_e_from_molcas_out(outname, nstate, sa_cas_e, ci_ssquare)
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_sa_cas_energies_from_output: &
                   &cas_prog cannot be recognized.'
  write(6,'(A)') 'cas_prog='//TRIM(cas_prog)
  stop
 end select
end subroutine read_sa_cas_energies_from_output

! read SA-CASSCF energies from PySCF output
subroutine read_sa_cas_e_from_pyscf_out(outname, nstate, sa_cas_e, ci_ssquare)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: sa_cas_e(0:nstate), ci_ssquare(0:nstate)

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
  write(6,'(/,A)') "ERROR in subroutine read_sa_cas_e_from_pyscf_out: no 'CASCI&
                   & energy for'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  stop
 end if

 do i = 0, nstate, 1
  read(fid,'(A)') buf
  j = INDEX(buf, '=')
  read(buf(j+1:),*) sa_cas_e(i)
  buf(j:j) = ' '
  j = INDEX(buf, '=')
  read(buf(j+1:),*) ci_ssquare(i)
 end do ! for i

 close(fid)
end subroutine read_sa_cas_e_from_pyscf_out

! read SA-CASSCF energies from ORCA output
subroutine read_sa_cas_e_from_orca_out(outname, nstate, sa_cas_e, ci_ssquare)
 implicit none
 integer :: i, j, k, m, mult, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: s
 real(kind=8), intent(out) :: sa_cas_e(0:nstate), ci_ssquare(0:nstate)

 sa_cas_e = 0d0; ci_ssquare = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:24) == 'CAS-SCF STATES FOR BLOCK') exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  write(6,'(/,A)') "ERROR in subroutine read_sa_cas_e_from_orca_out: no 'CAS-SC&
                  &F STATES FOR BLOCK'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  stop
 end if

 ! CAS-SCF STATES FOR BLOCK  1 MULT= 3 NROOTS= 3
 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) k

 j = INDEX(buf, '=')
 read(buf(j+1:),*) mult
 s = 0.5d0*DBLE(mult-1)
 ci_ssquare(0:k-1) = s*(s+1d0)

 m = -1
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == 'ROOT') then
   m = m + 1
   i = INDEX(buf, '=')
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

 i = INDEX(buf, '=', back=.true.)
 read(buf(i+1:),*) k

 j = INDEX(buf, '=')
 read(buf(j+1:),*) mult
 s = 0.5d0*DBLE(mult-1)
 ci_ssquare(m+1:) = s*(s+1d0)

 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:4) == 'ROOT') then
   m = m + 1
   i = INDEX(buf, '=')
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
    s = ci_ssquare(i); ci_ssquare(i) = ci_ssquare(j); ci_ssquare(j) = s
   end if
  end do ! for i
 end do ! for j
end subroutine read_sa_cas_e_from_orca_out

! read SA-CASSCF energies from OpenMolcas output
subroutine read_sa_cas_e_from_molcas_out(outname, nstate, sa_cas_e, ci_ssquare)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: s
 real(kind=8), intent(out) :: sa_cas_e(0:nstate), ci_ssquare(0:nstate)

 sa_cas_e = 0d0; ci_ssquare = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(7:19) == 'Final state e') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_sa_cas_e_from_molcas_out: no 'Fina&
                   &l state e' found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf ! skip 2 lines
 read(fid,'(A)') buf

 do i = 0, nstate, 1
  read(fid,'(A)') buf
  j = INDEX(buf, ':', back=.true.)
  read(buf(j+1:),*) sa_cas_e(i)
 end do ! for i

 i = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(4:12) == 'SPIN MULT') then
   read(buf(40:),*) j
   s = 0.5d0*DBLE(j-1)
   ci_ssquare(i) = s*(s+1d0)
   i = i + 1
   if(i == nstate+1) exit
  end if
 end do ! for while

 close(fid)
end subroutine read_sa_cas_e_from_molcas_out

! read multi-root CASCI-based NEVPT2 energies from a PySCF out file
subroutine read_multiroot_nevpt2_from_pyscf(outname, dmrg, nstate, nevpt2_e)
 use phys_cons, only: au2ev
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), allocatable :: casci_e(:)
 real(kind=8), intent(out) :: nevpt2_e(0:nstate)
 logical :: no_conv
 logical, intent(in) :: dmrg

 no_conv = .false.
 allocate(casci_e(0:nstate), source=0d0)

 if(dmrg) then ! DMRG-NEVPT2 based on multi-root DMRG-CASCI
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:12) == 'CASCI state ') then
    BACKSPACE(fid)
    exit
   end if
  end do ! for while

  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_multiroot_nevpt2_from_pyscf: no '&
                    &CASCI state ' found"
   write(6,'(A)') 'in file '//TRIM(outname)
   close(fid)
   stop
  end if
 else          ! NEVPT2 based on multi-root CASCI
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
   write(6,'(/,A)') "ERROR in subroutine read_multiroot_nevpt2_from_pyscf: no '&
                    &CASCI conv' found in"
   write(6,'(A)') 'file '//TRIM(outname)
   close(fid)
   stop
  end if

  if(no_conv) then
   write(6,'(A)') REPEAT('-',79)
   write(6,'(A)') 'Warning in subroutine read_multiroot_nevpt2_from_pyscf: CASC&
                  &I not converged.'
   write(6,'(A)') 'You should be cautious about the CASCI and NEVPT2 energies. &
                  &You are recommended'
   write(6,'(A)') 'to open file '//TRIM(outname)//' and check.'
   write(6,'(A)') 'The program will not stop but continue...'
   write(6,'(A)') REPEAT('-',79)
  end if
 end if

 do i = 0, nstate, 1
  read(fid,'(A)') buf
  j = INDEX(buf,'=')
  read(buf(j+1:),*) casci_e(i)
 end do ! for i

 i = 0
 do while(.true.)
  if(i == nstate+1) exit
  read(fid,'(A)') buf
  if(buf(1:8) == 'Nevpt2 E') then
   j = INDEX(buf,'=')
   read(buf(j+1:),*) nevpt2_e(i)
   i = i + 1
  end if
 end do ! for while

 close(fid)
 nevpt2_e = casci_e + nevpt2_e
 deallocate(casci_e)
end subroutine read_multiroot_nevpt2_from_pyscf

! read multi-root CASCI-based NEVPT2 energies from an ORCA out file
subroutine read_multiroot_nevpt2_from_orca(outname, nstate, nevpt2_e)
 use phys_cons, only: au2ev
 implicit none
 integer :: i, j, k, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: nevpt2_e(0:nstate)

 nevpt2_e = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 i = 0
 do while(.true.)
  read(fid,'(A)',iostat=j) buf
  if(j /= 0) exit
  if(index(buf,'Total Energy (') > 0) then
   k = INDEX(buf, '=')
   read(buf(k+1:),*) nevpt2_e(i)
   i = i + 1
   if(i == nstate+1) exit
  end if
 end do ! for while

 close(fid)
end subroutine read_multiroot_nevpt2_from_orca

! read oscillator strengths from various output file
subroutine read_fosc_from_output(cas_prog, outname, nstate, fosc)
 implicit none
 integer, intent(in) :: nstate
 character(len=10), intent(in) :: cas_prog
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: fosc(nstate)

 select case(TRIM(cas_prog))
 case('pyscf')
  call read_fosc_from_pyscf_out(outname, nstate, fosc)
 case('orca')
  call read_fosc_from_orca_out(outname, nstate, fosc)
 case('openmolcas')
  call read_fosc_from_molcas_out(outname, nstate, fosc)
 case default
  write(6,'(/,A)') 'ERROR in subroutine read_fosc_from_output: program cannot b&
                   &e recognized.'
  write(6,'(A)') 'CASCI_prog or CASSCF_prog='//TRIM(cas_prog)
  stop
 end select
end subroutine read_fosc_from_output

! read oscillator strengths from PySCF output file
subroutine read_fosc_from_pyscf_out(outname, nstate, fosc)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: fosc(nstate)

 fosc = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == '|0>') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_fosc_from_pyscf_out: no '|0>' foun&
                   &d in file "//TRIM(outname)
  close(fid)
  stop
 end if

 BACKSPACE(fid)
 do i = 1, nstate, 1
  read(fid,'(A)') buf
  j = INDEX(buf, '=')
  read(buf(j+1:),*) fosc(i)
 end do ! for i

 close(fid)
end subroutine read_fosc_from_pyscf_out

! read oscillator strengths from ORCA output file
subroutine read_fosc_from_orca_out(outname, nstate, fosc)
 implicit none
 integer :: i, k, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8), intent(out) :: fosc(nstate)
 logical :: orca6

 orca6 = .false.; fosc = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(36:48) == 'O   R   C   A') then
   i = -1
   exit
  end if
  if(buf(1:8) == '  States') exit
  if(buf(6:55) == 'Transition      Energy     Energy  Wavelength fosc') then
   orca6 = .true.
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_fosc_from_orca_out: failed to read&
                   & fosc in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf ! skip 2 lines
 read(fid,'(A)') buf

 k = 39
 if(orca6) k = 49
 do i = 1, nstate, 1
  read(fid,'(A)') buf
  read(buf(k:),*) fosc(i)
 end do ! for i

 close(fid)
end subroutine read_fosc_from_orca_out

! read oscillator strengths from OpenMolcas output file
subroutine read_fosc_from_molcas_out(outname, nstate, fosc)
 implicit none
 integer :: i, j, fid
 integer, intent(in) :: nstate
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 real(kind=8) :: rtmp
 real(kind=8), intent(out) :: fosc(nstate)

 fosc = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:11) == '++ Dipole t') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_fosc_from_molcas_out: no '++ Dipol&
                   &e t' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, 4   ! skip 4 lines
  read(fid,'(A)') buf
 end do ! for i

 ! it is possible that all oscillator strengths < 1e-5, simply return
 if(index(buf,'is only') > 0) then
  close(fid)
  return
 else
  read(fid,'(A)') buf ! skip 1 line
 end if

 !do i = 1, nstate, 1
 ! read(fid,*) j, j, fosc(i)
 !end do ! for i
 do while(.true.)
  read(fid,'(A)') buf
  if(index(buf,'---') > 0) exit
  read(buf,*) i, j, rtmp
  if(i == 1) fosc(j-1) = rtmp
 end do ! for while

 close(fid)
end subroutine read_fosc_from_molcas_out

! copy NTOs from OpenMolcas orbital files to .fch files
subroutine copy_nto_from_orb2fch(hf_fch, nstate)
 use util_wrapper, only: orb2fch_wrap
 implicit none
 integer :: i, k
 integer, intent(in) :: nstate
 character(len=240) :: part_orb, hole_orb, part_fch, hole_fch
 character(len=240), intent(in) :: hf_fch

 call find_specified_suffix(hf_fch, '.fch', k)

 do i = 2, nstate+1, 1
  write(part_orb,'(A,I0,A)') hf_fch(1:k-1)//'_SA-CAS.NTOrb.1_',i,'.a.PART'
  write(hole_orb,'(A,I0,A)') hf_fch(1:k-1)//'_SA-CAS.NTOrb.1_',i,'.a.HOLE'
  write(part_fch,'(A,I0,A)') hf_fch(1:k-1)//'_SA-CAS_NTO_P0',i-1,'.fch'
  write(hole_fch,'(A,I0,A)') hf_fch(1:k-1)//'_SA-CAS_NTO_H0',i-1,'.fch'
  call copy_file(hf_fch, part_fch, .false.)
  call copy_file(hf_fch, hole_fch, .false.)

  call check_molcas_nto_file_exist(part_orb)
  call orb2fch_wrap(part_orb, part_fch, .true.)

  call check_molcas_nto_file_exist(hole_orb)
  call orb2fch_wrap(hole_orb, hole_fch, .true.)

  call delete_files(2, [part_orb, hole_orb])
 end do ! for i
end subroutine copy_nto_from_orb2fch

! Check whether a specied NTO orbital file exists. If not, change the filename
!  and check again. If still not found, copy it from the $MOLCAS_WORKDIR/
!  directory.
! If you are running on a Cluster and the $MOLCAS_WORKDIR/ is network-mounted,
!  this could cause latency or timestamp mismatches. In such case OpenMolcas
!  might fail to copy some of the NTO files back into current directory. So we
!  need to do the copy (see https://gitlab.com/Molcas/OpenMolcas/-/issues/425)
subroutine check_molcas_nto_file_exist(orbname)
 implicit none
 integer :: i, k
 character(len=600) :: scr_path
 character(len=240) :: orbname1
 character(len=240), intent(inout) :: orbname ! allowed to be changed
 logical :: alive

 scr_path = ' '
 k = LEN_TRIM(orbname)
 inquire(file=orbname(1:k),exist=alive)
 if(alive) return

! For OpenMolcas < 23.10, NTO filenames are like xx.NTOrb.1_2.a.HOLE
! For OpenMolcas >=23.10, NTO filenames are like xx.NTOrb.SF.1_2.a.HOLE
! 'SF' for Spin-Free. Here we have to check two types of filenames.
 call find_specified_suffix(orbname, '.NTOrb', i)
 orbname1 = orbname(1:i+5)//'.SF'//orbname(i+6:k)
 inquire(file=orbname1(1:k+3),exist=alive)
 if(alive) then
  orbname = orbname1
  return
 end if

 call getenv('MOLCAS_WORKDIR', scr_path)
 scr_path = TRIM(scr_path)//'/'//orbname(1:i-1)

 inquire(file=TRIM(scr_path)//'/'//orbname(1:k),exist=alive)
 if(alive) then
  call copy_file_from_a_path(scr_path, orbname)
 else
  inquire(file=TRIM(scr_path)//'/'//orbname1(1:k+3),exist=alive)
  if(alive) then
   call copy_file_from_a_path(scr_path, orbname1)
   orbname = orbname1
  else
   write(6,'(/,A)') 'ERROR in subroutne check_molcas_nto_file_exist: NTO files &
                    &not found.'
   write(6,'(A)') 'Filename='//orbname(1:k)
   write(6,'(A)') 'Filename='//orbname1(1:k+3)
  end if
 end if
end subroutine check_molcas_nto_file_exist

