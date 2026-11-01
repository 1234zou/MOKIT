! written by jxzou at 20261028: moved here from file automr.f90

! perform the RASCI/RASSCF calculation
subroutine do_ras(scf)
 use mol, only: charge, mult, nbf, nif, ndb, nopen, nacta, nactb, nacte, nacto,&
  npair, ras1_nacto, ras3_nacto, ras_nhole, ras_npart, ptchg_e, nuc_pt_e, &
  rasci_e, rasscf_e
 use mr_keyword, only: mem, nproc, ist, iroot, nstate, new_mult, RI, DKH2, X2C,&
  rasci, rasscf, nacte_wish, nacto_wish, ras1_nacto_wish, ras3_nacto_wish, &
  ras_nhole_wish, ras_npart_wish, rasci_force, rasscf_force, molcas_omp, &
  dalton_mpi, bgchg, chgname, rasci_prog, rasscf_prog, molcas_path, orca_path, &
  molpro_path, gau_path, hf_fch, datname, casnofch
 use util_wrapper, only: add_bgcharge2inp_wrap, fch2inporb_wrap, orb2fch_wrap,&
  fch2dal_wrap, fch2mkl_wrap, mkl2gbw, gbw2molden, molden2fch_wrap, fch2com_wrap,&
  xml2fch_wrap, formchk, unfchk
 implicit none
 integer :: i, na, nb, nvir, n_pocc, nacte1, nacto1, int_tmp(6), RENAME
 real(kind=8) :: unpaired_e ! unpaired electrons
 real(kind=8) :: e(2) ! e(1) is RASCI energy, e(2) is RASSCF energy
 real(kind=8), allocatable :: noon(:), ras_on(:)
 character(len=10) :: ras_prog
 character(len=24) :: data_string = ' '
 character(len=28), parameter :: error_warn = 'ERROR in subroutine do_ras: '
 character(len=240) :: chkname, fchname, proname, inpname, outname, orbname, &
  mklname, molden, xmlname, gradname
 character(len=300) :: buf
 logical :: uhf, ras_force
 logical, intent(in) :: scf

 if(scf) then
  if(.not. rasscf) return
  if(TRIM(rasscf_prog) == 'orca') then
   write(6,'(/,A)') error_warn//'RASSCF not implemented in ORCA so far.'
   write(6,'(A)') 'If you insist on using ORCA, you can perform the RASCI calcu&
                  &lation instead.'
   write(6,'(A)') 'If you are willing to use another RASSCF_prog, you can speci&
                  &fy any of'
   write(6,'(A)') 'OpenMolcas, Dalton, Molpro, Gaussian.'
   stop
  end if
 else
  if(.not. rasci) return
 end if
 write(6,'(//,A)') 'Enter subroutine do_ras...'

 ras_prog = rasci_prog
 if(scf) ras_prog = rasscf_prog
 if(X2C .and. TRIM(ras_prog)=='gaussian') then
  write(6,'(/,A)') error_warn//'X2C is not supported when RASCI_prog/RASSCF_pro&
                  &g=Gaussian.'
  write(6,'(A)') 'You can specify OpenMolcas instead.'
  stop
 end if
 ras_force = (rasci_force .or. rasscf_force)

 if(ist == 5) then
  call check_uhf_in_fch(hf_fch, uhf)
  if(uhf) then
   write(6,'(/,A)') error_warn//'R(O)HF-type .fch(k) file is expected when'
   write(6,'(A)') 'using the `readno` keyword. But you provide a UHF-type file.&
                  & Please check'
   write(6,'(A)') TRIM(hf_fch)
   stop
  end if

  ! if the user wants to perform RASCI/RASSCF calculations using a new spin
  ! multiplicity, we copy a new .fch file and update mult as well as hf_fch
  if(new_mult>0 .and. new_mult/=mult) then
   if(MOD(IABS(new_mult-mult),2) /= 0) then
    write(6,'(/,A)') error_warn//'inconsistency found between mult and NewMult.'
    write(6,'(2(A,I0))') 'mult=', mult, ', NewMult=', new_mult
    stop
   end if
   write(6,'(A)') REPEAT('-', 79)
   write(6,'(A)') 'Remark: a new spin multiplicity is specified by the user. Up&
                  &dating...'
   call find_specified_suffix(hf_fch, '.fch', i)
   fchname = hf_fch(1:i-1)//'_M.fch'
   call sys_copy_file(TRIM(hf_fch), TRIM(fchname), .false.)
   call modify_charge_and_mult_in_fch(fchname, charge, new_mult)
   mult = new_mult
   hf_fch = fchname
   write(6,'(A)') 'Done. New filenames will be augmented with `_M`.'
   write(6,'(A)') REPEAT('-', 79)
  end if

  write(6,'(A)') 'Radical index for input NOs:'
  call calc_unpaired_from_fch(hf_fch, 3, .false., unpaired_e)
 else ! ist /= 5
  if(new_mult > 0) then
   write(6,'(/,A)') error_warn//'NewMult is supposed to be used along with ist=5.'
   write(6,'(A,I0)') 'But current ist=', ist
   stop
  end if
 end if

 call read_nbf_and_nif_from_fch(hf_fch, nbf, nif)
 call read_na_and_nb_from_fch(hf_fch, na, nb)
 nopen = na - nb
 call determine_cas_ras_proname(ist, npair, .true., scf, datname, hf_fch, &
                                fchname, proname)

 int_tmp = [nacte_wish, ras_nhole_wish, ras_npart_wish, ras1_nacto_wish, &
            nacto_wish, ras3_nacto_wish]
 if(ALL(int_tmp == 0)) then
  ras_nhole = 2; ras_npart = 2
  allocate(noon(nif))
  call read_eigenvalues_from_fch(fchname, nif, 'a', noon)
  call determine_ras_active_space_from_noon(nif, noon, ndb, ras1_nacto, nacto, &
                                            ras3_nacto)
  nacte = SUM(noon) - 2d0*(ndb+ras1_nacto)
 else
  nacte = nacte_wish ! active electrons in RAS2
  nacto = nacto_wish ! active orbitals in RAS2
  ras_nhole = ras_nhole_wish
  ras_npart = ras_npart_wish
  ras1_nacto = ras1_nacto_wish
  ras3_nacto = ras3_nacto_wish
  ndb = (na + nb - nacte - 2*ras1_nacto)/2
 end if

 nacte1 = 2*ras1_nacto + nacte
 nacto1 = ras1_nacto + nacto + ras3_nacto
 n_pocc = ndb + nacto1
 nvir = nif - n_pocc
 nacta = (nacte + nopen)/2
 nactb = nacte - nacta

 if(scf) then
  write(6,'(A)',advance='no') 'RASSCF('
 else
  write(6,'(A)',advance='no') 'RASCI('
 end if
 write(6,'(6(I0,A))') nacte1,',',ras_nhole,',',ras_npart,';',ras1_nacto,',', &
                       nacto,',',ras3_nacto,') using program '//TRIM(ras_prog)

 write(6,'(6(A,I0),2(A,L1))') 'ndb=',ndb,', RAS1=',ras1_nacto,', RAS2=',nacto,&
  ', RAS3=',ras3_nacto,', nvir=',nvir,', Root=',iroot,', RI=',RI,', Force=', &
  ras_force
 casnofch = TRIM(proname)//'_NO.fch'

 select case(TRIM(ras_prog))
 case('openmolcas')
  call check_exe_exist(molcas_path)
  inpname = TRIM(proname)//'.input'
  outname = TRIM(proname)//'.out'
  call fch2inporb_wrap(fchname, .false., inpname)
  write(orbname,'(A,I0)') TRIM(proname)//'.RasOrb.', iroot+1
  call prt_molcas_ras_inp(inpname, nacte, nacto, ras1_nacto, ras3_nacto, &
                          ras_nhole, ras_npart, scf, ras_force, .false.)
  if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
  call submit_molcas_job(inpname, mem, nproc, molcas_omp)
  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  call orb2fch_wrap(orbname, casnofch, .true., .false., 0)
  ! OpenMolcas CASSCF NOs may not be in ascending order, sort them
  call sort_no_by_noon(casnofch, ndb+1, n_pocc)

 case('dalton')
  inpname = TRIM(proname)//'.dal'
  outname = TRIM(proname)//'.out'
  call fch2dal_wrap(fchname, inpname, .true.)
  call prt_dalton_ras_inp(inpname, ndb, nacte, nacto, ras1_nacto, ras3_nacto, &
                          ras_nhole, ras_npart, scf, ras_force)
  if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
  call submit_dalton_job(proname, mem, nproc, dalton_mpi, .false., .false., &
                         .false.)

  ! untar/unzip the compressed package to get DALTON.MOPUN
  buf = 'tar -xpf '//TRIM(proname)//'.tar.gz DALTON.MOPUN'
  call run_command(TRIM(buf), .false., .true.)

  ! make a copy to save NOs
  call copy_file(fchname, casnofch, .false.)

  ! transfer NOs (and NOONs, if any) from DALTON.MOPUN to .fch file
  buf = 'dal2fch DALTON.MOPUN '//TRIM(casnofch)
  if(.not. ras_force) buf = TRIM(buf)//' -no'
  ! when force=.T., NOONs are not in DALTON.MOPUN, so we only transfer NOs here
  call run_command(TRIM(buf), .false., .true.)

  ! Now transfer NOONs when force=.T.
  if(ras_force) then
   allocate(noon(nif))
   call read_on_from_dalton_out(outname, n_pocc, noon(1:n_pocc))
   call write_eigenvalues_to_fch(casnofch, nif, 'a', noon, .true.)
   deallocate(noon)
  end if
  call delete_file('DALTON.MOPUN')

 case('orca')
  call check_exe_exist(orca_path)
  inpname = TRIM(proname)//'.inp'
  mklname = TRIM(proname)//'.mkl'
  outname = TRIM(proname)//'.out'
  gradname = TRIM(proname)//'.engrad'
  orbname = TRIM(proname)//'.gbw'
  molden = TRIM(proname)//'.molden'
  call fch2mkl_wrap(fchname, mklname, REPEAT(' ',30), .true.)
  call prt_orca_ras_inp(inpname, nacte, nacto, ras1_nacto, ras3_nacto, &
                        ras_nhole, ras_npart, .false.)
  if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
  ! if bgchg = .True., .inp and .mkl file will be updated
  call mkl2gbw(mklname)
  call delete_file(TRIM(mklname))
  if(ras_force) call add_force_key2orca_inp(inpname)

  call submit_orca_job(orca_path, inpname, .true., .false., .false.)
  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs

  call gbw2molden(orbname)
  call delete_file(TRIM(orbname))
  orbname = TRIM(proname)//'.fch'
  call molden2fch_wrap(molden, orbname, 'orca', .true.)
  allocate(ras_on(nacto1))
  call read_rasci_noon_from_orca_out(outname, nacto1, ras_on)
  allocate(noon(nif), source=0d0)
  if(ndb > 0) noon(1:ndb) = 2d0
  noon(ndb+1:ndb+nacto1) = ras_on
  deallocate(ras_on)
  call write_eigenvalues_to_fch(orbname, nif, 'a', noon, .true.)
  deallocate(noon)
  call copy_ev_and_mo_between_fch(orbname, casnofch, '-aa')
  call delete_files(3, [mklname, orbname, molden])

 case('molpro')
  call check_exe_exist(molpro_path)
  molden = fchname
  inpname = TRIM(proname)//'.com'
  outname = TRIM(proname)//'.out'
  xmlname = TRIM(proname)//'.xml'
  orbname = TRIM(proname)//'.fch'
  call fch2com_wrap(fchname, inpname)
  call convert2molpro_fname(molden, '.a')
  call convert2molpro_fname(orbname, '.a')
  i = RENAME(TRIM(molden), TRIM(orbname))
  call prt_molpro_ras_inp(inpname, ndb, nacte, nacto, ras1_nacto, ras3_nacto, &
                          ras_nhole, ras_npart, scf, ras_force)
  if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
  call submit_molpro_job(inpname, mem, nproc)
  call copy_file(fchname, casnofch, .false.) ! make a copy to save NOs
  call xml2fch_wrap(xmlname, casnofch, .true.)

 case('gaussian')
  call check_exe_exist(gau_path)
  inpname = TRIM(proname)//'.gjf'
  outname = TRIM(proname)//'.log'
  chkname = TRIM(proname)//'.chk'
  call prt_ras_gjf(inpname, nacte, nacto, ras1_nacto, ras3_nacto, ras_nhole, &
                   ras_npart, scf, ras_force, DKH2, .false.)
  if(bgchg) call add_bgcharge2inp_wrap(chgname, inpname)
  call unfchk(fchname, chkname)
  call submit_gau_job(gau_path, inpname, .true.)
  call formchk(chkname, casnofch)
  call delete_file(TRIM(chkname))
  call modify_irohf_in_fch(casnofch, 0)

 case default
  write(6,'(/,A)') error_warn//'invalid RASCI_prog/RASSCF_prog='//TRIM(ras_prog)
  write(6,'(A)') 'Currently only OpenMolcas/Gaussian are supported.'
  stop
 end select

 ! generate RASCI/RASSCF density from NOs and NOONs
 call update_density_using_no_and_on(casnofch)

 ! read energy, check convergence and check spin
 select case(TRIM(ras_prog))
 case('orca')
  ! TODO: consider background point charges
  ! TODO: consider multiple roots
  allocate(noon(0:nstate))
  call read_rasci_e_from_orca_out(outname, nstate, noon)
  e(1) = noon(0) ! RASCI energy
  deallocate(noon)
 case default
  call read_cas_energy_from_output(ras_prog, outname, mult, scf, .false., &
                                   ptchg_e, nuc_pt_e, e)
 end select

 rasci_e = e(1)
 write(6,'(/,A,F18.8,1X,A4)') 'E(RASCI)  = ', e(1), 'a.u.'
 if(scf) then
  rasscf_e = e(2)
  write(6,'(A,F18.8,1X,A4)') 'E(RASSCF) = ', e(2), 'a.u.'
 end if

 call calc_unpaired_from_fch(casnofch, 3, .false., unpaired_e)
 if(scf .and. TRIM(ras_prog)=='dalton') then
  write(6,'(/,A)') REPEAT('-',79)
  write(6,'(A)') 'Remark: in Dalton RASSCF, all off-diagonal elements of the on&
                 &e-electron density'
  write(6,'(A)') 'matrix which couple different RAS spaces are ignored. So the &
                 &natural orbital'
  write(6,'(A)') 'occupation numbers in xxx_RASSCF_NO.fch will be slightly diff&
                 &erent from those'
  write(6,'(A)') 'calculated by OpenMolcas/Gaussian. Usually it does not cause &
                 &any problem.'
  write(6,'(A)') REPEAT('-',79)
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_ras at '//TRIM(data_string)
end subroutine do_ras

subroutine determine_ras_active_space_from_noon(nif, noon, ndb, ras1_nacto, &
                                                nacto, ras3_nacto)
 implicit none
 integer :: i
 integer, intent(in) :: nif
 integer, intent(out) :: ndb, ras1_nacto, nacto, ras3_nacto
 real(kind=8), intent(in) :: noon(nif)

 ndb = 0; ras1_nacto = 0; nacto = 0; ras3_nacto = 0

 do i = 1, nif, 1
  if(noon(i)>1.98d0) then
   ndb = ndb + 1
  else if(noon(i)>1.9d0 .and. noon(i)<=1.98d0) then
   ras1_nacto = ras1_nacto + 1
  else if(noon(i)>=0.1d0 .and. noon(i)<=1.9d0) then
   nacto = nacto + 1
  else if(noon(i)>=0.02d0 .and. noon(i)<0.1d0) then
   ras3_nacto = ras3_nacto + 1
  end if
 end do ! for i
end subroutine determine_ras_active_space_from_noon

! print RASCI/RASSCF and RI (if required) keywords in to a given (Open)Molcas
! input file
subroutine prt_molcas_ras_inp(inpname, nacte, nacto, ras1_nacto, ras3_nacto, &
                              ras_nhole, ras_npart, scf, force, polar)
 use mr_keyword, only: dmrgci, dmrgscf, RI, hardwfn, crazywfn, xmult, iroot, &
  nstate
 use mol, only: mult, charge
 implicit none
 integer :: i, nacte1, nroots, lroots, fid1, fid2, RENAME
 integer, intent(in) :: nacte, nacto, ras1_nacto,ras3_nacto,ras_nhole,ras_npart
 integer, allocatable :: weight(:)
 character(len=40), parameter :: error_warn = 'ERROR in subroutine prt_molcas_r&
                                              &as_inp: '
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf, force, polar

 nacte1 = nacte + 2*ras1_nacto
 if(dmrgci .or. dmrgscf) then
  write(6,'(/,A)') error_warn//'DMRG cannot be combined with RASCI/RASSCF'
  write(6,'(A)') 'currently. Please check your input settings.'
  stop
 end if

 ! Since MOKIT v1.2.7rc4, we use RICD rather than RIJK in OpenMolcas, if RI is
 ! required.
 ! TODO: check whether we need to use a larger RICD auxbasis
 if(RI) call add_RI_kywd_into_molcas_inp(inpname, .true.)

 call find_specified_suffix(inpname, '.in', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4)=="&SCF" .or. buf(1:7)=="&RASSCF") exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') error_warn//'`&SCF`/`&RASSCF` not found in file'
  write(6,'(A)') TRIM(inpname)
  close(fid1)
  close(fid2,status='delete')
  stop
 end if

 if(iroot>0 .and. mult/=xmult) then
  write(6,'(/,A)') 'ERROR in subroutine prt_molcas_cas_para: Xmult/=mult!'
  write(6,'(A)') 'When calling OpenMolcas, the ground state spin must be equal &
                 &to the excited'
  write(6,'(A)') 'state spin.'
  stop
 end if

 write(fid2,'(/,A)') "&RASSCF"
 write(fid2,'(A)') 'CIMX= 200'
 write(fid2,'(A)') 'Tight= 5d-8 5d-6'
 if(crazywfn) then
  write(fid2,'(A)') 'SDav= 500'
 else if(hardwfn) then
  write(fid2,'(A)') 'SDav= 300'
 end if

 write(fid2,'(A,I0)') 'Spin= ', mult
 write(fid2,'(A,I0)') 'Charge= ', charge
 write(fid2,'(A,3(1X,I0))') 'NACTEL=', nacte1, ras_nhole, ras_npart
 write(fid2,'(A,I0)') 'RAS1= ', ras1_nacto
 write(fid2,'(A,I0)') 'RAS2= ', nacto
 write(fid2,'(A,I0)') 'RAS3= ', ras3_nacto
 !TODO: we can also consider RAS4 in the future
 if(iroot > 9) then
  close(fid2)
  write(6,'(/,A)') 'ERROR in subroutine prt_molcas_cas_para: iroot>9. Please c&
                   &ontact the MOKIT'
  write(6,'(A)') 'developers to modify the input file print format.'
  stop
 end if

 ! SS-RASCI
 if(iroot > 0) then
  nroots = iroot + 1
  lroots = iroot + 3
  allocate(weight(nroots), source=0)
  weight(nroots) = 1
  write(fid2,'(A,2(1X,I0),A)',advance='no') 'CIroot=',nroots,lroots,';'
  write(fid2,'(9I2)',advance='no') (i,i=1,nroots)
  write(fid2,'(A,9I2)') ';', weight
  deallocate(weight)
 end if

 ! multi-root RASCI
 if(nstate > 0) then
  nroots = nstate + 1
  lroots = nstate + 3
  write(fid2,'(A,2(1X,I0),A)') 'CIroot=',nroots,lroots,' 1'
 end if

 call find_specified_suffix(inpname, '.in', i)
 write(fid2,'(A)') 'FILEORB= '//inpname(1:i-1)//'.INPORB'
 if(.not. scf) write(fid2,'(A)') 'CIonly'

 if(nstate > 0) then   ! generate hole and particle NTOs
  write(fid2,'(/)',advance='no')
  ! these CASCI states share the same basis set and the same MOs, so copy from
  ! one common $Project.JobIph
  do i = 1, nstate+1, 1
   write(fid2,'(A,I3.3)') '>> COPY $Project.JobIph JOB',i
  end do ! for i
  write(fid2,'(/,A)') "&RASSI"
  write(fid2,'(A)') 'NTOC'
  write(fid2,'(A,I0,99I2)') 'Nr of JobIphs=',nstate+1,(1,i=1,nstate+1)
  do i = 1, nstate, 1
   write(fid2,'(I0,A1)',advance='no') i,';'
  end do ! for i
  write(fid2,'(I0)') nstate+1
 end if

 if(force) write(fid2,'(/,A)') '&ALASKA'
 write(fid2,'(/)',advance='no')
 close(fid2)
 close(fid1,status='delete')
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_molcas_ras_inp

! print RASCI/RASSCF keywords into a given Dalton input file
subroutine prt_dalton_ras_inp(dalname, ndb, nacte, nacto, ras1_nacto, ras3_nacto,&
                              ras_nhole, ras_npart, scf, force)
 use mol, only: mult
 use mr_keyword, only: DKH2, hardwfn, crazywfn
 implicit none
 integer :: i, nacte1, fid, fid1, RENAME
 integer, intent(in) :: ndb,nacte,nacto,ras1_nacto,ras3_nacto,ras_nhole,ras_npart
 character(len=240) :: buf, dalname1
 character(len=240), intent(in) :: dalname
 logical, intent(in) :: scf, force

 nacte1 = nacte + 2*ras1_nacto
 call find_specified_suffix(dalname, '.dal', i)
 dalname1 = dalname(1:i-1)//'.t'
 open(newunit=fid1,file=TRIM(dalname1),status='replace')

 write(fid1,'(A)') '**DALTON INPUT'
 if(DKH2) write(fid1,'(A)') '.DOUGLAS-KROLL'
 if(force) then
  write(fid1,'(A)') '.RUN PROPERTIES'
 else
  write(fid1,'(A)') '.RUN WAVE FUNCTIONS'
 end if
 write(fid1,'(A,/,A)') '**WAVE FUNCTIONS','.HF'
 if(scf) then
  write(fid1,'(A)') '.MCSCF'
 else
  write(fid1,'(A)') '.CI'
 end if
 write(fid1,'(A,/,A)') '*SCF INPUT','.NODIIS'
 write(fid1,'(A,/,A)') '.NOQCSCF','.NONCANONICAL'
 write(fid1,'(A,/,A)') '*CONFIGURATION INPUT', '.SPIN MULTIPLICITY'
 write(fid1,'(I0)') mult
 write(fid1,'(A,/,I0)') '.INACTIVE ORBITALS', ndb
 write(fid1,'(A,/,I0)') '.RAS1 SPACE', ras1_nacto
 write(fid1,'(A,/,I0)') '.RAS2 SPACE', nacto
 write(fid1,'(A,/,I0)') '.RAS3 SPACE', ras3_nacto
 write(fid1,'(A,/,A1,1X,I0)') '.RAS1 HOLES', '0', ras_nhole
 write(fid1,'(A,/,A1,1X,I0)') '.RAS3 ELECTRONS', '0', ras_npart
 write(fid1,'(A,/,I0)') '.ELECTRONS', nacte1
 if(scf) then
  write(fid1,'(A)') '*OPTIMIZATION'
  write(fid1,'(A,/,A)') '.MAX CI','500'
  write(fid1,'(A,/,A)') '.MAX MACRO ITERATIONS','50'
  write(fid1,'(A,/,A)') '.MAX MICRO ITERATIONS','200'
  if(hardwfn) then
   write(fid1,'(A,/,A)') '.CI PHP MATRIX','500'
  else if(crazywfn) then
   ! Dalton MCSCF uses CSF by default, a large figure is not needed
   write(fid1,'(A,/,A)') '.CI PHP MATRIX','1000'
  else
   write(fid1,'(A,/,A)') '.CI PHP MATRIX','200'
  end if
 else
  write(fid1,'(A)') '*CI INPUT'
  write(fid1,'(A,/,A)') '.MAX ITERATIONS','500'
 end if
 write(fid1,'(A,/,A)') '*PRINT LEVELS', '.CANONI'

 open(newunit=fid,file=TRIM(dalname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:12) == '*ORBITAL INP') exit
 end do ! for while
 BACKSPACE(fid)

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12) == '**END OF INP') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine prt_cas_dalton_inp: no "**END OF INP" f&
                   &ound in file'
  write(6,'(A)') TRIM(dalname)
  close(fid1,status='delete')
  close(fid)
  stop
 end if

 if(force) write(fid1,'(A,/,A)') '**PROPERTIES', '.MOLGRA'
 write(fid1,'(A)') '**END OF INPUT'
 close(fid,status='delete')
 close(fid1)

 i = RENAME(TRIM(dalname1), TRIM(dalname))
end subroutine prt_dalton_ras_inp

! Print RASCI keywords in to a given ORCA .inp file. We will import
! `RIJK_bas` and `RIC_bas` from the module mr_keyword, but not import `RI`
! from that module. In such case, we have better control of `RI`.
subroutine prt_orca_ras_inp(inpname, nacte, nacto, ras1_nacto, ras3_nacto, &
                            ras_nhole, ras_npart, polar)
 use mol, only: mult
 use mr_keyword, only: mem, nproc, RI, RIJK_bas, hardwfn, crazywfn, iroot, &
  nstate, xmult
 implicit none
 integer :: i, k, mem1, nproc1, nacte1, target_mult, fid1, fid2, RENAME
 integer, intent(in) :: nacte, nacto, ras1_nacto,ras3_nacto,ras_nhole,ras_npart
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: polar

 if(iroot > 19) then
  write(6,'(/,A)') 'ERROR in subroutine prt_orca_ras_inp: please contact the MO&
                   &KIT developers'
  write(6,'(A)') 'to modify the printing format.'
  write(6,'(A,I0)') 'iroot=', iroot
  stop
 end if

 target_mult = mult
 if(iroot>0 .and. mult/=xmult) target_mult = xmult
 nacte1 = nacte + 2*ras1_nacto

 call reduce_nproc_and_enlarge_mem(nproc, mem, 2, nproc1, mem1)
 call find_specified_suffix(inpname, '.inp', i)
 inpname1 = inpname(1:i-1)//'.t'
 open(newunit=fid1,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(inpname1),status='replace')

 read(fid1,'(A)') buf   ! skip %pal
 read(fid1,'(A)') buf   ! skip %maxcore
 write(fid2,'(A,I0,A)') '%pal nprocs ', nproc1, ' end'
 write(fid2,'(A,I0)') '%maxcore ', mem1

 read(fid1,'(A)') buf   ! skip '!' line
 write(fid2,'(A)',advance='no') '! NoIter'
 ! RIJK in CASSCF must be combined with CONVentional
 if(RI) write(fid2,'(A)',advance='no') ' RIJK conv '//TRIM(RIJK_bas)
 if(polar) write(fid2,'(A)',advance='no') ' def2-TZVPPD/C'
 ! TODO: automatically determine RIC_bas in this case

 ! CASCI reference
 write(fid2,'(/,A)') '%casscf'
 write(fid2,'(A,I0)') ' nel ', nacte
 write(fid2,'(A,I0)') ' norb ', nacto
 if(iroot > 0) call prt_orca_ss_cas_weight(fid2, mult, xmult, iroot)
 call prt_hard_or_crazy_casci_orca(1, fid2, hardwfn, crazywfn)
 if(RI) write(fid2,'(A)') ' TrafoStep RI'
 write(fid2,'(A)') 'end'

 ! RASCI
 write(fid2,'(A)') '%rasci'
 write(fid2,'(6(A,I0),A)') ' refs RAS (',nacte1,': ',ras1_nacto,' ',ras_nhole, &
                           ' / ',nacto,' / ',ras3_nacto,' ',ras_npart,') end'
 write(fid2,'(A,I0)') ' mult ', target_mult
 if(iroot > 1) write(fid2,'(A,I0)') ' nroots ', iroot
 if(nstate > 1) write(fid2,'(A,I0)') ' nroots ', nstate
 write(fid2,'(A)') ' CIstep CSFCI'
 !write(fid2,'(A)') ' ActOrbs NatOrbs'
 if(crazywfn) then
  write(fid2,'(A)') ' MaxIter 1000'
 else if(hardwfn) then
  write(fid2,'(A)') ' MaxIter 500'
 else
  write(fid2,'(A)') ' MaxIter 100'
 end if
 write(fid2,'(A)') 'end'

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  k = LEN_TRIM(buf)
  if(k > 0) write(fid2,'(A)') buf(1:k)
  ! there is no need to write any blank line in ORCA input
 end do ! for while

 if(polar) then
  write(fid2,'(A,/,A,/,A)') '%casresp',' MaxIter 200','end'
  write(fid2,'(A,/,A,/,A)') '%elprop',' polar true','end'
 end if

 close(fid1,status='delete')
 close(fid2)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_orca_ras_inp

! print RASCI/RASSCF keywords into a given Molpro input file
subroutine prt_molpro_ras_inp(inpname, ndb, nacte, nacto, ras1_nacto, &
                              ras3_nacto, ras_nhole, ras_npart, scf, force)
 use mr_keyword, only: RI, RIJK_bas
 implicit none
 integer :: i, k, nacto1, fid
 integer, intent(in):: ndb, nacte, nacto, ras1_nacto, ras3_nacto, ras_nhole, &
                       ras_npart
 character(len=21) :: RIJK_bas1
 character(len=80) :: sbuf
 character(len=240) :: buf, orbfile, put
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: scf, force

 orbfile = inpname
 call convert2molpro_fname(orbfile, '.a')
 nacto1 = ndb + ras1_nacto + nacto + ras3_nacto

 put = ' '
 open(newunit=fid,file=TRIM(inpname),status='old',position='append')
 BACKSPACE(fid)
 read(fid,'(A)') put

 if(put(1:4) /= '{put') then
  close(fid)
  write(6,'(/,A)') 'ERROR in subroutine prt_molpro_ras_inp: wrong content found&
                   & in the final'
  write(6,'(A)') 'line of file '//TRIM(inpname)
  stop
 end if

 rewind(fid)
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:2) == 'wf') exit
 end do ! for while

 write(fid,'(A)') '{matrop;'
 write(fid,'(A)') 'read,mo,ORB,file='//TRIM(orbfile)//';'
 write(fid,'(A)') 'save,mo,2140.2,ORBITALS}'

 write(fid,'(A)',advance='no') '{'
 if(RI) write(fid,'(A)',advance='no') 'DF-'
 write(fid,'(A)',advance='no') 'CASSCF'
 if(RI) then
  call auxbas_convert(RIJK_bas, RIJK_bas1, 2)
  write(fid,'(A)',advance='no') ',df_basis='//TRIM(RIJK_bas1)//',df_basis_exch='&
                                //TRIM(RIJK_bas1)
 end if

 write(fid,'(2(A,I0))',advance='no') ';closed,', ndb, ';occ,', nacto1

 ! Note: we need 'NoExtra' to completely turn off symmetry.
 ! Otherwise the CASCI energy is slightly different to that of other programs
 if(scf) then
  write(fid,'(A)') ';NoExtra;MaxIter,300;'
 else
  write(fid,'(A)') ';NoExtra;DONT,ORBITAL;'
 end if
 write(fid,'(A)') 'config;'

 write(sbuf,'(A,I0,A)') '(A,',ras1_nacto,'(A1,I0,A2),A)'
 k = ndb
 write(fid,fmt=TRIM(sbuf)) 'restrict,2,4',(',',i,'.1',i=k+1,k+ras1_nacto),';'

 write(sbuf,'(A,I0,A)') '(A,',ras3_nacto,'(A1,I0,A2),A)'
 k = k + ras1_nacto + nacto
 write(fid,fmt=TRIM(sbuf)) 'restrict,0,2',(',',i,'.1',i=k+1,k+ras3_nacto),';'

 write(fid,'(A)') 'natorb}'
 if(force) write(fid,'(A)') 'Forces'
 write(fid,'(A)') TRIM(put)
 close(fid)
end subroutine prt_molpro_ras_inp

! print a RASCI/RASSCF .gjf
subroutine prt_ras_gjf(gjfname, nacte, nacto, ras1_nacto, ras3_nacto, ras_nhole,&
                       ras_npart, scf, force, dkh2, polar)
 use mr_keyword, only: mem, nproc, iroot
 implicit none
 integer :: i, nacte1, nacto1, fid
 integer, intent(in) :: nacte, nacto, ras1_nacto,ras3_nacto,ras_nhole,ras_npart
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: scf, force, dkh2, polar

 nacte1 = nacte + 2*ras1_nacto
 nacto1 = ras1_nacto + nacto + ras3_nacto

 call find_specified_suffix(gjfname, '.gjf', i)
 i = INDEX(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=',nproc
 write(fid,'(6(A,I0),A)',advance='no') '#p CASSCF(',nacte1,',',nacto1, &
  ',RAS=(',ras_nhole,',',ras1_nacto,',',ras_npart,',',ras3_nacto,')'
 if(iroot > 0) write(fid,'(A,I0)',advance='no') ',qc,nroot=',iroot+1
 write(fid,'(A)',advance='no') ') chkbasis nosymm guess=read geom=allcheck'

 if(dkh2) then
  write(fid,'(A)',advance='no') ' int(nobasistransform,DKH2) iop(3/93=1)'
 else
  write(fid,'(A)',advance='no') ' int=nobasistransform'
 end if
 if(force) write(fid,'(A)',advance='no') ' force'
 if(polar) write(fid,'(A)',advance='no') ' polar(numerical,step=10)'
 !TODO: it seems that Gaussian does not recognize `step=` as electric field
 ! setting.

 if(scf) then ! RASSCF
  write(fid,'(A)') ' scf(maxcycle=300,NoVarAcc)'
 else         ! RASCI
  write(fid,'(A)') ' scf(maxcycle=-2)'
  ! To obtain CASCI NOs we need -2, since -1 only calculates the CASCI energy
 end if

 if(polar) then
  write(fid,'(/)')
 else
  write(fid,'(/,A)') '--Link1--'
  write(fid,'(A)') '%chk='//TRIM(chkname)
  write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
  write(fid,'(A,I0)') '%nprocshared=',nproc
  write(fid,'(A)',advance='no') '#p chkbasis nosymm guess(read,only,save,&
                                &NaturalOrbitals) geom=allcheck'
  if(dkh2) then
   write(fid,'(A,/)') ' int(nobasistransform,DKH2) iop(3/93=1)'
  else
   write(fid,'(A,/)') ' int=nobasistransform'
  end if
 end if

 close(fid)
end subroutine prt_ras_gjf

! read RASCI Natural Orbital Occupation Numbers (generated by unrelaxed density)
! from a given ORCA output file
subroutine read_rasci_noon_from_orca_out(outname, nacto, ras_on)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nacto
!f2py intent(in) :: nacto
 real(kind=8), intent(out) :: ras_on(nacto)
!f2py intent(out) :: ras_on
!f2py depend(nacto) :: ras_on
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
!f2py intent(in) :: outname

 ras_on = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:17) == 'Natural Orbital O') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_rasci_noon_from_orca_out: failed t&
                   &o locate'
  write(6,'(A)') '`Natural Orbital O` in file '//TRIM(outname)
  close(fid)
  stop
 end if

 do i = 1, nacto, 1
  read(fid,'(A)') buf
  call get_dpv_after_flag(buf, '=', .true., ras_on(i))
 end do ! for i

 close(fid)
end subroutine read_rasci_noon_from_orca_out

