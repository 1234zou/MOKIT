! written by jxzou at 20210305: move HF subroutines here
! updated by jxzou at 20210305: add interfaces with ORCA and PSI4 (for RI-JK HF)

! perform RHF/UHF computation using Gaussian/PySCF/PSI4/ORCA
subroutine do_hf(prt_mr_strategy)
 use mol, only: natom, atom2frag, nfrag, frag_char_mult, coor, elem, nuc, &
  charge, mult, rhf_e, uhf_e
 use mr_keyword, only: hf_prog, readuhf, readrhf, skiphf, gau_path, hf_fch, &
  ist, mo_rhf, bgchg, read_bgchg_from_gjf, gjfname, chgname, uno, vir_proj, &
  prt_strategy, gau_path, orca_path, psi4_path, frag_guess, HFonly
 implicit none
 integer :: i, system
 real(kind=8) :: ssquare = 0d0
 real(kind=8), parameter :: r_u_diff = 1d-5 ! a.u.
 character(len=24) :: data_string = ' '
 character(len=240) :: rhf_gjfname, uhf_gjfname, hf_prog_path
 logical :: eq, noiter
 logical, intent(in) :: prt_mr_strategy

 write(6,'(//,A)') 'Enter subroutine do_hf...'

 if(skiphf) then
  if(ist == 6) then
   write(6,'(A)') 'fch file will be provided. Skip the RHF/UHF step...'
   call read_natom_from_gjf(gjfname, natom)
   allocate(coor(3,natom), elem(natom), nuc(natom))
   call read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
   if(bgchg) call read_bgchg_from_gjf(.false.)
   return
  else
   write(6,'(A)') 'Provided fch(k) file. Skip the RHF/UHF step...'
   call read_natom_from_fch(hf_fch, natom)
   allocate(coor(3,natom), elem(natom), nuc(natom))
   call read_elem_and_coor_from_fch(hf_fch, natom, elem, nuc, coor, charge, mult)
  end if

  if(readuhf .and. mult==1) then
   write(6,'(A)') 'Check whether provided UHF wavefunction is equivalent to RHF...'
   call check_if_uhf_equal_rhf(hf_fch, eq)
   if(eq) then
    write(6,'(A)') 'This is actually a RHF wave function. Alpha=Beta. Switching&
                   & to ist=3.'
    i = system('fch_u2r '//TRIM(hf_fch))
    i = index(hf_fch, '.fch', back=.true.)
    hf_fch = hf_fch(1:i-1)//'_r.fch'
    readuhf = .false.; readrhf = .true.; ist = 3
    vir_proj = .true.; mo_rhf = .true. ; uno = .false.
    if(prt_mr_strategy) then
     write(6,'(A)') 'Strategy updated:'
     call prt_strategy()
    end if
   else
    write(6,'(A)') 'This seems a truly UHF wave function.'
   end if
  end if

  if(bgchg) call read_bgchg_from_gjf(.true.)
  call fdate(data_string)
  write(6,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)
  return
 end if

 call read_natom_from_gjf(gjfname, natom)
 allocate(coor(3,natom), elem(natom), nuc(natom))
 call read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)

 if(frag_guess) then
  allocate(atom2frag(natom), frag_char_mult(2,nfrag))
  call read_frag_guess_from_gjf(gjfname, natom, atom2frag, nfrag, frag_char_mult)
  if(mult == 1) write(6,'(A)') 'Fragment guess required. Only UHF will be performed.'
 else
  call check_frag_guess_in_gjf(gjfname)
 end if
 if(bgchg) call read_bgchg_from_gjf(.false.)

 i = index(gjfname, '.gjf', back=.true.)
 rhf_gjfname = gjfname(1:i-1)//'_rhf.gjf'
 uhf_gjfname = gjfname(1:i-1)//'_uhf.gjf'

 noiter = .true. ! For HF_prog=PySCF/PSI4/ORCA, we need Gaussian .fch file to
                 ! generate PySCF/PSI4/ORCA input file
 select case(TRIM(hf_prog))
 case('gaussian')
  hf_prog_path = gau_path
  noiter = .false.
 case('pyscf')
  hf_prog_path = 'python'
 case('psi4')
  hf_prog_path = psi4_path
 case('orca')
  hf_prog_path = orca_path
 case default
  write(6,'(A)') 'ERROR in subroutine do_hf: invalid HF_prog='//TRIM(hf_prog)
  write(6,'(A)') 'HF_prog can only be one of Gaussian/PySCF/PSI4/ORCA.'
  stop
 end select

 write(6,'(A)') 'HF using program '//TRIM(hf_prog)
 if(TRIM(hf_prog) /= 'pyscf') call check_exe_exist(hf_prog_path)

 ! 1) If HF_prog /= 'Gaussian', in the three 'call do_scf_and_read_e' below,
 !  only initial guess will be generated, i.e. guess(save,only).
 ! 2) Then the input file of PySCF/PSI4/ORCA will be generated from this initial
 !  guess .fch file.
 ! 3) Then the input file will be submitted to PySCF/PSI4/ORCA, respectively

 if(mult==1 .and. (.not.frag_guess)) then
  call generate_hf_gjf(rhf_gjfname, .false., noiter)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(rhf_gjfname))
  call do_scf_and_read_e(gau_path, hf_prog_path, rhf_gjfname, .false., rhf_e, ssquare)
  write(6,'(/,A,F18.8,1X,A,F7.3)') 'E(RHF) = ',rhf_e,'a.u., <S**2>=',0.0

  call generate_hf_gjf(uhf_gjfname, .true., noiter)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_gjfname))
  call do_scf_and_read_e(gau_path, hf_prog_path, uhf_gjfname, .false., uhf_e, ssquare)
  write(6,'(A,F18.8,1X,A,F7.3)')   'E(UHF) = ',uhf_e,'a.u., <S**2>=',ssquare

  if(rhf_e - uhf_e > r_u_diff) then
   write(6,'(A)') 'UHF energy is lower, choose UHF wave function.'
   ist = 1
   mo_rhf = .false.
   i = index(gjfname, '.gjf', back=.true.)
   hf_fch = gjfname(1:i-1)//'_uhf.fch'
  else
   write(6,'(A)') 'RHF/UHF energy is equal, or has little difference, choose RHF.'
   ist = 3
   vir_proj = .true.; mo_rhf = .true.; uno = .false.
   i = index(gjfname, '.gjf', back=.true.)
   hf_fch = gjfname(1:i-1)//'_rhf.fch'
  end if

 else ! only perform UHF

  call generate_hf_gjf(uhf_gjfname, .true., noiter)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_gjfname))
  call do_scf_and_read_e(gau_path, hf_prog_path, uhf_gjfname, .false., uhf_e, ssquare)
  write(6,'(/,A,F18.8,1X,A,F7.3)') 'E(UHF) = ',uhf_e,'a.u., <S**2>=',ssquare
  ist = 1
  mo_rhf = .false.
  i = index(gjfname, '.gjf', back=.true.)
  hf_fch = gjfname(1:i-1)//'_uhf.fch'
 end if

 if(prt_mr_strategy) then
  write(6,'(A)') 'Strategy updated:'
  call prt_strategy()
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)

 if(HFonly) then
  write(6,'(/,A)') 'HFonly keyword is specified. Now stop the program.'
  call fdate(data_string)
  write(6,'(/,A)') 'Normal termination of AutoMR at '//TRIM(data_string)
  stop
 end if
end subroutine do_hf

! generate a RHF/UHF .gjf file (DKH, guess=fragment can be taken into account)
subroutine generate_hf_gjf(gjfname, uhf, noiter)
 use mol, only: charge, mult, natom, nuc, elem, coor, nfrag, atom2frag, frag_char_mult
 use mr_keyword, only: mem, nproc, basis, cart, dkh2_or_x2c, frag_guess, basname,&
                       origin_gjf=>gjfname
 implicit none
 integer :: i, fid
 character(len=21) :: basis1
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: uhf, noiter
 logical :: rel, create

 create = .false.
 if(frag_guess .and. (.not.uhf)) then
  write(6,'(A)') 'ERROR in subroutine generate_hf_gjf: both frag_guess and RHF&
                 & are .True.'
  write(6,'(A)') 'Fragment guess can only be used with UHF wave function.'
  stop
 end if

 call upper(basis)
 i = index(basis, 'MA')
 if(i > 0) basis(i:i+1) = 'ma'
 i = index(basis, 'DEF')
 if(i > 0) basis(i:i+2) = 'def'
 i = index(basis, 'CC-P')
 if(i > 0) basis(i:i+3) = 'cc-p'
 i = index(basis, 'PCSSEG')
 if(i > 0) basis(i:i+5) = 'pcSseg'
 i = index(basis, 'GENECP')
 if(i > 0) basis(i:i+5) = 'genecp'
 i = index(basis, 'GEN')
 if(i > 0) basis(i:i+2) = 'gen'

 rel = .false.
 select case(TRIM(basis))
 case('DKH-def2-SV(P)','DKH-def2-SVP','DKH-def2-TZVP','DKH-def2-TZVP(-f)',&
      'DKH-def2-TZVPP','DKH-def2-QZVPP','ZORA-def2-SV(P)','ZORA-def2-SVP',&
      'ZORA-def2-TZVP','ZORA-def2-TZVP(-f)','ZORA-def2-TZVPP','ZORA-def2-QZVPP',&
      'ma-ma-DKH-def2-SV(P)','ma-DKH-def2-SVP','ma-DKH-def2-TZVP',&
      'ma-DKH-def2-TZVP(-f)','ma-DKH-def2-TZVPP','ma-DKH-def2-QZVPP',&
      'ma-ZORA-def2-SVP','ma-ZORA-def2-SV(P)','ma-ZORA-def2-TZVP',&
      'ma-ZORA-def2-TZVP(-f)','ma-ZORA-def2-TZVPP','ma-ZORA-def2-QZVPP')
  basis1 = 'gen'
  rel = .true.

  do i = 1, natom, 1
   if(nuc(i) > 36) then
    write(6,'(/,A)') 'ERROR in subroutine generate_hf_gjf: basis sets DKH-def2-&
                     & or ZORA-def2- series'
    write(6,'(A)') "have no definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i
  create = .true.

 case('ANO-RCC-MB','ANO-RCC-VDZP','ANO-RCC-VTZP','ANO-RCC-VQZP') 
  basis1 = 'gen'
  rel = .true.

  do i = 1, natom, 1
   if(nuc(i) > 96) then
    write(6,'(/,A)') 'ERROR in subroutine generate_hf_gjf: basis sets ANO-RCC-V&
                     &nZP series have no'
    write(6,'(A)') "definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i
  create = .true.

 case('pcSseg-1','pcSseg-2')
  basis1 = 'gen'

  do i = 1, natom, 1
   if(nuc(i) > 36) then ! H~Kr
    write(6,'(/,A)') 'ERROR in subroutine generate_hf_gjf: basis sets pcSseg-n&
                     & series have no'
    write(6,'(A)') "definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i
  create = .true.

 case default
  basis1 = basis
 end select

 if(rel .and. (.not. dkh2_or_x2c)) then
  write(6,'(/,A61)') REPEAT('-',61)
  write(6,'(A)') ' Warning: you are using relativistic all-electron basis set.'
  write(6,'(A)') " But you did not specify 'DKH2' or 'X2C' keyword in mokit{}."
  write(6,'(A61)') REPEAT('-',61)
 end if

 i = index(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p scf(xqc,maxcycle=512) nosymm int(nobasistransform'

 if(dkh2_or_x2c) then
  write(fid,'(A)',advance='no') ',DKH2)'
 else
  write(fid,'(A)',advance='no') ')'
 end if

 if(frag_guess) then
  if(uhf) write(fid,'(A,I0,A)',advance='no') ' guess(fragment=',nfrag,')'
  if(nfrag < 2) then
   write(fid,'(A)') 'ERROR in subroutine generate_hf_gjf: frag_guess is acitva&
                    &ted. But got nfrag < 2.'
   write(fid,'(A,I0)') 'nfrag=', nfrag
   close(fid)
   stop
  end if
  if(.not. (allocated(atom2frag) .and. allocated(frag_char_mult))) then
   write(fid,'(A)') 'ERROR in subroutine generate_hf_gjf: frag_guess is acitva&
                    &ted. But arrays atom2frag'
   write(fid,'(A)') 'and/or frag_char_mult are not allocated.'
   close(fid)
   stop
  end if
 end if

 if(uhf) then
  write(fid,'(A)',advance='no') ' UHF/'//TRIM(basis1)
  if(.not. frag_guess) then
   if(mult == 1) then
    if(noiter) then
     write(fid,'(A)',advance='no') ' guess(mix,only,save)'
    else
     write(fid,'(A)',advance='no') ' guess=mix stable=opt'
    end if
   else ! mult /= 1
    if(noiter) then
     write(fid,'(A)',advance='no') ' guess(only,save)'
    else
     write(fid,'(A)',advance='no') ' stable=opt'
    end if
   end if
  end if
 else ! RHF
  write(fid,'(A)',advance='no') ' RHF/'//TRIM(basis1)
  if(noiter) write(fid,'(A)',advance='no') ' guess(only,save)'
  if(mult /= 1) then
   write(6,'(A)') 'ERROR in subroutine generate_hf_gjf: this molecule is not sp&
                  &in singlet, but RHF is specified.'
   close(fid)
   stop
  end if
 end if

 if(cart) then
  write(fid,'(A)') ' 6D 10F'
 else
  write(fid,'(A)') ' 5D 7F'
 end if

 write(fid,'(/,A,/)') 'HF file generated by AutoMR of MOKIT'

 if(frag_guess .and. uhf) then
  write(fid,'(I0,1X,I0)',advance='no') charge, mult
  do i = 1, nfrag-1, 1
   write(fid,'(2(1X,I0))',advance='no') frag_char_mult(1,i), frag_char_mult(2,i)
  end do ! for i
  write(fid,'(2(1X,I0))') frag_char_mult(1,nfrag), frag_char_mult(2,nfrag)
  do i = 1, natom, 1
   write(fid,'(A,I0,A1,2X,3F15.8)') TRIM(elem(i))//'(fragment=',atom2frag(i),')', coor(1:3,i)
  end do ! for i
  deallocate(atom2frag, frag_char_mult)

 else
  write(fid,'(I0,1X,I0)') charge, mult
  do i = 1, natom, 1
   write(fid,'(A2,3X,3F15.8)') elem(i), coor(1:3,i)
  end do ! for i
 end if

 if(create) then
  i = index(origin_gjf, '.gjf')
  basname = origin_gjf(1:i-1)//'.bas'
  call create_basfile(basname, TRIM(basis))
 end if

 if(create .or. ((.not.create) .and. basis(1:3)=='gen')) then
  write(fid,'(/)',advance='no')
  close(fid)
  call copy_gen_basis_bas2gjf(basname, gjfname)
  open(newunit=fid,file=TRIM(gjfname),status='old',position='append')
 end if

 ! If DKH Hamiltonian is used,
 ! Gaussian default    : Gaussian function distribution
 ! Gaussian iop(3/93=1): point nuclei charge distribution
 ! GAMESS default      : point nuclei charge distribution
 ! I found that if iop(3/93=1) is used initially, SCF sometimes converges slowly,
 ! so I use a --Link1-- to add iop(3/93=1) later
 if(dkh2_or_x2c) then
  if(.not. noiter) then
   write(fid,'(/,A)') '--Link1--'
   write(fid,'(A)') '%chk='//TRIM(chkname)
   write(fid,'(A,I0,A)') '%mem=',mem,'GB'
   write(fid,'(A,I0)') '%nprocshared=', nproc
   write(fid,'(A)',advance='no') '#p scf(xqc,maxcycle=512)'
   if(uhf) then
    write(fid,'(A)',advance='no') ' UHF'
   else
    write(fid,'(A)',advance='no') ' RHF'
   end if
   write(fid,'(A)') ' chkbasis nosymm guess=read geom=allcheck iop(3/93=1)&
                    & int(nobasistransform,DKH2)'
  end if
 else
  if(frag_guess .and. (.not. noiter)) then
   write(fid,'(/,A)') '--Link1--'
   write(fid,'(A)') '%chk='//TRIM(chkname)
   write(fid,'(A,I0,A)') '%mem=',mem,'GB'
   write(fid,'(A,I0)') '%nprocshared=', nproc
   write(fid,'(A)') '#p scf(xqc,maxcycle=512) UHF chkbasis stable=opt nosymm&
                   & guess=read geom=allcheck int(nobasistransform)'
  end if
 end if

 write(fid,'(/)',advance='no')
 close(fid)
end subroutine generate_hf_gjf

! perform SCF computaton using Gaussian/PySCF/PSI4/ORCA, then read electronic
! energy and spin square
! Note: parameters {nproc, bgchg, chgname} are taken from
!  module mr_keyword. You need to initilize them before calling this subroutine.
subroutine do_scf_and_read_e(gau_path, hf_prog_path, gjfname, noiter, e, ssquare)
 use mr_keyword, only: nproc, bgchg, chgname, orca_path, psi4_path, DKH2, X2C
 use mol, only: ptchg_e
 use util_wrapper, only: formchk, bas_fch2py_wrap, fch2psi_wrap, fch2mkl_wrap,&
  mkl2gbw, gbw2mkl, mkl2fch_wrap
 implicit none
 integer :: i, hf_type, system, RENAME
 real(kind=8), intent(out) :: e, ssquare
 character(len=20) :: prog_name
 character(len=240) :: chkname, fchname, inpname, outname1, outname2, mklname
 character(len=240) :: gbwname, prpname1, prpname2, inpname1
 character(len=240), intent(in) :: gau_path, hf_prog_path, gjfname
 ! gau_path is always the path of Gaussian
 ! when Gaussian is used to compute SCF, gau_path = hf_prog_path
 ! when PySCF/PSI4/ORCA is used to compute SCF, gau_path /= hf_prog_path
 logical, intent(in) :: noiter ! skipped SCF, no iteration

 e = 0d0; ssquare = 0d0

 i = index(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'
 fchname = gjfname(1:i-1)//'.fch'
 mklname = gjfname(1:i-1)//'.mkl'
 gbwname = gjfname(1:i-1)//'.gbw'
 prpname1 = gjfname(1:i-1)//'.prop'
 prpname2 = gjfname(1:i-1)//'_property.txt'
 inpname = gjfname(1:i-1)//'.inp'
 outname2 = gjfname(1:i-1)//'.out' ! PySCF/PSI4/ORCA output file
#ifdef _WIN32
 outname1 = gjfname(1:i-1)//'.out' ! Gaussian output file under Windows
#else
 outname1 = gjfname(1:i-1)//'.log' ! Gaussian output file under Linux
#endif

 i = system(TRIM(gau_path)//' '//TRIM(gjfname))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine do_scf_and_read_e: Gaussian SCF job fai&
                   &iled.'
  write(6,'(A)') 'You can open file '//TRIM(outname1)//' and check why.'
  stop
 end if

 call delete_file('fort.7')
 if(noiter) return ! no energy read

 call formchk(chkname, fchname)
 call delete_file(chkname)

 if(TRIM(hf_prog_path) == TRIM(gau_path)) then
  ! For g09 or older, add DKH2/X2C into Route Section if needed
  if(index(gau_path,'g03')>0 .or. index(gau_path,'g09')>0) then
   if(DKH2) then
    call add_DKH2_into_fch(fchname)
   else if(X2C) then
    call add_X2C_into_fch(fchname)
   end if
  end if

  call read_hf_e_and_ss_from_gau_out(outname1, e, ssquare)
!  call delete_file(gjfname)
  return
 end if
 call delete_files(2, [gjfname, outname1])

 if(TRIM(hf_prog_path) == 'python') then
  prog_name = 'pyscf'
  i = index(gjfname, '.gjf', back=.true.)
  inpname = gjfname(1:i-1)//'.py'
 else
  i = index(hf_prog_path, '/', back=.true.)
  if(i == 0) then
   write(6,'(A)') "ERROR in subroutine do_scf_and_read_e: no '/' symbol&
                  & found in string '"//TRIM(hf_prog_path)//"'."
   stop
  end if
  prog_name = TRIM(hf_prog_path(i+1:))
 end if

 select case(prog_name)
 case('pyscf')
  call bas_fch2py_wrap(fchname, .true.)
  call read_hf_type_from_pyscf_inp(inpname, hf_type)
  call prt_hf_pyscf_inp(inpname, hf_type)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call submit_pyscf_job(inpname)
  call read_hf_e_and_ss_from_pyscf_out(outname2, hf_type, e, ssquare)
  e = e + ptchg_e

 case('psi4')
  call fch2psi_wrap(fchname, inpname)
  call read_hf_type_from_psi4_inp(inpname, hf_type)
  call prt_hf_psi4_inp(inpname, hf_type)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  call submit_psi4_job(psi4_path, inpname, nproc)
  call read_hf_e_and_ss_from_psi4_out(outname2, hf_type, e, ssquare)
  e = e + ptchg_e

  call delete_file(inpname)
  i = index(fchname, '.fch', back=.true.)
  prpname1 = fchname(1:i-1)//'2.fch'
  call copy_orb_and_den_in_fch(prpname1, fchname, .true.)

 case('orca')
  call fch2mkl_wrap(fchname, mklname)
  i = index(gjfname, '.gjf', back=.true.)
  inpname1 = gjfname(1:i-1)//'_o.inp'
  i = RENAME(TRIM(inpname1), TRIM(inpname))
  call read_hf_type_from_orca_inp(inpname, hf_type)
  call prt_hf_orca_inp(inpname, hf_type)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call mkl2gbw(mklname)
  ! mkl2gbw should be called after add_bgcharge_to_inp, since add_bgcharge_to_inp
  ! will modify both .inp and .mkl file
  call submit_orca_job(orca_path, inpname)
  call read_hf_e_and_ss_from_orca_out(outname2, hf_type, e, ssquare)
  call gbw2mkl(gbwname)
  call mkl2fch_wrap(mklname, fchname, .false.)
  call update_density_using_mo_in_fch(fchname)
  call delete_files(5, [inpname, mklname, gbwname, prpname1, prpname2])
 case default
  write(6,'(A)') 'ERROR in subroutine do_scf_and_read_e: invalid prog_name = '&
                //TRIM(prog_name)
  stop
 end select
end subroutine do_scf_and_read_e

! read HF electronic energy from a Gaussian .log/.out file
subroutine read_hf_e_and_ss_from_gau_out(logname, e, ss)
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e, ss ! HF energy and spin square
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 e = 0d0; ss = 0d0
 open(newunit=fid,file=TRIM(logname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(2:18) == 'Entering Gaussian') then
   i = -1
   exit
  end if
  if(index(buf,'SCF Done') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_gau_out: 'SCF Don&
                   &e' not found"
  write(6,'(A)') 'in file '//TRIM(logname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) e

 ! We do not read <S**2> below 'SCF Done' because when the spin is very high,
 !  the format here would become <S**2>=******* due to Fortran features.
 ! Instead, we search the 'S**2 before annihilation' below 'SCF Done'
 do i = 1, 7
  read(fid,'(A)') buf
  if(buf(2:14) == 'S**2 before a') then
   read(buf(26:),*) ss
   exit
  end if
 end do ! for i

 close(fid)
end subroutine read_hf_e_and_ss_from_gau_out

! read HF electronic energy from a PySCF .out file
subroutine read_hf_e_and_ss_from_pyscf_out(outname, wfn_type, e, ss)
 implicit none
 integer :: i, mult, fid
 integer, intent(in) :: wfn_type
 real(kind=8), intent(out) :: e, ss
 character(len=240) :: buf, inpname
 character(len=240), intent(in) :: outname

 e = 0d0; ss = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 write(fid,'(/)',advance='no')
 ! add a blank line, in case 'converged SCF e' is already in the last line

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:15) == 'converged SCF e') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_pyscf_out: no 'c&
                   &onverged SCF e' found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) e
 close(fid)

 select case(wfn_type)
 case(1,2) ! R(O)HF
  i = index(outname, '.out', back=.true.)
  inpname = outname(1:i-1)//'.py'
  call read_mult_from_pyscf_inp(inpname, mult)
  ss = DBLE((mult-1))*0.5d0
  ss = DBLE(ss*(ss+1))
 case(3)   ! UHF
  i = index(buf,'=')
  buf(i:i) = ' '
  i = index(buf,'=')
  read(buf(i+1:),*) ss
 case default
  write(6,'(A,I0)') 'ERROR in subroutine read_hf_e_and_ss_from_pyscf_out: inva&
                    &lid wfn_type=', wfn_type
  stop
 end select
end subroutine read_hf_e_and_ss_from_pyscf_out

! read HF electronic energy from a PSI4 .out file
subroutine read_hf_e_and_ss_from_psi4_out(outname, hf_type, e, ss)
 implicit none
 integer :: i, mult, fid
 integer, intent(in) :: hf_type
 real(kind=8), intent(out) :: e, ss
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 e = 0d0; ss = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')
 ! There are two HF energies in this file, we should begin from the end of file
 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(5:14) == 'HF Final E') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_psi4_out: no&
                     & 'HF Final E' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, ':')
 read(buf(i+1:),*) e

 select case(hf_type)
 case(1) ! RHF
  close(fid)
  call read_mult_from_psi4_out(outname, mult)
  ss = DBLE(mult*(mult+1))

 case(3)   ! UHF
  do while(.true.)
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   BACKSPACE(fid,iostat=i)
   if(i /= 0) exit
   read(fid,'(A)') buf
   if(buf(5:16) == 'S^2 Observed') exit
  end do ! for while

  close(fid)
  if(i /= 0) then
   write(6,'(A)') "ERROR in subroutine read_hf_e_and_ss_from_psi4_out: no&
                    & 'S^2 Observed' found in"
   write(6,'(A)') 'file '//TRIM(outname)
   stop
  end if
  i = index(buf, ':', back=.true.)
  read(buf(i+1:),*) ss

 case default
  write(6,'(A,I0)') 'ERROR in subroutine read_hf_e_and_ss_from_psi4_out:&
                      & invalid hf_type=', hf_type
  stop
 end select

end subroutine read_hf_e_and_ss_from_psi4_out

! read HF electronic energy from a ORCA .out file
subroutine read_hf_e_and_ss_from_orca_out(outname, hf_type, e, ss)
 implicit none
 integer :: i, mult, fid
 integer, intent(in) :: hf_type
 real(kind=8), intent(out) :: e, ss
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 e = 0d0; ss = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  BACKSPACE(fid,iostat=i)
  if(i /= 0) exit
  read(fid,'(A)') buf
  if(buf(1:12) == 'Total Energy') exit
  if(buf(2:12) == 'With contri') then
   i = -1
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_orca_out: no 'To&
                   &tal Energy' found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, ':')
 read(buf(i+1:),*) e

 select case(hf_type)
 case(1,2) ! R(O)HF
  close(fid)
  ! pure spin state, simply calculate the spin square
  call read_mult_from_orca_out(outname, mult)
  ss = DBLE(mult*(mult+1))

 case(3)   ! UHF
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:15) == 'Expectation val') exit
  end do ! for while

  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_orca_out: no 'Ex&
                    &pectation val' found"
   write(6,'(A)') 'in file '//TRIM(outname)
   stop
  end if
  i = index(buf, ':', back=.true.)
  read(buf(i+1:),*) ss

 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine read_hf_e_and_ss_from_orca_out: inva&
                      &lid hf_type=', hf_type
  write(6,'(A)') 'outname='//TRIM(outname)
  close(fid)
  stop
 end select
end subroutine read_hf_e_and_ss_from_orca_out

! read spin multiplicity from a PySCF .py file
subroutine read_mult_from_pyscf_inp(inpname, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == 'mol.spin') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mult_from_pyscf_inp: no 'mol.spin' &
                & found in"
  write(6,'(A)') 'file '//TRIM(inpname)
  stop
 end if

 i = index(buf,'=')
 read(buf(i+1:),*) mult ! this is No.(alpha-beta)
 mult = mult + 1
end subroutine read_mult_from_pyscf_inp

! read spin multiplicity from a PSI4 output file
subroutine read_mult_from_psi4_out(outname, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(3:14) == 'Multiplicity') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mult_from_psi4_out: no 'Multiplicit&
                 &y' found in file "//TRIM(outname)
  stop
 end if
 i = index(buf,'=',back=.true.)
 read(buf(i+1:),*) mult
end subroutine read_mult_from_psi4_out

! read spin multiplicity from an ORCA output file
subroutine read_mult_from_orca_out(outname, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:13) == 'Multiplicity') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mult_from_orca_out: no 'Multip&
                    &licity' found in file "//TRIM(outname)
  stop
 end if
 i = index(buf,'.',back=.true.)
 read(buf(i+1:),*) mult
end subroutine read_mult_from_orca_out

! print HF job of PySCF input file
subroutine prt_hf_pyscf_inp(inpname, hf_type)
 use mr_keyword, only: mem, nproc, dkh2_or_x2c
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: hf_type
 character(len=240) :: buf, inpname1, fchname
 character(len=240), intent(in) :: inpname

 i = index(inpname, '.py', back=.true.)
 fchname = inpname(1:i-1)//'.fch'
 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine prt_hf_pyscf_inp: no blank line found i&
                   &n file "//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(A)') 'from pyscf import lib'
 if(hf_type > 2) write(fid1,'(A)') 'from mokit.lib.rwwfn import update_density_&
                                   &using_mo_in_fch'
 write(fid1,'(A)') 'from mokit.lib.py2fch import py2fch'
 write(fid1,'(A,I0,A,/)') 'lib.num_threads(',nproc,')'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:8) == 'mf = scf') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine prt_hf_pyscf_inp: no 'mf = scf' found i&
                   &n file "//TRIM(inpname)
  stop
  close(fid)
  close(fid1,status='delete')
 end if

 if(dkh2_or_x2c) then
  write(fid1,'(A)') TRIM(buf)//'.x2c1e()'
 else
  write(fid1,'(A)') TRIM(buf)
 end if
 write(fid1,'(A,I0,A)') 'mf.max_memory = ',mem*1000,' # MB'

 i = 0
 do while(.true.)
  read(fid,'(A)') buf
  if(buf(1:9) == 'mf.kernel') then
   i = i + 1
   if(i == 2) exit
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 write(fid1,'(A,I0,A)') 'mf.max_memory = ',mem*1000,' # MB'
 write(fid1,'(A)') 'old_e = mf.kernel(dm0=dm)'
 ! If SCF is not converged, use the Newton method to continue
 write(fid1,'(/,A)') 'if mf.converged is False:'
 write(fid1,'(A)') '  mf = mf.newton()'
 write(fid1,'(A,/)') '  old_e = mf.kernel()'

 select case(hf_type)
 case(1,2) ! R(O)HF
  write(fid1,'(A)') '# save R(O)HF MOs into .fch file'
  write(fid1,'(A)') "py2fch('"//TRIM(fchname)//"',nbf,nif,mf.mo_coeff,'a',mf.mo&
                    &_energy,False,True)"
 case(3)   ! UHF
  ! loop to check wave function stability
  write(fid1,'(A)') 'new_e = old_e + 2e-5'
  write(fid1,'(A)') 'for i in range(0,10):'
  write(fid1,'(A)') '  mo1 = mf.stability()[0]'
  write(fid1,'(A)') '  dm1 = mf.make_rdm1(mo1, mf.mo_occ)'
  write(fid1,'(A)') '  mf = mf.newton()'
  write(fid1,'(A)') '  new_e = mf.kernel(dm0=dm1)'
  write(fid1,'(A)') '  if(abs(new_e-old_e) < 1e-5):'
  write(fid1,'(A)') '    break # cannot find lower solution'
  write(fid1,'(A)') '  old_e = new_e'

  write(fid1,'(/,A)') '# save UHF MOs into .fch file'
  write(fid1,'(A)') "py2fch('"//TRIM(fchname)//"',nbf,nif,mf.mo_coeff[0],'a',mf&
                    &.mo_energy[0],False,False)"
  write(fid1,'(A)') "py2fch('"//TRIM(fchname)//"',nbf,nif,mf.mo_coeff[1],'b',mf&
                    &.mo_energy[1],False,False)"
  write(fid1,'(A)') "update_density_using_mo_in_fch('"//TRIM(fchname)//"')"
 case default
  write(6,'(A,I0)') 'ERROR in subroutine prt_hf_pyscf_inp: invalid hf_type=',&
                     hf_type
  stop
 end select

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_hf_pyscf_inp

! print HF job of PSI4 input file
subroutine prt_hf_psi4_inp(inpname, hf_type)
 use mr_keyword, only: mem
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: hf_type
 character(len=240) :: buf, inpname1, fchname
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 i = index(inpname, '.inp', back=.true.)
 fchname = inpname(1:i-1)//'2.fch'
 ! Why append '2' in the filename: the generated .fch file will not be directly
 ! used. The utility fch_mo_copy will be called to copy alpha (and beta) MOs

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 read(fid,'(A)') buf
 write(fid1,'(A)') TRIM(buf)
 read(fid,'(A)') buf
 write(fid1,'(A,I0,A)') 'memory ',mem,' GB'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(buf(2:11) == 'guess read') exit
 end do ! for while

 ! UHF wfn stability check
 if(hf_type == 3) then
  write(fid1,'(A)') ' stability_analysis follow'
  write(fid1,'(A)') ' max_attempts 5'
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 BACKSPACE(fid1)
 write(fid1,'(A)') "scfenergy, scf_wfn = energy('scf', return_wfn=True)"
 write(fid1,'(A)') "fchk(scf_wfn,'"//TRIM(fchname)//"')"
 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_hf_psi4_inp

! print HF job of ORCA input file
subroutine prt_hf_orca_inp(inpname, hf_type)
 use mr_keyword, only: mem, nproc, dkh2_or_x2c
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: hf_type
 character(len=240) :: buf, inpname1
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(inpname1),status='replace')

 read(fid,'(A)') buf
 write(fid1,'(A,I0,A)') '%pal nprocs ',nproc,' end'
 read(fid,'(A)') buf
 write(fid1,'(A,I0)') '%maxcore ',CEILING(DBLE(1000*mem)/DBLE(nproc))

 if(dkh2_or_x2c .and. hf_type==3) then
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:3) == 'end') exit
   write(fid1,'(A)') TRIM(buf)
  end do ! for while
  write(fid1,'(A)') ' OneCenter True'
  write(fid1,'(A)') 'end'
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
  if(buf(1:4) == '%scf') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine prt_hf_orca_inp: no '%scf' found in file&
                 & "//TRIM(inpname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 if(hf_type == 3) then
  write(fid1,'(A)') ' STABPerform true'
  write(fid1,'(A)') ' STABRestartUHFifUnstable true'
  write(fid1,'(A)') ' STABMaxIter 500'
  write(fid1,'(A)') ' STABDTol 1e-5'
  write(fid1,'(A)') ' STABRTol 1e-5'
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(inpname1), TRIM(inpname))
end subroutine prt_hf_orca_inp

! read hf_type from a PySCF .inp file
subroutine read_hf_type_from_pyscf_inp(inpname, hf_type)
 implicit none
 integer :: i, fid
 integer, intent(out) :: hf_type ! 0/1/2/3 for undetermined/RHF/ROHF/UHF
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:2) == 'mf') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_hf_type_from_pyscf_inp: no 'mf' fou&
                 &nd in file "//TRIM(inpname)
  stop
 end if

 if(index(buf,'RHF') > 0) then
  hf_type = 1
 else if(index(buf,'ROHF') > 0) then
  hf_type = 2
 else if(index(buf,'HF') > 0) then
  hf_type = 3
 else 
  write(6,'(A)') "ERROR in subroutine read_hf_type_from_pyscf_inp: no 'HF' fou&
                 &nd in file "//TRIM(inpname)
  close(fid)
  stop
 end if
end subroutine read_hf_type_from_pyscf_inp

! read hf_type from a PSI4 .inp file
subroutine read_hf_type_from_psi4_inp(inpname, hf_type)
 implicit none
 integer :: i, fid
 integer, intent(out) :: hf_type ! 0/1/2/3 for undetermined/RHF/ROHF/UHF
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:10) == 'reference') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_hf_type_from_psi4_inp: no 'refe&
                    &rence' found in file "//TRIM(inpname)
  stop
 end if

 if(index(buf,'rhf') > 0) then
  hf_type = 1
 else if(index(buf,'rohf') > 0) then
  hf_type = 2
 else if(index(buf,'uhf') > 0) then
  hf_type = 3
 else if(index(buf,'hf') > 0) then
  hf_type = 0
 else 
  write(6,'(A)') "ERROR in subroutine read_hf_type_from_psi4_inp: no 'hf'&
                   & found in file "//TRIM(inpname)
  close(fid)
  stop
 end if
end subroutine read_hf_type_from_psi4_inp

! read hf_type from an ORCA .inp file
subroutine read_hf_type_from_orca_inp(inpname, hf_type)
 implicit none
 integer :: i, fid
 integer, intent(out) :: hf_type ! 0/1/2/3 for undetermined/RHF/ROHF/UHF
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname

 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '!') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_hf_type_from_orca_inp: no '!' found&
                 & in file "//TRIM(inpname)
  stop
 end if

 if(index(buf,'RHF') > 0) then
  hf_type = 1
 else if(index(buf,'ROHF') > 0) then
  hf_type = 2
 else if(index(buf,'UHF') > 0) then
  hf_type = 3
 else if(index(buf,'HF') > 0) then
  hf_type = 0
 else 
  write(6,'(A)') "ERROR in subroutine read_hf_type_from_orca_inp: no 'HF' found&
                & in file "//TRIM(inpname)
  close(fid)
  stop
 end if
end subroutine read_hf_type_from_orca_inp

! Check whether fragment information can be found in Cartesian coordinate
!  section when 'guess(fragment=N)' is not specified.
! If fragment information is found, stop and print error
subroutine check_frag_guess_in_gjf(gjfname)
 implicit none
 integer :: i, fid, nblank
 character(len=240) :: buf
 character(len=240), intent(in) :: gjfname

 open(newunit=fid,file=TRIM(gjfname),status='old',position='rewind')
 nblank = 0
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(LEN_TRIM(buf) == 0) nblank = nblank + 1
  if(nblank == 2) exit
 end do ! for while

 if(i /= 0) then
  close(fid)
  return
 end if

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 close(fid)
 call lower(buf)
 if(index(buf,'fragment') > 0) then
  write(6,'(/,A)') 'ERROR in subroutine check_frag_guess_in_gjf: fragment infor&
                   &mation found in'
  write(6,'(A)') "coordinate section. But 'guess(fragment=N)' keyword not found&
                & in file "//TRIM(gjfname)//'.'
  write(6,'(A)') 'It seems that you forgot to write necessary keywords.'
  stop
 end if
end subroutine check_frag_guess_in_gjf

