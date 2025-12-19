! written by jxzou at 20210305: move HF subroutines here
! updated by jxzou at 20210305: add interfaces with ORCA and PSI4 (for RI-JK HF)

! perform RHF/UHF computation using Gaussian/PySCF/PSI4/ORCA
subroutine do_hf(prt_mr_strategy)
 use mol, only: natom, atom2frag, nfrag, frag_char_mult, coor, elem, nuc, &
  charge, mult, rhf_e, uhf_e, uhf_ssquare
 use mr_keyword, only: hf_prog, readuhf, readrhf, skiphf, gau_path, hf_fch, &
  ist, mo_rhf, bgchg, read_bgchg_from_gjf, gjfname, chgname, uno, vir_proj, &
  prt_strategy, gau_path, orca_path, psi4_path, frag_guess
 use util_wrapper, only: fch_u2r_wrap
 implicit none
 integer :: i, SYSTEM
 real(kind=8) :: ssquare = 0d0
 real(kind=8), parameter :: r_u_diff = 1d-4 ! a.u.
 character(len=24) :: data_string = ' '
 character(len=240) :: proname, rhf_inp, uhf_inp, hf_prog_path
 logical :: eq, noiter
 logical, intent(in) :: prt_mr_strategy

 write(6,'(//,A)') 'Enter subroutine do_hf...'

 call find_specified_suffix(gjfname, '.gjf', i)
 proname = gjfname(1:i-1)

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
    call fch_u2r_wrap(hf_fch)
    i = INDEX(hf_fch, '.fch', back=.true.)
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

 rhf_inp = TRIM(proname)//'_rhf.gjf'
 uhf_inp = TRIM(proname)//'_uhf.gjf'

 ! For HF_prog=Gaussian/PySCF/ORCA, we simply perform HF calculations using them.
 ! For HF_prog=PSI4, we need Gaussian .fch file to generate PSI4 input file.
 select case(TRIM(hf_prog))
 case('gaussian')
  hf_prog_path = gau_path
  noiter = .false.
 case('pyscf')
  hf_prog_path = 'python'
  noiter = .false.
  rhf_inp = TRIM(proname)//'_rhf.py'
  uhf_inp = TRIM(proname)//'_uhf.py'
 case('psi4')
  hf_prog_path = psi4_path
  noiter = .false.
 case('orca')
  hf_prog_path = orca_path
  noiter = .false.
  rhf_inp = TRIM(proname)//'_rhf.inp'
  uhf_inp = TRIM(proname)//'_uhf.inp'
 case default
  write(6,'(/,A)') 'ERROR in subroutine do_hf: invalid HF_prog='//TRIM(hf_prog)
  write(6,'(A)') 'HF_prog can only be one of Gaussian/PySCF/PSI4/ORCA.'
  stop
 end select

 write(6,'(A)') 'HF using program '//TRIM(hf_prog)
 if(TRIM(hf_prog) /= 'pyscf') call check_exe_exist(hf_prog_path)

 ! 1) If HF_prog is ORCA/PSI4, in the three 'call do_scf_and_read_e' below,
 !  only initial guess will be generated, i.e. guess(save,only).
 ! 2) Then the input file of PSI4/ORCA will be generated from this initial
 !  guess .fch file.
 ! 3) Then the input file will be submitted to PSI4/ORCA, respectively.

 if(mult==1 .and. (.not.frag_guess) .and. (ist/=7)) then
  call gen_hf_inp(hf_prog, rhf_inp, .false., noiter)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(rhf_inp))
  call do_scf_and_read_e(gau_path, hf_prog_path, rhf_inp, .false., rhf_e, ssquare)
  write(6,'(/,A,F18.8,1X,A,F7.3)') 'E(RHF) = ',rhf_e,'a.u., <S**2>=',0.0

  call gen_hf_inp(hf_prog, uhf_inp, .true., noiter)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_inp))
  call do_scf_and_read_e(gau_path, hf_prog_path, uhf_inp, .false., uhf_e, ssquare)
  write(6,'(A,F18.8,1X,A,F7.3)') 'E(UHF) = ',uhf_e,'a.u., <S**2>=',ssquare
  uhf_ssquare = ssquare

  if(rhf_e - uhf_e > r_u_diff) then
   write(6,'(A)') 'UHF energy is lower, choose UHF wave function.'
   ist = 1
   mo_rhf = .false.
   hf_fch = TRIM(proname)//'_uhf.fch'
  else
   write(6,'(A)') 'RHF/UHF energies are too close, choose RHF.'
   ist = 3
   vir_proj = .true.; mo_rhf = .true.; uno = .false.
   hf_fch = TRIM(proname)//'_rhf.fch'
  end if

 else ! only perform UHF

  call gen_hf_inp(hf_prog, uhf_inp, .true., noiter)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_inp))
  call do_scf_and_read_e(gau_path, hf_prog_path, uhf_inp, .false., uhf_e, ssquare)
  write(6,'(/,A,F18.8,1X,A,F7.3)') 'E(UHF) = ',uhf_e,'a.u., <S**2>=',ssquare
  uhf_ssquare = ssquare
  if(ist /= 7) ist = 1
  mo_rhf = .false.
  hf_fch = TRIM(proname)//'_uhf.fch'
 end if

 if(prt_mr_strategy .and. (ist /= 7)) then
  write(6,'(A)') 'Strategy updated:'
  call prt_strategy()
 end if

 call fdate(data_string)
 write(6,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)
end subroutine do_hf

subroutine process_basis_set(iprog, basis, basis1, dkh2_or_x2c, natom, nuc, &
                             elem, create)
 implicit none
 integer :: i
 integer, intent(in) :: iprog, natom
 integer, intent(in) :: nuc(natom)
 integer, parameter :: x2c_idx(38) = (/19,20,37,38,55,56,57,58,59,60,61,62,63,&
  64,65,66,67,68,69,70,71,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,&
  103/)
 character(len=2), intent(in) :: elem(natom)
 character(len=21), intent(inout) :: basis
 character(len=21), intent(out) :: basis1
 logical :: rel
 logical, intent(in) :: dkh2_or_x2c
 logical, intent(out) :: create

 rel = .false.; create = .false.; basis1 = ' '
 call upper(basis)

 i = INDEX(basis, 'CC-P')
 if(i > 0) basis(i:i+3) = 'cc-p'
 i = INDEX(basis, 'DEF')
 if(i > 0) basis(i:i+2) = 'def'
 i = INDEX(basis, 'MA-')
 if(i > 0) basis(i:i+2) = 'ma-'
 i = INDEX(basis, 'PCSSEG')
 if(i > 0) basis(i:i+5) = 'pcSseg'
 i = INDEX(basis, 'GENECP')
 if(i > 0) basis(i:i+5) = 'genecp'
 i = INDEX(basis, 'GEN')
 if(i > 0) basis(i:i+2) = 'gen'
 i = INDEX(basis, 'X2C-')
 if(i > 0) basis(i:i+3) = 'x2c-'
 i = INDEX(basis, 'ALL')
 if(i > 0) basis(i:i+2) = 'all'
 i = INDEX(basis, '(D,P)')
 if(i > 0) basis(i:i+4) = '(d,p)'
 i = INDEX(basis, '(-F)')
 if(i > 0) basis(i:i+3) = '(-f)'

 select case(TRIM(basis))
 case('ma-def2SV(P)','ma-def2SVP','ma-def2TZVP(-f)','ma-def2TZVP', &
      'ma-def2TZVPP','ma-def2QZVP','ma-def2QZVPP')
  if(iprog == 3) then
   create = .false.
   i = INDEX(basis, 'def2')
   basis1 = basis(1:i+3)//'-'//TRIM(basis(i+4:))
  else
   create = .true.
   if( ANY(nuc>36) ) then
    basis1 = 'genecp'
   else
    basis1 = 'gen'
   end if
  end if

  do i = 1, natom, 1
   if(nuc(i) > 86) then
    write(6,'(/,A)') 'ERROR in subroutine process_basis_set: basis sets ma-def2&
                     & series have'
    write(6,'(A)') "no definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i

 case('DKH-def2SV(P)','DKH-def2SVP','DKH-def2TZVP','DKH-def2TZVP(-f)', &
      'DKH-def2TZVPP','DKH-def2QZVPP','ma-DKH-def2SV(P)','ma-DKH-def2SVP', &
      'ma-DKH-def2TZVP','ma-DKH-def2TZVP(-f)','ma-DKH-def2TZVPP', &
      'ma-DKH-def2QZVPP')
  rel = .true.
  if(iprog == 3) then
   create = .false.
   i = INDEX(basis, 'def2')
   basis1 = basis(1:i+3)//'-'//TRIM(basis(i+4:))
  else
   basis1 = 'gen'; create = .true.
  end if

  do i = 1, natom, 1
   if(nuc(i) > 36) then
    write(6,'(/,A)') 'ERROR in subroutine process_basis_set: basis sets DKH-def&
                     &2- series have'
    write(6,'(A)') "no definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i

 case('x2c-SV(P)all','x2c-SVPall','x2c-TZVPall','x2c-TZVPPall','x2c-QZVPall',&
      'x2c-QZVPPall')
  rel = .true.
  if(iprog == 3) then
   basis1 = basis; create = .false.
  else
   basis1 = 'gen'; create = .true.
  end if

  do i = 1, natom, 1
   if(nuc(i) > 86) then
    write(6,'(/,A)') 'ERROR in subroutine process_basis_set: basis sets x2c- se&
                     &ries have no'
    write(6,'(A)') "definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i

 case('ANO-RCC-MB','ANO-RCC-VDZP','ANO-RCC-VTZP','ANO-RCC-VQZP') 
  basis1 = 'gen'; rel = .true.; create = .true.

  do i = 1, natom, 1
   if(nuc(i) > 96) then
    write(6,'(/,A)') 'ERROR in subroutine process_basis_set: basis sets ANO-RCC&
                     &-VnZP series have'
    write(6,'(A)') "no definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i

 case('ANO-R1','ANO-R2','ANO-R3') 
  basis1 = 'gen'; rel = .true.; create = .true.

  do i = 1, natom, 1
   if(nuc(i) > 86) then
    write(6,'(/,A)') 'ERROR in subroutine process_basis_set: basis sets ANO-Rn &
                     &series have no'
    write(6,'(A)') "definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i

 case('cc-pVDZ-X2C','cc-pVTZ-X2C') 
  basis1 = 'gen'; rel = .true.; create = .true.

  do i = 1, natom, 1
   if(ANY(x2c_idx /= nuc(i))) then
    write(6,'(/,A)') 'ERROR in subroutine process_basis_set: basis sets cc-pVnZ&
                     &-X2C series have'
    write(6,'(A)') "no definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i

 case('pcSseg-1','pcSseg-2')
  if(iprog == 3) then
   basis1 = basis; create = .false.
  else
   basis1 = 'gen'; create = .true.
  end if

  do i = 1, natom, 1
   if(nuc(i) > 36) then ! H~Kr
    write(6,'(/,A)') 'ERROR in subroutine process_basis_set: basis sets pcSseg-&
                     &n series have'
    write(6,'(A)') "no definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i

 case('cc-pVDZ-F12','cc-pVTZ-F12')
  if(iprog == 3) then
   basis1 = basis; create = .false.
  else
   basis1 = 'gen'; create = .true.
  end if

  do i = 1, natom, 1
   if(nuc(i) > 18) then ! H~Kr
    write(6,'(/,A)') 'ERROR in subroutine process_basis_set: basis sets cc-pVnZ&
                     &-F12 series have'
    write(6,'(A)') "no definition on element '"//TRIM(elem(i))//"'."
    stop
   end if
  end do ! for i

 case default
  basis1 = basis
 end select

 if(rel .and. (.not. dkh2_or_x2c)) then
  write(6,'(/,A)') REPEAT('-',61)
  write(6,'(A)') ' Warning: you are using relativistic all-electron basis set.'
  write(6,'(A)') " But you did not specify 'DKH2' or 'X2C' keyword in mokit{}."
  write(6,'(A)') REPEAT('-',61)
 end if
end subroutine process_basis_set

! generate RHF/UHF input files
subroutine gen_hf_inp(hf_prog, inpname, uhf, noiter)
 implicit none
 character(len=10), intent(in) :: hf_prog
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: uhf, noiter

 select case(TRIM(hf_prog))
 case('gaussian')
  call gen_hf_gjf(inpname, uhf, noiter)
 case('pyscf')
  call gen_hf_pyscf_inp(inpname, uhf)
 case('orca')
  call gen_hf_orca_inp(inpname, uhf)
 case default
  write(6,'(/,A)') 'ERROR in subroutine gen_hf_inp: invalid HF_prog='//&
                   TRIM(hf_prog)
  stop
 end select
end subroutine gen_hf_inp

! generate an RHF/UHF .gjf file (DKH2, guess=fragment can be taken into account)
subroutine gen_hf_gjf(gjfname, uhf, noiter)
 use mol, only: charge, mult, natom, nuc, elem, coor, nfrag, atom2frag, &
  frag_char_mult
 use mr_keyword, only: mem, nproc, basis, cart, dkh2_or_x2c, frag_guess, basname,&
  origin_gjf=>gjfname
 implicit none
 integer :: i, fid
 character(len=21) :: basis1
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: uhf, noiter
 logical :: create, def2_ecp

 if(frag_guess .and. (.not.uhf)) then
  write(6,'(/,A)') 'ERROR in subroutine gen_hf_gjf: both frag_guess and RHF are&
                   & .True.'
  write(6,'(A)') 'Fragment guess can only be used with UHF wave function.'
  stop
 end if

 call process_basis_set(1, basis, basis1, dkh2_or_x2c, natom, nuc, elem, create)
 def2_ecp = .false.
 if(basis(1:6)=='ma-def' .and. ANY(nuc>36)) def2_ecp = .true.

 i = INDEX(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=', mem, 'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p scf(xqc,maxcycle=512) nosymm int(nobasistran&
                               &sform'
 if(dkh2_or_x2c) write(fid,'(A)',advance='no') ',DKH2'
 write(fid,'(A)',advance='no') ')'

 if(frag_guess) then
  if(uhf) write(fid,'(A,I0,A)',advance='no') ' guess(fragment=',nfrag,')'
  if(nfrag < 2) then
   write(fid,'(/,A)') 'ERROR in subroutine gen_hf_gjf: frag_guess is activated.&
                      & But got nfrag<2.'
   write(fid,'(A,I0)') 'nfrag=', nfrag
   close(fid)
   stop
  end if
  if(.not. (allocated(atom2frag) .and. allocated(frag_char_mult))) then
   write(fid,'(/,A)') 'ERROR in subroutine gen_hf_gjf: frag_guess is activated.&
                      & But arrays atom2frag'
   write(fid,'(A)') 'and/or frag_char_mult are not allocated.'
   close(fid)
   stop
  end if
 end if

 if(uhf) then ! UHF
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
 else         ! RHF/ROHF
  write(fid,'(A)',advance='no') ' R'
  if(mult > 1) write(fid,'(A)',advance='no') 'O'
  write(fid,'(A)',advance='no') 'HF/'//TRIM(basis1)
  if(noiter) write(fid,'(A)',advance='no') ' guess(only,save)'
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
   write(fid,'(A,I0,A1,2X,3F15.8)') TRIM(elem(i))//'(fragment=', atom2frag(i), &
                                    ')', coor(1:3,i)
  end do ! for i
  deallocate(atom2frag, frag_char_mult)
 else
  write(fid,'(I0,1X,I0)') charge, mult
  do i = 1, natom, 1
   write(fid,'(A2,3X,3F15.8)') elem(i), coor(1:3,i)
  end do ! for i
 end if

 if(create) then
  i = INDEX(origin_gjf, '.gjf')
  basname = origin_gjf(1:i-1)//'.bas'
  call create_basfile(basname, TRIM(basis), def2_ecp)
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
    if(mult > 1) then
      write(fid,'(A)',advance='no') ' ROHF'
     else
      write(fid,'(A)',advance='no') ' RHF'
    end if
   end if
   write(fid,'(A)') ' chkbasis nosymm guess=read geom=allcheck iop(3/93=1) int(&
                    &nobasistransform,DKH2)'
  end if
 else
  if(frag_guess .and. (.not. noiter)) then
   write(fid,'(/,A)') '--Link1--'
   write(fid,'(A)') '%chk='//TRIM(chkname)
   write(fid,'(A,I0,A)') '%mem=',mem,'GB'
   write(fid,'(A,I0)') '%nprocshared=', nproc
   write(fid,'(A)') '#p UHF chkbasis nosymm int=nobasistransform guess=read geo&
                    &m=allcheck scf(xqc,maxcycle=512) stable=opt'
  end if
 end if

 write(fid,'(/)',advance='no')
 close(fid)
end subroutine gen_hf_gjf

! generate a PySCF RHF/UHF .py file
subroutine gen_hf_pyscf_inp(pyname, uhf)
 use mol, only: charge, mult, natom, nuc, elem, coor
 use mr_keyword, only: mem, nproc, basis, cart, dkh2_or_x2c, frag_guess, basname,&
  gjfname
 implicit none
 integer :: i, fid
 character(len=21) :: basis1
 character(len=240) :: fchname
 character(len=240), intent(in) :: pyname
 logical, intent(in) :: uhf
 logical :: create

 if(frag_guess) then
  write(6,'(/,A)') 'ERROR in subroutine gen_hf_pyscf_inp: frag_guess is current&
                   &ly not supported'
  write(6,'(A)') 'when HF_prog=PySCF. You can use HF_prog=Gaussian.'
  stop
 end if

 call process_basis_set(2, basis, basis1, dkh2_or_x2c, natom, nuc, elem, create)
 call find_specified_suffix(gjfname, '.gjf', i)
 basname = gjfname(1:i-1)//'.bas'
 if(create) call create_basfile(basname, TRIM(basis), .false.)

 call find_specified_suffix(pyname, '.py', i)
 fchname = pyname(1:i-1)//'.fch'

 open(newunit=fid,file=TRIM(pyname),status='replace')
 write(fid,'(A)') 'from pyscf import gto, scf, lib'
 write(fid,'(A)') 'from mokit.lib.py2fch_direct import fchk'
 write(fid,'(A)') 'from mokit.lib.rwwfn import get_nmo_from_ao_ovlp'
 if(uhf) then
  write(fid,'(A)') 'from mokit.lib.rwwfn import update_density_using_mo_in_fch'
  write(fid,'(A)') 'from mokit.lib.stability import uhf_stable_opt_internal'
  if(mult == 1) then
   write(fid,'(A)') 'from mokit.lib.rwwfn import gen_no_from_dm_and_ao_ovlp, \'
   write(fid,'(A)') '                            calc_dm_using_mo_and_on'
   write(fid,'(A)') 'from math import sqrt'
   write(fid,'(A)') 'import numpy as np'
  end if
 else
  if(mult > 1) then
   write(fid,'(A)') 'from mokit.lib.stability import rohf_stable_opt_internal'
  end if
 end if
 write(fid,'(/,A,I0,A)') 'lib.num_threads(', nproc, ')'

 write(fid,'(A)') 'mol = gto.M()'
 write(fid,'(A,I0,A)') '# ', natom, ' atom(s)'
 write(fid,'(A)') "mol.atom = '''"
 do i = 1, natom, 1
  write(fid,'(A2,1X,3(1X,F17.8))') elem(i), coor(:,i)
 end do ! for i
 write(fid,'(A)') "'''"
 write(fid,'(A)') "mol.basis = '"//TRIM(basis1)//"'"

 if(basis1(1:4)=='def2' .and. ANY(nuc>36)) then
  write(fid,'(A)') 'mol.ecp = {'
  do i = 1, natom, 1
   if(nuc(i) > 36) then
    write(fid,'(A)') " '"//TRIM(elem(i))//"': '"//TRIM(basis1)//"',"
   end if
  end do ! for i
  write(fid,'(A)') '}'
 end if

 write(fid,'(A)') '# Remember to check the charge and spin'
 write(fid,'(A,I0)') 'mol.charge = ', charge
 write(fid,'(A,I0)') 'mol.spin = ', mult-1
 if(cart) write(fid,'(A)') 'mol.cart = True'
 write(fid,'(A)') 'mol.verbose = 4'
 write(fid,'(A)') 'mol.build(parse_arg=False)'

 write(fid,'(/,A)',advance='no') 'mf = scf.'
 if(uhf) then
  write(fid,'(A)',advance='no') 'UHF'
 else
  if(mult > 1) then
   write(fid,'(A)',advance='no') 'ROHF'
  else
   write(fid,'(A)',advance='no') 'RHF'
  end if
 end if

 if(dkh2_or_x2c) then
  write(fid,'(A)') '(mol).x2c1e()'
 else
  write(fid,'(A)') '(mol)'
 end if
 write(fid,'(A,I0,A)') 'mf.max_memory = ',mem*1000,' # MB'
 write(fid,'(A)') 'mf.max_cycle = 128'

 ! check if there is basis set linear dependency
 write(fid,'(A)') "S = mol.intor_symmetric('int1e_ovlp')"
 write(fid,'(A)') 'nbf = S.shape[0]'
 write(fid,'(A)') 'nif = get_nmo_from_ao_ovlp(nbf, S)'
 write(fid,'(A)') 'if(nif < nbf):'
 write(fid,'(A)') '  mf2 = mf.copy()'
 write(fid,'(A)') '  mf = scf.remove_linear_dep_(mf2, threshold=1e-6, lindep=1e-6)'

 if(uhf .and. mult==1) then
  ! For singlet UHF, construct broken symm initial guess
  write(fid,'(A)') 'dm_a, dm_b = mf.get_init_guess()'
  write(fid,'(A)') 'occ_a, mo_a = gen_no_from_dm_and_ao_ovlp(nbf, nif, dm_a, S)'
  write(fid,'(A)') 'ndb = np.count_nonzero(occ_a > 0.5)'
  write(fid,'(A)') 'occ_a[0:ndb] = 1.0'
  write(fid,'(A)') 'occ_a[ndb:] = 0.0'
  write(fid,'(A)') 'homo = mo_a[:,ndb-1].copy()'
  write(fid,'(A)') 'lumo = mo_a[:,ndb].copy()'
  write(fid,'(A)') 'mo_a[:,ndb-1] = lumo.copy()'
  write(fid,'(A)') 'mo_a[:,ndb]   = homo.copy()'
  write(fid,'(A)') 'dm_b = calc_dm_using_mo_and_on(nbf, nif, mo_a, occ_a)'
  write(fid,'(A)') 'dm = (dm_a, dm_b)'
  write(fid,'(A)') 'mf.kernel(dm0=dm)'
 else
  write(fid,'(A)') 'mf.kernel()'
 end if

 ! If normal SCF is unconverged, use the Newton method to continue
 write(fid,'(/,A)') 'if mf.converged is False:'
 write(fid,'(A)')   '  mf = mf.newton()'
 write(fid,'(A,/)') '  mf.kernel()'

 if(uhf) then       ! UHF
  write(fid,'(A)') '# stable=opt'
  write(fid,'(A)') 'mf = uhf_stable_opt_internal(mf)'
  write(fid,'(/,A)') '# save UHF MOs into .fch file'
  write(fid,'(A)') "uhf_fch = '"//TRIM(fchname)//"'"
  write(fid,'(A)') 'fchk(mf, uhf_fch)'
  write(fid,'(A)') 'update_density_using_mo_in_fch(uhf_fch)'
 else
  if(mult > 1) then ! ROHF
   write(fid,'(A)') '# stable=opt'
   write(fid,'(A)') 'mf = rohf_stable_opt_internal(mf)'
   write(fid,'(/,A)') '# save ROHF MOs into .fch file'
   write(fid,'(A)') "rohf_fch = '"//TRIM(fchname)//"'"
   write(fid,'(A)') 'fchk(mf, rohf_fch)'
   write(fid,'(A)') 'update_density_using_mo_in_fch(rohf_fch)'
  else              ! RHF
   write(fid,'(A)') 'if mf.converged is False:'
   write(fid,'(A)') "  raise OSError('PySCF R(O)HF job failed.')"
   write(fid,'(/,A)') '# save R(O)HF MOs into .fch file'
   write(fid,'(A)') "fchk(mf, '"//TRIM(fchname)//"', density=True)"
  end if
 end if

 close(fid)
 if(TRIM(basis1) == 'gen') call add_gen_bas2pyscf_inp(basname,pyname,natom,elem)
end subroutine gen_hf_pyscf_inp

! generate an ORCA RHF/UHF .inp file
subroutine gen_hf_orca_inp(inpname, uhf)
 use mol, only: charge, mult, natom, nuc, elem, coor
 use mr_keyword, only: mem, nproc, basis, cart, DKH2, X2C, dkh2_or_x2c, RI, &
  frag_guess
 implicit none
 integer :: i, fid
 character(len=21) :: basis1
 character(len=240), intent(in) :: inpname
 logical, intent(in) :: uhf
 logical :: create

 if(frag_guess) then
  write(6,'(/,A)') 'ERROR in subroutine gen_hf_orca_inp: frag_guess is currentl&
                   &y not supported'
  write(6,'(A)') 'when HF_prog=ORCA. Please use the default HF_prog=Gaussian.'
  stop
 end if

 if(cart) then
  write(6,'(/,A)') 'ERROR in subroutine gen_hf_orca_inp: Cartesian-type basis f&
                   &unctions (6D 10F)'
  write(6,'(A)') 'are not supported in ORCA. Please use spherical harmonic type&
                 & basis functions.'
  stop
 end if

 call process_basis_set(3, basis, basis1, dkh2_or_x2c, natom, nuc, elem, create)
 if(basis1(1:3) == 'gen') then
  write(6,'(/,A)') 'ERROR in subroutine gen_hf_orca_inp: gen/genecp is not supp&
                   &orted in this'
  write(6,'(A)') 'functionality. Please use an ORCA built-in basis set.'
  stop
 end if

 i = INDEX(basis1, 'def2')
 if(i > 0) then
  if(basis1(i+4:i+4)/='-') basis1 = basis1(1:i+3)//'-'//TRIM(basis(i+4:))
 end if

 open(newunit=fid,file=TRIM(inpname),status='replace')
 write(fid,'(A,I0,A)') '%pal nprocs ', nproc, ' end'
 write(fid,'(A,I0)') '%maxcore ', FLOOR(1d3*DBLE(mem)/DBLE(nproc))

 if(uhf) then
  write(fid,'(A)',advance='no') '! UHF '
 else
  write(fid,'(A)',advance='no') '! RHF '
 end if
 write(fid,'(A)',advance='no') TRIM(basis1)
 if(RI) then
  if(DKH2) then
   write(fid,'(A)',advance='no') ' RIJCOSX SARC/J defgrid3'
  else if(X2C) then
   write(fid,'(A)',advance='no') ' RIJCOSX x2c/J defgrid3'
  else
   select case(basis1(1:6))
   case('cc-pVT','cc-pVQ','cc-pV5','aug-cc')
    write(fid,'(A)',advance='no') ' RIJK '//TRIM(basis1)//'/JK'
   case default
    write(fid,'(A)',advance='no') ' RIJK def2/JK'
   end select
  end if
 else
  write(fid,'(A)',advance='no') ' noRI'
 end if
 write(fid,'(A)') ' VeryTightSCF'

 write(fid,'(A)') '%scf'
 write(fid,'(A)') ' Thresh 1e-12'
 write(fid,'(A)') ' Tcut 1e-14'
 write(fid,'(A)') ' MaxIter 512'
 write(fid,'(A)') ' sthresh 1e-6'
 if(uhf .and. (.not.dkh2_or_x2c)) then
  ! There is no need to use BrokenSym/FlipSpin for singlet UHF, since a .gbw
  ! file includes broken symmetry MOs will be generated from RHF calculation
  ! result. In fact, 'BrokenSym 1,1' + 'stable=opt' functionality seems not
  ! correct in ORCA 6.0.0.
  ! ORCA stable=opt has many features not implemented when DKH2/X2C is invoked.
  ! So stability check is not done here.
  write(fid,'(A)') ' STABPerform true'
  write(fid,'(A)') ' STABRestartUHFifUnstable true'
  write(fid,'(A)') ' STABMaxIter 500'
  write(fid,'(A)') ' STABDTol 1e-5'
  write(fid,'(A)') ' STABRTol 1e-5'
 end if
 write(fid,'(A)') 'end'

 if(dkh2_or_x2c) then
  write(fid,'(A)') '%rel'
  if(DKH2) then
   write(fid,'(A)') ' method DKH'
   write(fid,'(A)') ' order 2'
  else if(X2C) then
   write(fid,'(A)') ' method X2C'
  end if
  write(fid,'(A)') 'end'
 end if

 write(fid,'(A,I0,1X,I0)') '* xyz ', charge, mult
 do i = 1, natom, 1
  write(fid,'(A2,1X,3(1X,F17.8))') elem(i), coor(:,i)
 end do ! for i
 write(fid,'(A)') '*'
 close(fid)
end subroutine gen_hf_orca_inp

! perform SCF computation using Gaussian/PySCF/PSI4/ORCA, then read electronic
! energy and spin square
! Note: parameters {nproc, bgchg, chgname} are taken from
!  module mr_keyword. You need to initialize them before calling this subroutine.
subroutine do_scf_and_read_e(gau_path, hf_prog_path, inpname, noiter, e, ssquare)
 use mr_keyword, only: nproc, bgchg, chgname, orca_path, psi4_path, DKH2, X2C
 use mol, only: ptchg_e
 use util_wrapper, only: formchk, bas_fch2py_wrap, fch2psi_wrap, fch2mkl_wrap,&
  mkl2gbw, gbw2mkl, mkl2fch_wrap
 implicit none
 integer :: i, irel, hf_type, SYSTEM, RENAME
 real(kind=8), intent(out) :: e, ssquare
 character(len=20) :: prog_name
 character(len=240) :: proname, chkname, fchname, outname, mklname, gbwname, &
  prpname1, prpname2, prpname3, inpname1, denf1, denf2
 character(len=240), intent(in) :: gau_path, hf_prog_path, inpname
 ! gau_path is always the path of Gaussian
 ! when Gaussian is used to compute SCF, gau_path = hf_prog_path
 ! when PySCF/PSI4/ORCA is used to compute SCF, gau_path /= hf_prog_path
 logical, intent(in) :: noiter ! skipped SCF, no iteration

 e = 0d0; ssquare = 0d0
 call find_specified_suffix(inpname, '.', i)
 proname = inpname(1:i-1)
 outname = inpname(1:i-1)//'.out'
 fchname  = TRIM(proname)//'.fch'
 i = LEN_TRIM(hf_prog_path)

 if(hf_prog_path(1:6) == 'python') then
  if(noiter) then
   write(6,'(/,A)') 'ERROR in subroutine do_scf_and_read_e: internal inconsiste&
                    &ncy.'
   stop
  end if
  call read_hf_type_from_pyscf_inp(inpname, hf_type)
  call submit_pyscf_job(inpname, .false.)
  call read_hf_e_and_ss_from_pyscf_out(outname, hf_type, e, ssquare)
  e = e + ptchg_e
  return
 else if(hf_prog_path(i-3:i) == 'orca') then
  if(noiter) then
   write(6,'(/,A)') 'ERROR in subroutine do_scf_and_read_e: internal inconsiste&
                    &ncy.'
   stop
  end if
  mklname  = TRIM(proname)//'.mkl'
  gbwname  = TRIM(proname)//'.gbw'
  prpname1 = TRIM(proname)//'.prop'
  prpname2 = TRIM(proname)//'_property.txt'
  prpname3 = TRIM(proname)//'.property.txt'
  inpname1 = TRIM(proname)//'_o.inp'
  denf1    = TRIM(proname)//'.densities'
  denf2    = TRIM(proname)//'.densitiesinfo'
  call read_hf_type_from_orca_inp(inpname, hf_type)
  if(hf_type == 3) call gen_brokensym_guess_orca(mklname)
  call submit_orca_job(orca_path, inpname, .true., .false., .false.)
  call read_scf_e_and_ss_from_orca_out(outname, hf_type, e, ssquare)
  e = e + ptchg_e
  call gbw2mkl(gbwname)
  irel = -1
  if(DKH2) irel = 2
  if(X2C) irel = -3
  call mkl2fch_wrap(mklname=mklname,fchname=fchname,irel=irel)
  call update_density_using_mo_in_fch(fchname)
  call delete_files(7,[inpname,gbwname,prpname1,prpname2,prpname3,denf1,denf2])
  if(hf_type == 3) call delete_file(TRIM(mklname))
  return
 end if

 call submit_gau_job(gau_path, inpname, .false.)
 if(noiter) return ! no energy to read

 chkname = TRIM(proname)//'.chk'
#ifdef _WIN32
 outname = TRIM(proname)//'.out' ! Gaussian output file under Windows
#else
 outname = TRIM(proname)//'.log' ! Gaussian output file under Linux
#endif

 call formchk(chkname)
 call delete_file(chkname)

 if(TRIM(hf_prog_path) == TRIM(gau_path)) then
  call simplify_fch(fchname)

  ! For g09 or older, add DKH2/X2C into Route Section if needed
  if(INDEX(gau_path,'g03')>0 .or. INDEX(gau_path,'g09')>0) then
   if(DKH2) then
    call add_DKH2_into_fch(fchname)
   else if(X2C) then
    call add_X2C_into_fch(fchname)
   end if
  end if

  call read_hf_e_and_ss_from_gau_log(outname, e, ssquare)
!  call delete_file(gjfname)
  return
 end if
 call delete_files(2, [inpname, outname])

 call find_specified_suffix(hf_prog_path, '/', i)
 prog_name = TRIM(hf_prog_path(i+1:))
 if(X2C) call add_X2C_into_fch(fchname)

 select case(prog_name)
 case('psi4')
  call fch2psi_wrap(fchname, inpname)
  call read_hf_type_from_psi4_inp(inpname, hf_type)
  call prt_hf_psi4_inp(inpname, hf_type)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  call submit_psi4_job(psi4_path, inpname, nproc)
  call read_hf_e_and_ss_from_psi4_out(outname, hf_type, e, ssquare)
  e = e + ptchg_e

  call delete_file(inpname)
  i = INDEX(fchname, '.fch', back=.true.)
  prpname1 = fchname(1:i-1)//'2.fch'
  call copy_orb_and_den_in_fch(prpname1, fchname, .true.)

 case('orca')
  mklname  = TRIM(proname)//'.mkl'
  gbwname  = TRIM(proname)//'.gbw'
  prpname1 = TRIM(proname)//'.prop'
  prpname2 = TRIM(proname)//'_property.txt'
  inpname1 = TRIM(proname)//'_o.inp'
  call fch2mkl_wrap(fchname, mklname)
  i = RENAME(TRIM(inpname1), TRIM(inpname))

  call read_hf_type_from_orca_inp(inpname, hf_type)
  !call prt_hf_orca_inp(inpname, hf_type)
  if(bgchg) i = SYSTEM('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))
  call mkl2gbw(mklname)
  ! mkl2gbw should be called after add_bgcharge_to_inp, since add_bgcharge_to_inp
  ! will modify both .inp and .mkl file
  call submit_orca_job(orca_path, inpname, .true., .false., .false.)
  call read_scf_e_and_ss_from_orca_out(outname, hf_type, e, ssquare)
  call gbw2mkl(gbwname)
  call mkl2fch_wrap(mklname, fchname)
  call update_density_using_mo_in_fch(fchname)
  call delete_files(5, [inpname, mklname, gbwname, prpname1, prpname2])

 case default
  write(6,'(/,A)') 'ERROR in subroutine do_scf_and_read_e: invalid prog_name='&
                  //TRIM(prog_name)
  stop
 end select
end subroutine do_scf_and_read_e

! read HF electronic energy from a Gaussian .log/.out file
subroutine read_hf_e_and_ss_from_gau_log(logname, e, ss)
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
  if(INDEX(buf,'SCF Done') > 0) exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_gau_log: 'SCF Don&
                   &e' not found"
  write(6,'(A)') 'in file '//TRIM(logname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
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
end subroutine read_hf_e_and_ss_from_gau_log

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
                   &onverged SCF e'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, '=')
 read(buf(i+1:),*) e
 close(fid)

 select case(wfn_type)
 case(1,2) ! R(O)HF
  i = INDEX(outname, '.out', back=.true.)
  inpname = outname(1:i-1)//'.py'
  call read_mult_from_pyscf_inp(inpname, mult)
  ss = DBLE((mult-1))*0.5d0
  ss = DBLE(ss*(ss+1))
 case(3)   ! UHF
  i = INDEX(buf,'=')
  buf(i:i) = ' '
  i = INDEX(buf,'=')
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
  write(6,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_psi4_out: no 'HF &
                   &Final E' found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
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
   write(6,'(A)') "ERROR in subroutine read_hf_e_and_ss_from_psi4_out: no 'S^2 &
                  &Observed' found"
   write(6,'(A)') 'in file '//TRIM(outname)
   stop
  end if
  i = INDEX(buf, ':', back=.true.)
  read(buf(i+1:),*) ss

 case default
  write(6,'(A,I0)') 'ERROR in subroutine read_hf_e_and_ss_from_psi4_out: invali&
                    &d hf_type=', hf_type
  stop
 end select

end subroutine read_hf_e_and_ss_from_psi4_out

! read spin square <S^2> from a specified ORCA output file
subroutine read_spin_square_from_orca_out(outname, hf_type, ss)
 implicit none
 integer :: i, mult, fid
 integer, intent(in) :: hf_type
 real(kind=8), intent(out) :: ss
 real(kind=8), external :: mult2ssquare
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 ss = 0d0
 select case(hf_type)
 case(1,2) ! R(O)HF
  ! pure spin state, simply calculate the spin square
  call read_mult_from_orca_out(outname, mult)
  ss = mult2ssquare(mult)

 case(3)   ! UHF
  open(newunit=fid,file=TRIM(outname),status='old',position='rewind')
  do while(.true.)
   read(fid,'(A)',iostat=i) buf
   if(i /= 0) exit
   if(buf(1:15) == 'Expectation val') exit
  end do ! for while

  close(fid)
  if(i /= 0) then
   write(6,'(/,A)') "ERROR in subroutine read_spin_square_from_orca_out: no 'Ex&
                    &pectation val'"
   write(6,'(A)') 'found in file '//TRIM(outname)
   stop
  end if
  i = INDEX(buf, ':', back=.true.)
  read(buf(i+1:),*) ss

 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine read_spin_square_from_orca_out: inva&
                      &lid hf_type=', hf_type
  write(6,'(A)') 'outname='//TRIM(outname)
  stop
 end select
end subroutine read_spin_square_from_orca_out

! read HF energy from an ORCA .out file
subroutine read_hf_e_and_ss_from_orca_out(outname, hf_type, e, ss)
 implicit none
 integer :: i, fid
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
  if(buf(1:15) == 'Total Energy   ') exit
  if(buf(2:12) == 'With contri') then
   i = -1
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_orca_out: no 'To&
                   &tal Energy' found"
  write(6,'(A)') 'in file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = INDEX(buf, ':')
 read(buf(i+1:),*) e
 close(fid)
 call read_spin_square_from_orca_out(outname, hf_type, ss)
end subroutine read_hf_e_and_ss_from_orca_out

! read KS-DFT energy from an ORCA .out file
subroutine read_dft_e_and_ss_from_orca_out(outname, hf_type, e, ss)
 implicit none
 integer :: i, fid
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
  if(buf(1:15) == 'FINAL SINGLE PO') exit
  if(buf(2:12) == 'With contri') then
   i = -1
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_dft_e_and_ss_from_orca_out: no 'FI&
                   &NAL SINGLE PO'"
  write(6,'(A)') 'found in file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(buf(26:),*) e
 close(fid)
 call read_spin_square_from_orca_out(outname, hf_type, ss)
end subroutine read_dft_e_and_ss_from_orca_out

! check whether DFT or HF is used in an ORCA output file
subroutine check_dft_in_orca_out(outname, dft)
 implicit none
 integer :: i, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: outname
 logical, intent(out) :: dft

 dft = .false.
 open(newunit=fid,file=TRIM(outname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12) == 'SCF SETTINGS') exit
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine check_dft_in_orca_out: 'SCF SETTINGS' n&
                   &ot found in"
  write(6,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 read(fid,'(A)') buf 
 close(fid)
 if(buf(46:48) == 'DFT') dft = .true.
end subroutine check_dft_in_orca_out

! Read SCF energy from an ORCA .out file. DFT may have extra terms like D3
! correction, so we have to use two different subroutines.
subroutine read_scf_e_and_ss_from_orca_out(outname, hf_type, e, ss)
 implicit none
 integer, intent(in) :: hf_type
 real(kind=8), intent(out) :: e, ss
 character(len=240), intent(in) :: outname
 logical :: dft

 call check_dft_in_orca_out(outname, dft)
 if(dft) then
  call read_dft_e_and_ss_from_orca_out(outname, hf_type, e, ss)
 else
  call read_hf_e_and_ss_from_orca_out(outname, hf_type, e, ss)
 end if
end subroutine read_scf_e_and_ss_from_orca_out

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

 i = INDEX(buf,'=')
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
 i = INDEX(buf,'=',back=.true.)
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
 i = INDEX(buf,'.',back=.true.)
 read(buf(i+1:),*) mult
end subroutine read_mult_from_orca_out

! print HF job of PSI4 input file
subroutine prt_hf_psi4_inp(inpname, hf_type)
 use mr_keyword, only: mem
 implicit none
 integer :: i, fid, fid1, RENAME
 integer, intent(in) :: hf_type
 character(len=240) :: buf, inpname1, fchname
 character(len=240), intent(in) :: inpname

 inpname1 = TRIM(inpname)//'.t'
 i = INDEX(inpname, '.inp', back=.true.)
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
  write(fid1,'(A)') ' max_attempts 8'
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
  write(6,'(/,A)') "ERROR in subroutine read_hf_type_from_pyscf_inp: no 'mf' fo&
                   &und in file "//TRIM(inpname)
  stop
 end if

 if(INDEX(buf,'RHF') > 0) then
  hf_type = 1
 else if(INDEX(buf,'ROHF') > 0) then
  hf_type = 2
 else if(INDEX(buf,'HF') > 0) then
  hf_type = 3
 else 
  write(6,'(/,A)') "ERROR in subroutine read_hf_type_from_pyscf_inp: no 'HF' fo&
                   &und in file "//TRIM(inpname)
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
  write(6,'(/,A)') "ERROR in subroutine read_hf_type_from_psi4_inp: no 'referen&
                   &ce' found in file"
  write(6,'(A)') TRIM(inpname)
  stop
 end if

 if(INDEX(buf,'rhf') > 0) then
  hf_type = 1
 else if(INDEX(buf,'rohf') > 0) then
  hf_type = 2
 else if(INDEX(buf,'uhf') > 0) then
  hf_type = 3
 else if(INDEX(buf,'hf') > 0) then
  hf_type = 0
 else 
  write(6,'(/,A)') "ERROR in subroutine read_hf_type_from_psi4_inp: no 'hf' fou&
                   &nd in file "//TRIM(inpname)
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
  write(6,'(/,A)') "ERROR in subroutine read_hf_type_from_orca_inp: no '!' foun&
                  &d in file "//TRIM(inpname)
  stop
 end if

 if(INDEX(buf,'RHF') > 0) then
  hf_type = 1
 else if(INDEX(buf,'ROHF') > 0) then
  hf_type = 2
 else if(INDEX(buf,'UHF') > 0) then
  hf_type = 3
 else if(INDEX(buf,'HF') > 0) then
  hf_type = 0
 else 
  write(6,'(/,A)') "ERROR in subroutine read_hf_type_from_orca_inp: no 'HF' fou&
                   &nd in file "//TRIM(inpname)
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

 if(INDEX(buf,'fragment') > 0) then
  write(6,'(/,A)') 'ERROR in subroutine check_frag_guess_in_gjf: fragment infor&
                   &mation found in'
  write(6,'(A)') "coordinate section. But 'guess(fragment=N)' keyword not found&
                & in file "//TRIM(gjfname)//'.'
  write(6,'(A)') 'It seems that you forgot to write necessary keywords.'
  stop
 end if
end subroutine check_frag_guess_in_gjf

subroutine add_gen_bas2pyscf_inp(basname, pyname, natom, elem)
 implicit none
 integer :: i, j, fid, fid1, fid2, RENAME
 integer, intent(in) :: natom
 character(len=2), intent(in) :: elem(natom)
 character(len=3) :: str3
 character(len=4) :: str4
 character(len=300) :: buf
 character(len=240) :: pyname1, bas_copy
 character(len=240), intent(in) :: basname, pyname

 str3 = ' '; buf = ' '; bas_copy = ' '
 call find_specified_suffix(pyname, '.py', i)
 pyname1 = pyname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(pyname1),status='replace')

 do while(.true.)
  read(fid,'(A)') buf
  if(LEN_TRIM(buf) == 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while
 write(fid1,'(A,/)') 'from pyscf.gto.basis.parse_gaussian import load'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:17) == "mol.basis = 'gen'") exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine add_gen_bas2pyscf_inp: no ""mol.basis = &
                   &'gen'"" found"
  write(6,'(A)') 'in file '//TRIM(pyname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 write(fid1,'(A)') 'mol.basis = {'
 open(newunit=fid2,file=TRIM(basname),status='old',position='rewind')
 read(fid2,'(A)') buf

 if(buf(1:2) == '@/') then
  i = LEN_TRIM(buf)
  j = INDEX(buf(1:i-2), '/', back=.true.)
  bas_copy = buf(j+1:i-2)
  call sys_copy_file(buf(2:i-2), TRIM(bas_copy), .false.)
  call del_hyphen_for_elem_in_basfile(bas_copy)
  do i = 1, natom, 1
   if(i > 1) then
    if(ANY(elem(1:i-1) == elem(i))) cycle
   end if
   write(fid1,'(A)') " '"//TRIM(elem(i))//"': load('"//TRIM(bas_copy)//"','"//&
                     TRIM(elem(i))//"'),"
  end do ! for i
 else
  rewind(fid2)
  do while(.true.)
   read(fid2,*,iostat=i) str3
   if(i /= 0) exit
   if(LEN_TRIM(str3) == 0) exit
   read(fid2,'(A)') buf
   if(buf(1:2) == '@/') then
    i = LEN_TRIM(buf)
    j = INDEX(buf(1:i-2), '/', back=.true.)
    bas_copy = buf(j+1:i-2)
    call sys_copy_file(buf(2:i-2), TRIM(bas_copy), .false.)
    call del_hyphen_for_elem_in_basfile(bas_copy)
    write(fid1,'(A)') " '"//TRIM(str3(2:3))//"': load('"//TRIM(bas_copy)//&
                      "','"//TRIM(str3(2:3))//"'),"
   else
    write(fid1,'(A)') " '"//TRIM(str3(2:3))//"': '"//TRIM(buf)//"',"
   end if
   read(fid2,'(A)') str4
  end do ! for while
 end if

 close(fid2)
 write(fid1,'(A)') '}'

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
end subroutine add_gen_bas2pyscf_inp

! find the version of ORCA
subroutine find_orca_ver(orca_path, iver)
 implicit none
 integer :: i, fid, SYSTEM
 integer, intent(out) :: iver ! 4/5/6
 character(len=20) :: txtname
 character(len=240) :: buf
 character(len=240), intent(in) :: orca_path

 iver = 0 ! initialization
 call get_a_random_int(i)
 write(txtname,'(I0,A)') i, '.txt'
 i = SYSTEM(TRIM(orca_path)//' -v >'//TRIM(txtname)//" 2>&1")

 open(newunit=fid,file=TRIM(txtname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(26:40) == 'Program Version') exit
 end do ! for while

 close(fid,status='delete')
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine find_orca_ver: no 'Program Version' found."
  write(6,'(A)') 'Please check whether ORCA is installed correctly.'
  stop
 end if

 i = INDEX(buf,'.')
 read(buf(41:i-1),*) iver

 if(iver < 4) then
  write(6,'(/,A)') 'ERROR in subroutine find_orca_ver: ORCA older than 4.0 is u&
                   &nsupported.'
  stop
 end if
end subroutine find_orca_ver

! Generate a singlet UHF broken symmetry initial guess. The input .mkl should be
! xxx_uhf.mkl, in order to find the corresponding xxx_rhf.mkl and xxx_uhf.inp.
subroutine gen_brokensym_guess_orca(uhf_mkl)
 use util_wrapper, only: mkl2gbw
 implicit none
 integer :: i, mult, RENAME
 character(len=240) :: uhf_inp, rhf_mkl
 character(len=240), intent(in) :: uhf_mkl

 i = LEN_TRIM(uhf_mkl)
 if(uhf_mkl(i-7:i) /= '_uhf.mkl') return

 uhf_inp = uhf_mkl(1:i-8)//'_uhf.inp'
 call read_mult_from_orca_inp(uhf_inp, mult)
 if(mult /= 1) return ! only singlet UHF is considered here

 rhf_mkl = uhf_mkl(1:i-8)//'_rhf.mkl'
 call require_file_exist(rhf_mkl)

 i = RENAME(TRIM(rhf_mkl), TRIM(uhf_mkl))
 call mkl_r2u(uhf_mkl, .true.)
 call mkl2gbw(uhf_mkl)
end subroutine gen_brokensym_guess_orca

