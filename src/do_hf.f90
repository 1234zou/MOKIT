! written by jxzou at 20210305: move HF subroutines here
! updated by jxzou at 20210305: add interfaces with ORCA and PSI4 (for RI-JK HF)

! perform RHF/UHF computation using Gaussian/PySCF/PSI4/ORCA
subroutine do_hf()
 use print_id, only: iout
 use mol, only: natom, atom2frag, nfrag, frag_char_mult, coor, elem, nuc, &
  charge, mult, rhf_e, uhf_e
 use mr_keyword, only: hf_prog, readuhf, readrhf, skiphf, mem, nproc, basis, &
  cart, gau_path, hf_fch, ist, mo_rhf, bgchg, read_bgchg_from_gjf, gjfname, &
  chgname, uno, dkh2_or_x2c, vir_proj, prt_strategy, gau_path, orca_path, &
  psi4_path, frag_guess
 use util_wrapper, only: formchk
 implicit none
 integer :: i, system
 real(kind=8) :: ssquare = 0d0
 character(len=24) :: data_string = ' '
 character(len=240) :: rhf_gjfname, uhf_gjfname, chkname, hf_prog_path
 logical :: eq

 write(iout,'(//,A)') 'Enter subroutine do_hf...'

 if(skiphf) then
  write(iout,'(A)') 'Provided .fch(k) file. Skip the RHF/UHF step...'
  call read_natom_from_fch(hf_fch, natom)
  allocate(coor(3,natom), elem(natom), nuc(natom))
  call read_elem_and_coor_from_fch(hf_fch, natom, elem, nuc, coor, charge, mult)

  if(readuhf .and. mult==1) then
   write(iout,'(A)') 'Check whether provided UHF is equivalent to RHF...'
   call check_if_uhf_equal_rhf(hf_fch, eq)
   if(eq) then
    write(iout,'(A)') 'This is actually a RHF wave function. Alpha=Beta.&
                     & Switching to ist=3.'
    i = system('fch_u2r '//TRIM(hf_fch))
    i = index(hf_fch, '.fch', back=.true.)
    hf_fch = hf_fch(1:i-1)//'_r.fch'
    readuhf = .false.; readrhf = .true.; ist = 3
    vir_proj = .true.; mo_rhf = .true. ; uno = .false.
    write(iout,'(A)') 'Strategy updated:'
    call prt_strategy()
   else
    write(iout,'(A)') 'This seems a truly UHF wave function.'
   end if
  end if

  if(bgchg) call read_bgchg_from_gjf(.true.)
  call fdate(data_string)
  write(iout,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)
  return
 end if

 call read_natom_from_gjf(gjfname, natom)
 allocate(coor(3,natom), elem(natom), nuc(natom))
 call read_elem_and_coor_from_gjf(gjfname, natom, elem, nuc, coor, charge, mult)
 if(frag_guess) then
  allocate(atom2frag(natom), frag_char_mult(2,nfrag))
  call read_frag_guess_from_gjf(gjfname, natom, atom2frag, nfrag, frag_char_mult)
 end if
 if(bgchg) call read_bgchg_from_gjf(.false.)

 i = index(gjfname, '.gjf', back=.true.)
 rhf_gjfname = gjfname(1:i-1)//'_rhf.gjf'
 uhf_gjfname = gjfname(1:i-1)//'_uhf.gjf'

 select case(TRIM(hf_prog))
 case('gaussian')
  hf_prog_path = gau_path
 case('pyscf')
  hf_prog_path = 'python'
 case('psi4')
  hf_prog_path = psi4_path
 case('orca')
  hf_prog_path = orca_path
 case default
  write(iout,'(A)') 'ERROR in subroutine do_hf: invalid HF_prog='//TRIM(hf_prog)
  write(iout,'(A)') 'HF_prog can only be one of Gaussian/PySCF/PSI4/ORCA.'
  stop
 end select

 if(TRIM(hf_prog) /= 'pyscf') call check_exe_exist(hf_prog_path)

 ! 1) If HF_prog /= 'Gaussian', in the three 'call do_scf_and_read_e' below,
 !  only initial guess will be generated, i.e. guess(save,only).
 ! 2) Then the input file of PySCF/PSI4/ORCA will be generated from this initial
 !  guess .fch file.
 ! 3) Then the input file will be submitted to PySCF/PSI4/ORCA, respectively

 if(mult == 1) then ! singlet, perform RHF and UHF
  call generate_hf_gjf(rhf_gjfname, .false.)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(rhf_gjfname))
  call do_scf_and_read_e(gau_path, hf_prog_path, rhf_gjfname, .false., rhf_e, ssquare)

  call generate_hf_gjf(uhf_gjfname, .true.)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_gjfname))
  call do_scf_and_read_e(gau_path, hf_prog_path, uhf_gjfname, .false., uhf_e, ssquare)

 else               ! not singlet, only perform UHF
  call generate_hf_gjf(uhf_gjfname, .true.)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(uhf_gjfname))
  call do_scf_and_read_e(gau_path, hf_prog_path, uhf_gjfname, .false., uhf_e, ssquare)
 end if

 if(mult == 1) then
  write(iout,'(/,A,F18.8,1X,A,F7.3)') 'E(RHF) = ',rhf_e,'a.u., <S**2>=',0.0
  write(iout,'(A,F18.8,1X,A,F7.3)')   'E(UHF) = ',uhf_e,'a.u., <S**2>=',ssquare
 else
  write(iout,'(/,A,F18.8,1X,A,F7.3)') 'E(UHF) = ',uhf_e,'a.u., <S**2>=',ssquare
 end if

 if(rhf_e - uhf_e > 1D-4) then
  write(iout,'(A)') 'UHF energy is lower, choose UHF wave function.'
  ist = 1
  mo_rhf = .false.
  i = index(gjfname, '.gjf', back=.true.)
  chkname = gjfname(1:i-1)//'_uhf.chk'
  hf_fch = gjfname(1:i-1)//'_uhf.fch'
 else
  write(iout,'(A)') 'RHF/UHF energy is equal, or has little difference, choose RHF.'
  ist = 3
  vir_proj = .true.; mo_rhf = .true.; uno = .false.
  i = index(gjfname, '.gjf', back=.true.)
  chkname = gjfname(1:i-1)//'_rhf.chk'
  hf_fch = gjfname(1:i-1)//'_rhf.fch'
 end if

 call formchk(chkname, hf_fch)

 write(iout,'(A)') 'Strategy updated:'
 call prt_strategy()

 call fdate(data_string)
 write(iout,'(A)') 'Leave subroutine do_hf at '//TRIM(data_string)
 return
end subroutine do_hf

! generate a RHF/UHF .gjf file (DKH, guess=fragment can be taken into account)
subroutine generate_hf_gjf(gjfname, uhf)
 use mol, only: charge, mult, natom, elem, coor, nfrag, atom2frag, frag_char_mult
 use mr_keyword, only: iout, mem, nproc, basis, cart, dkh2_or_x2c, frag_guess,&
  mokit_root
 implicit none
 integer :: i, fid
 character(len=21) :: basis1
 character(len=240) :: chkname
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: uhf
 logical :: rel

 call upper(basis)
 i = index(basis, 'DEF')
 if(i > 0) basis(i:i+2) = 'def'
 i = index(basis, 'MA')
 if(i > 0) basis(i:i+1) = 'ma'

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
 case default
  basis1 = basis
 end select

 i = index(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A)') '%chk='//TRIM(chkname)
 write(fid,'(A,I0,A)') '%mem=',mem,'GB'
 write(fid,'(A,I0)') '%nprocshared=', nproc
 write(fid,'(A)',advance='no') '#p scf(xqc,maxcycle=128) nosymm int(nobasistransform'

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
  if(frag_guess) then
   write(fid,'(A)',advance='no') ' guess=mix'
  else
   if(mult == 1) write(fid,'(A)',advance='no') ' guess=mix'
   write(fid,'(A)',advance='no') ' stable=opt'
  end if
 else
  write(fid,'(A)',advance='no') ' RHF/'//TRIM(basis1)
  if(mult /= 1) then
   write(iout,'(A)') 'ERROR in subroutine generate_hf_gjf: this molecule is&
                    & not spin singlet, but RHF is specified.'
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

 if(rel .and. TRIM(basis1)=='gen') then
  write(fid,'(/,A,/)') '@'//TRIM(mokit_root)//'/basis/'//TRIM(basis)
 end if

 ! If DKH Hamiltonian is used,
 ! Gaussian default    : Gaussian function distribution
 ! Gaussian iop(3/93=1): point nuclei charge distribution
 ! GAMESS default      : point nuclei charge distribution
 ! I found that if iop(3/93=1) is used initially, SCF sometimes converges slowly,
 ! so I use a --Link1-- to add iop(3/93=1) later
 if(dkh2_or_x2c) then
  write(fid,'(/,A)') '--Link1--'
  write(fid,'(A)') '%chk='//TRIM(chkname)
  write(fid,'(A,I0,A)') '%mem=',mem,'GB'
  write(fid,'(A,I0)') '%nprocshared=', nproc
  write(fid,'(A)',advance='no') '#p scf(xqc,maxcycle=128)'
  if(uhf) then
   write(fid,'(A)',advance='no') ' UHF stable=opt'
  else
   write(fid,'(A)',advance='no') ' RHF'
  end if
  write(fid,'(A)') ' chkbasis nosymm guess=read geom=allcheck iop(3/93=1)&
                   & int(nobasistransform,DKH2)'
 else
  if(frag_guess .and. uhf) then
   write(fid,'(/,A)') '--Link1--'
   write(fid,'(A)') '%chk='//TRIM(chkname)
   write(fid,'(A,I0,A)') '%mem=',mem,'GB'
   write(fid,'(A,I0)') '%nprocshared=', nproc
   write(fid,'(A)') '#p scf(xqc,maxcycle=128) UHF chkbasis stable=opt nosymm&
                   & guess=read geom=allcheck int(nobasistransform)'
  end if
 end if

 write(fid,'(/)',advance='no')
 close(fid)
 return
end subroutine generate_hf_gjf

! perform SCF computaton using Gaussian/PySCF/PSI4/ORCA, then read electronic
! energy and spin square
! Note: parameters {mem, nproc, RI, RIJK_bas, bgchg, chgname} are taken from
!  module mr_keyword. You need to initilize them before calling this subroutine.
subroutine do_scf_and_read_e(gau_path, hf_prog_path, gjfname, noiter, e, ssquare)
 use print_id, only: iout
 use mr_keyword, only: mem, nproc, RI, RIJK_bas, bgchg, chgname
 use util_wrapper, only: formchk
 implicit none
 integer :: i, j, system, RENAME
 real(kind=8), intent(out) :: e, ssquare
 character(len=20) :: prog_name
 character(len=240) :: buf, chkname, fchname, inpname, outname1, outname2
 character(len=240), intent(in) :: gau_path, hf_prog_path, gjfname
 ! gau_path is always the path of Gaussian
 ! when Gaussian is used to compute SCF, gau_path = hf_prog_path
 ! when PySCF/PSI4/ORCA is used to compute SCF, gau_path /= hf_prog_path
 logical, intent(in) :: noiter ! skipped SCF, no iteration

 e = 0d0; ssquare = 0d0

 i = index(gjfname, '.gjf', back=.true.)
 chkname = gjfname(1:i-1)//'.chk'
 fchname = gjfname(1:i-1)//'.fch'
 inpname = gjfname(1:i-1)//'.inp'
 outname2 = gjfname(1:i-1)//'.out' ! PySCF/PSI4/ORCA output file
#ifdef _WIN32
 outname1 = gjfname(1:i-1)//'.out' ! Gaussian output file under Windows
#else
 outname1 = gjfname(1:i-1)//'.log' ! Gaussian output file under Linux
#endif

 i = system(TRIM(gau_path)//' '//TRIM(gjfname))
 if(i /= 0) then
  write(iout,'(/,A)') 'ERROR in subroutine do_scf_and_read_e: Gaussian SCF job failed.'
  write(iout,'(A)') 'You can open file '//TRIM(outname1)//' and check why.'
  stop
 end if

 if(noiter) return
 if(TRIM(hf_prog_path) == TRIM(gau_path)) then
  call read_hf_e_and_ss_from_gau_out(outname1, e, ssquare)
  return
 end if

 i = index(hf_prog_path, '/', back=.true.)
 prog_name = TRIM(hf_prog_path(i+1:))

 select case(prog_name)
 case('pyscf')
  call formchk(chkname, fchname)
  i = system('bas_fch2py '//TRIM(fchname))

  buf = 'python '//TRIM(inpname)//' >'//TRIM(outname2)//" 2>&1"
  write(iout,'(A)') '$'//TRIM(buf)
  i = system(TRIM(buf))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_scf_and_read_e: PySCF SCF job failed.'
   write(iout,'(A)') 'You can open file '//TRIM(outname2)//' and see why.'
   stop
  end if

  call read_hf_e_and_ss_from_pyscf_out(outname2, e, ssquare)

 case('psi4')
  call formchk(chkname, fchname)
  i = system('fch2psi '//TRIM(fchname))

  call modify_memory_in_psi4_inp(inpname, mem)
  call add_RIJK_bas_into_psi4_inp(inpname, RIJK_bas)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  i = LEN_TRIM(outname1)
  if(outname1(i-3:i) == '.out') then
   j = RENAME(TRIM(outname1), outname1(1:i-3)//'out1')
  end if

  write(buf,'(A,I0)') TRIM(hf_prog_path)//' '//TRIM(inpname)//' '//TRIM(outname2)//' -n ',nproc
  write(iout,'(A)') '$'//TRIM(buf)

  i = system(TRIM(buf))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_scf_and_read_e: PSI4 SCF job failed.'
   write(iout,'(A)') 'You can open file '//TRIM(outname2)//' and see why.'
   stop
  end if
  call read_hf_e_and_ss_from_psi4_out(outname2, e, ssquare)

 case('orca')
  call formchk(chkname, fchname)
  i = system('fch2mkl '//TRIM(fchname))

  call modify_memory_in_orca_inp(inpname, mem)
  call add_RIJK_bas_into_orca_inp(inpname, RIJK_bas)
  if(bgchg) i = system('add_bgcharge_to_inp '//TRIM(chgname)//' '//TRIM(inpname))

  i = LEN_TRIM(outname1)
  if(outname1(i-3:i) == '.out') then
   j = RENAME(TRIM(outname1), outname1(1:i-3)//'out1')
  end if

  write(buf,'(A)') TRIM(inpname)//' >'//TRIM(outname2)//" 2>&1"
  write(iout,'(A)') '$$ORCA '//TRIM(buf)

  i = system(TRIM(hf_prog_path)//' '//TRIM(buf))
  if(i /= 0) then
   write(iout,'(/,A)') 'ERROR in subroutine do_scf_and_read_e: ORCA SCF job failed.'
   write(iout,'(A)') 'You can open file '//TRIM(outname2)//' and see why.'
   stop
  end if
  call read_hf_e_and_ss_from_orca_out(outname2, e, ssquare)

 case default
  write(iout,'(A)') 'ERROR in subroutine do_scf_and_read_e: invalid prog_name&
                   & = '//TRIM(prog_name)
  stop
 end select

 return
end subroutine do_scf_and_read_e

! read HF electronic energy from a Gaussian .log/.out file
subroutine read_hf_e_and_ss_from_gau_out(logname, e, ss)
 use print_id, only: iout
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e, ss ! HF energy and spin square
 character(len=240) :: buf
 character(len=240), intent(in) :: logname

 e = 0d0; ss = 0d0
 open(newunit=fid,file=TRIM(logname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == 'Entering Gaussian') then
   i = -1
   exit
  end if
  if(index(buf,'SCF Done') /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_gau_out: no&
                     & 'SCF Done' found in"
  write(iout,'(A)') 'file '//TRIM(logname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) e

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 i = index(buf, 'S**2')
 if(i /= 0) read(buf(i+6:),*) ss

 close(fid)
 return
end subroutine read_hf_e_and_ss_from_gau_out

! read HF electronic energy from a PySCF .out file
subroutine read_hf_e_and_ss_from_pyscf_out(outname, e, ss)
 use print_id, only: iout
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e, ss
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 e = 0d0; ss = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == 'Entering Gaussian') then
   i = -1
   exit
  end if
  if(index(buf,'SCF Done') /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_pyscf_out: no&
                     & 'SCF Done' found in"
  write(iout,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) e

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 i = index(buf, 'S**2')
 if(i /= 0) read(buf(i+6:),*) ss

 close(fid)
 return
end subroutine read_hf_e_and_ss_from_pyscf_out

! read HF electronic energy from a PSI4 .out file
subroutine read_hf_e_and_ss_from_psi4_out(outname, e, ss)
 use print_id, only: iout
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e, ss
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 e = 0d0; ss = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == 'Entering Gaussian') then
   i = -1
   exit
  end if
  if(index(buf,'SCF Done') /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_psi4_out: no&
                     & 'SCF Done' found in"
  write(iout,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) e

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 i = index(buf, 'S**2')
 if(i /= 0) read(buf(i+6:),*) ss

 close(fid)
 return
end subroutine read_hf_e_and_ss_from_psi4_out

! read HF electronic energy from a ORCA .out file
subroutine read_hf_e_and_ss_from_orca_out(outname, e, ss)
 use print_id, only: iout
 implicit none
 integer :: i, fid
 real(kind=8), intent(out) :: e, ss
 character(len=240) :: buf
 character(len=240), intent(in) :: outname

 e = 0d0; ss = 0d0
 open(newunit=fid,file=TRIM(outname),status='old',position='append')

 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(2:18) == 'Entering Gaussian') then
   i = -1
   exit
  end if
  if(index(buf,'SCF Done') /= 0) exit
 end do ! for while

 if(i /= 0) then
  write(iout,'(/,A)') "ERROR in subroutine read_hf_e_and_ss_from_orca_out: no&
                     & 'SCF Done' found in"
  write(iout,'(A)') 'file '//TRIM(outname)
  close(fid)
  stop
 end if

 i = index(buf, '=')
 read(buf(i+1:),*) e

 read(fid,'(A)') buf
 read(fid,'(A)') buf
 i = index(buf, 'S**2')
 if(i /= 0) read(buf(i+6:),*) ss

 close(fid)
 return
end subroutine read_hf_e_and_ss_from_orca_out

