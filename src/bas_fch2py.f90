! written by jxzou at 20200512: generate PySCF input file (.py) from Gaussian .fch(k) file
! updated by jxzou at 20200809: combined with util_wrapper.f90

! Note: this subroutine is actually a wrapper of two utilities 'fch2inp' and
! 'bas_gms2py', thus they must be compiled as well.

!TODO: to save the file conversion time, we may need to skip calling bas_gms2py
! and directly generate .py file from .fch(k) file.

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i, j
 character(len=5) :: buf
 character(len=240) :: fchname
 logical :: prt_dft, pbc, obj_only, rest

 i = iargc()
 if(i<1 .or. i>4) then
  write(6,'(/,A)') ' ERROR in program bas_fch2py: wrong command line arguments!'
  write(6,'(A)')   ' Example 1  (R(O)HF, UHF): bas_fch2py a.fch'
  write(6,'(A)')   ' Example 2          (DFT): bas_fch2py a.fch -dft'
  write(6,'(A)')   ' Example 3         (REST): bas_fch2py a.fch -dft -rest'
  write(6,'(A)')   ' Example 4       (PBC-HF): bas_fch2py a.fch -pbc'
  write(6,'(A)')   ' Example 5      (PBC-DFT): bas_fch2py a.fch -pbc -dft'
  write(6,'(A)')   ' Example 6  (object only): bas_fch2py a.fch -obj'
  write(6,'(A,/)') ' Example 7 (PBC obj only): bas_fch2py a.fch -pbc -obj'
  stop
 end if

 fchname = ' '
 call getarg(1, fchname)
 call require_file_exist(fchname)

 prt_dft=.false.; pbc=.false.; obj_only=.false.; rest=.false.

 if(i > 1) then
  do j = 2, i, 1
   call getarg(j, buf)
   buf = ADJUSTL(buf)
   select case(TRIM(buf))
   case('-dft')
    prt_dft = .true.
   case('-pbc')
    pbc = .true.
   case('-obj')
    obj_only = .true.
   case('-rest')
    rest = .true.
   case default
    write(6,'(/,A)') 'ERROR in subroutine bas_fch2py: wrong command line arguments!'
    write(6,'(A)') "The arguments can only be -dft/-pbc/-obj/-rest"
    stop
   end select
  end do ! for j
 end if

 if(obj_only .and. (prt_dft .or. rest)) then
  write(6,'(/,A)') 'ERROR in subroutine bas_fch2py: -obj is incompatible with -&
                   &dft/-rest.'
  stop
 end if

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call bas_fch2py(fchname, prt_dft, pbc, obj_only, rest)
end program main

! generate PySCF format basis set (.py file) from Gaussian .fch(k) file
subroutine bas_fch2py(fchname, prt_dft, pbc, obj_only, rest)
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, system, RENAME
 character(len=15) :: dftname
 character(len=240) :: inpname, inpname1, pyname
 character(len=240), intent(in) :: fchname
 character(len=300) :: command
 logical :: alive, cart, is_hf, rotype, untype
 logical, intent(in) :: prt_dft, pbc, obj_only, rest

 call find_specified_suffix(fchname, '.fch', i)

 ! if the user provides a .fchk file, copy this file to .fch
 if(INDEX(fchname, '.fchk') > 0) then
  inpname = fchname(1:i-1)//'.fch'
  call copy_file(fchname, inpname, .false.)
 end if

 inpname = fchname(1:i-1)//'.inp'
 inpname1 = fchname(1:i-1)//'.t'
 pyname = fchname(1:i-1)//'.py'

 ! if inpname already exists, rename it
 inquire(file=TRIM(inpname),exist=alive)
 if(alive) i = RENAME(TRIM(inpname), TRIM(inpname1))

 call determine_sph_or_cart(fchname, cart) 
 call fch2inp_wrap(fchname, .false., 0, 0, .true., .false.)

 command = 'bas_gms2py '//TRIM(inpname)
 if(.not. cart) command = TRIM(command)//' -sph'
 if(pbc) command = TRIM(command)//' -pbc'
 if(obj_only) command = TRIM(command)//' -obj'
 if(rest) command = TRIM(command)//' -rest'

 i = SYSTEM(TRIM(command))
 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine bas_fch2py: failed to call utility bas_&
                   &gms2py.'
  write(6,'(A)') 'The file '//TRIM(fchname)//' may be incomplete.'
  stop
 end if

 if(alive) then
  i = RENAME(TRIM(inpname1), TRIM(inpname))
 else
  ! delete the inpname
  open(newunit=i,file=TRIM(inpname),status='old')
  close(unit=i,status='delete')
 end if

 if(prt_dft) then
  call find_dftname_in_fch(fchname, dftname, is_hf, rotype, untype)
  call prt_dft_key2pyscf_script(dftname, is_hf, rotype, untype, pyname, pbc, rest)
 else
  if(rest) then
   !call find_dftname_in_fch(fchname, dftname, is_hf, rotype, untype)
   ! let HF scf run
   ! rotype and untype do not matter in this case
   dftname = ' '
   call prt_dft_key2pyscf_script(dftname, .true., .false., .false., pyname, &
                                 pbc, rest)
  end if
 end if
end subroutine bas_fch2py

! find the DFT name in a Gaussian .fch file
subroutine find_dftname_in_fch(fchname, dftname, is_hf, rotype, untype)
 implicit none
 integer :: i, j, fid
 character(len=240) :: buf
 character(len=480) :: longbuf
 character(len=240), intent(in) :: fchname
 character(len=15), intent(out) :: dftname
 logical, intent(out) :: is_hf, rotype, untype

 is_hf = .false.; rotype = .false.; untype = .false.; dftname = ' '
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == 'Route') then
   read(fid,'(A)') buf
   ! in case that the Route Section is long and be divided into 2 lines, let's
   ! read one more line
   read(fid,'(A)') longbuf
   if(longbuf(1:6) == 'Charge') then
    longbuf = buf
   else
    longbuf = TRIM(buf)//TRIM(longbuf)
   end if
   j = INDEX(longbuf,'/')
   i = INDEX(longbuf(1:j-1), ' ', back=.true.)
   dftname = longbuf(i+1:j-1)
   exit
  end if
 end do ! for while

 close(fid)
 if(LEN_TRIM(dftname) == 0) then
  write(6,'(A)') REPEAT('-',79)
  write(6,'(A)') 'Warning from subroutine find_dftname_in_fch: DFT name is not &
                 &found in file'
  write(6,'(A)') TRIM(fchname)
  write(6,'(A)') 'Possible reason: you are using a fch file generated by old ve&
                 &rsion of Gaussian,'
  write(6,'(A)') 'e.g. g03, g09. The density functional will be set to PBEPBE.'
  write(6,'(A)') REPEAT('-',79)
  dftname = 'pbepbe'
 end if

 call lower(dftname)
 if(dftname(1:2) == 'ro') then
  rotype = .true.
  dftname = dftname(3:)
 else if(dftname(1:1) == 'u') then
  untype = .true.
  dftname = dftname(2:)
 else if(dftname(1:1) == 'r') then
  dftname = dftname(2:)
 end if

 if(TRIM(dftname) == 'hf') is_hf = .true.

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 read(fid,'(A)') buf
 read(fid,'(A)') buf
 close(fid)

 if(buf(11:11) == 'U') then
  untype = .true.
 else if(buf(11:12) =='RO') then
  rotype = .true.
 end if
end subroutine find_dftname_in_fch

! delete the key 'scf' from buf
!subroutine del_scf_in_buf(buf, deleted)
! implicit none
! integer :: i, k
! character(len=240), intent(inout) :: buf
! logical, intent(out) :: deleted
!
! deleted = .false.
! i = INDEX(buf,'scf')
!
! if(buf(1:17)=='from pyscf import' .and. i>0) then
!  if(buf(i+3:i+3) /= ',') then
!   write(6,'(/,A)') "ERROR in subroutine del_scf_in_buf: no 'scf,' key found in&
!                    & buf. This"
!   write(6,'(A)') 'subroutine is not designed for this case.'
!   write(6,'(A)') 'buf='//TRIM(buf)
!   stop
!  end if
!  ! now we have 'scf,'
!  k = LEN_TRIM(buf)
!  buf = TRIM(buf(1:i-1))//' '//ADJUSTL(buf(i+4:k))
!  deleted = .true.
! end if
!end subroutine del_scf_in_buf

! print DFT keywords into a PySCF .py script
subroutine prt_dft_key2pyscf_script(dftname, is_hf, rotype, untype, pyname, pbc,&
                                    rest)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=4) :: str4
 character(len=15) :: dftname1
 character(len=240) :: buf, pyname1, pchkname
 character(len=15), intent(in) :: dftname
 character(len=240), intent(in) :: pyname
 logical, intent(in) :: is_hf, rotype, untype, pbc, rest
 logical :: printed

 str4 = 'mol'
 if(pbc) str4 = 'cell'

 call find_specified_suffix(pyname, '.py', i)
 pyname1 = pyname(1:i-1)//'.t'
 pchkname = pyname(1:i-1)//'.pchk'
 open(newunit=fid,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(pyname1),status='replace')

 if(.not. is_hf) then ! KS-DFT
  if(pbc) then
   write(fid1,'(A)') 'from pyscf.pbc import dft'
  else
   printed = .false.
   do i = 1, 3
    read(fid,'(A)') buf
    if((.not. printed) .and. buf(1:17)=='from pyscf import') then
     buf = TRIM(buf)//', dft'
     printed = .true.
    end if
    write(fid1,'(A)') TRIM(buf)
   end do ! for i
  end if
 end if

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:3) == '#dm') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while
 close(fid,status='delete')

 dftname1 = dftname
 select case(TRIM(dftname1))
 case('b3lyp')
  dftname1 = 'b3lypg' ! Since PySCF-2.3, B3LYP in PySCF is the same to that in
  ! Gaussian. In older versions of PySCF, the functional is b3lypg. So here
  ! b3lypg is used for compatibility.
 case('pbepbe')
  dftname1 = 'pbe,pbe'
 case('pbe1pbe')
  dftname1 = 'pbe0' ! PySCF also recognizes pbe1pbe
 case('hseh1pbe')
  dftname1 = 'hse06'
 end select
 
 write(fid1,'(A)') 'dm = mf.make_rdm1()'

 if(.not. is_hf) then ! KS-DFT
  if(untype) then
   write(fid1,'(A)',advance='no') 'mf = dft.UKS('
  else if(rotype) then
   write(fid1,'(A)',advance='no') 'mf = dft.ROKS('
  else
   write(fid1,'(A)',advance='no') 'mf = dft.RKS('
  end if
  write(fid1,'(A)') TRIM(str4)//')'
  write(fid1,'(A)') "mf.xc = '"//TRIM(dftname1)//"'"
  write(fid1,'(A)') 'mf.grids.atom_grid = (99,590)' ! ultrafine
 end if

 write(fid1,'(A)') 'mf.verbose = 4'
 write(fid1,'(A)') 'mf.max_cycle = 128'
 if(rest) write(fid1,'(A)') "mf.chkfile = '"//TRIM(pchkname)//"'"
 write(fid1,'(A)') 'mf.kernel(dm0=dm)'

 close(fid1)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
end subroutine prt_dft_key2pyscf_script

