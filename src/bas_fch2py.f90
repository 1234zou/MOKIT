! written by jxzou at 20200512: generate PySCF format basis set (.py file) from Gaussian .fch(k) file
! updated by jxzou at 20200809: combined with util_wrapper.f90

! Note: this subroutine is actually a wrapper of two utilities 'fch2inp' and 'bas_gms2py',
!       thus 'fch2inp' and 'bas_gms2py' must be compiled as well.

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i
 character(len=4) :: str
 character(len=240) :: fchname
 logical :: prt_dft

 i = iargc()
 if(i<1 .or. i>2) then
  write(6,'(/,A)') ' ERROR in program bas_fch2py: wrong command line argument!'
  write(6,'(A)')   ' Example 1 (R(O)HF, UHF): bas_fch2py a.fch'
  write(6,'(A,/)') ' Example 2         (DFT): bas_fch2py a.fch -dft'
  stop
 end if

 call getarg(1, fchname)
 call require_file_exist(fchname)

 prt_dft = .false.
 if(i == 2) then
  call getarg(2, str)
  select case(str)
  case('-dft')
   prt_dft = .true.
  case default
   write(6,'(/,A)') "ERROR in program bas_fch2py: the 2nd argument can only be&
                    & '-dft'."
   write(6,'(A)') "But your input is '"//str//"'."
   stop
  end select
 end if

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call bas_fch2py(fchname, prt_dft)
end program main

! generate PySCF format basis set (.py file) from Gaussian .fch(k) file
subroutine bas_fch2py(fchname, prt_dft)
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, system, RENAME
 character(len=15) :: dftname
 character(len=240) :: inpname, inpname1, pyname
 character(len=240), intent(in) :: fchname
 logical :: alive, cart
 logical, intent(in) :: prt_dft

 i = index(fchname, '.fch', back=.true.)
 if(i == 0) then
  write(6,'(A)') "ERROR in subroutine bas_fch2py: '.fch' suffix not found in fi&
                 &lename "//TRIM(fchname)
  stop
 end if
 inpname = fchname(1:i-1)//'.inp'
 inpname1 = fchname(1:i-1)//'.t'
 pyname = fchname(1:i-1)//'.py'

 ! if inpname already exists, rename it
 inquire(file=TRIM(inpname),exist=alive)
 if(alive) i = RENAME(TRIM(inpname), TRIM(inpname1))

 call determine_sph_or_cart(fchname, cart) 
 call fch2inp_wrap(fchname, .false., 0, 0)

 if(cart) then ! Cartesian functions
  i = system('bas_gms2py '//TRIM(inpname))
 else          ! sperical harmonic functions
  i = system('bas_gms2py '//TRIM(inpname)//' -sph')
 end if

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine bas_fch2py: call utility bas_gms2py failed.'
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
  call find_dftname_in_fch(fchname, dftname)
  call prt_dft_key2pyscf_script(dftname, pyname)
 end if
end subroutine bas_fch2py

! find the DFT name in a Gaussian .fch file
subroutine find_dftname_in_fch(fchname, dftname)
 implicit none
 integer :: i, j, fid
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
 character(len=15), intent(out) :: dftname

 dftname = ' '
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:5) == 'Route') then
   read(fid,'(A)') buf
   j = index(buf,'/')
   i = index(buf(1:j-1), ' ')
   dftname = buf(i+1:j-1)
   exit
  end if
 end do ! for while

 close(fid)
 if(LEN_TRIM(dftname) == 0) then
  write(6,'(/,A)') 'Warning from subroutine find_dftname_in_fch: DFT name is no&
                   &t found in file '//TRIM(fchname)
  stop
 end if
end subroutine find_dftname_in_fch

! print DFT keywords into a PySCF .py script
subroutine prt_dft_key2pyscf_script(dftname, pyname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=15) :: dftname1
 character(len=240) :: buf, pyname1
 character(len=15), intent(in) :: dftname
 character(len=240), intent(in) :: pyname
 logical :: prt

 pyname1 = TRIM(pyname)//'.t'
 open(newunit=fid,file=TRIM(pyname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(pyname1),status='replace')
 prt = .false.

 do i = 1, 3
  read(fid,'(A)') buf
  if((.not.prt) .and. buf(1:17)=='from pyscf import') then
   buf = TRIM(buf)//', dft'
   prt = .true.
  end if
  write(fid1,'(A)') TRIM(buf)
 end do ! for i

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 close(fid,status='delete')
 dftname1 = dftname
 call lower(dftname1)
 if(TRIM(dftname1) == 'b3lyp') dftname1 = 'b3lypg'
 write(fid1,'(A)') 'dm = mf.make_rdm1()'
 write(fid1,'(A)') 'mf = dft.RKS(mol)'
 write(fid1,'(A)') "mf.xc = '"//TRIM(dftname1)//"'"
 write(fid1,'(A)') 'mf.grids.atom_grid = (99,590)'
 write(fid1,'(A)') 'mf.kernel(dm0=dm)'
 close(fid1)
 i = RENAME(TRIM(pyname1), TRIM(pyname))
end subroutine prt_dft_key2pyscf_script

