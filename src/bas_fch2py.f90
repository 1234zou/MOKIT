! written by jxzou at 20200512: generate PySCF format basis set (.py file) from Gaussian .fch(k) file
! updated by jxzou at 20200809: combined with util_wrapper.f90

! Note: this subroutine is actually a wrapper of two utilities 'fch2inp' and 'bas_gms2py',
!       thus 'fch2inp' and 'bas_gms2py' must be compiled as well.

program main
 use util_wrapper, only: formchk
 implicit none
 integer :: i
 character(len=240) :: fchname

 i = iargc()
 if(i /= 1) then
  write(6,'(/,A)') ' ERROR in subroutine bas_fch2py: wrong command line argument!'
  write(6,'(A,/)') ' Example (R(O)HF, UHF): bas_fch2py a.fch'
  stop
 end if

 call getarg(1, fchname)
 call require_file_exist(fchname)

 ! if .chk file provided, convert into .fch file automatically
 i = LEN_TRIM(fchname)
 if(fchname(i-3:i) == '.chk') then
  call formchk(fchname)
  fchname = fchname(1:i-3)//'fch'
 end if

 call bas_fch2py(fchname)
end program main

! generate PySCF format basis set (.py file) from Gaussian .fch(k) file
subroutine bas_fch2py(fchname)
 use util_wrapper, only: fch2inp_wrap
 implicit none
 integer :: i, system, RENAME
 character(len=240) :: inpname, inpname1
 character(len=240), intent(in) :: fchname
 logical :: alive, cart

 i = index(fchname, '.fch', back=.true.)
 if(i == 0) then
  write(6,'(A)') "ERROR in subroutine bas_fch2py: '.fch' not found in filename&
                & "//TRIM(fchname)
  stop
 end if
 inpname = fchname(1:i-1)//'.inp'
 inpname1 = fchname(1:i-1)//'.inp.t'

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
end subroutine bas_fch2py

