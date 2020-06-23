! written by jxzou at 20200512: generate PySCF format basis set (.py file) from Gaussian .fch(k) file

! Note: this subroutine is actually a wrapper of two utilities 'fch2inp' and 'bas_gms2py',
!       thus 'fch2inp' and 'bas_gms2py' must be compiled as well.

program main
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=4) :: str
 character(len=240) :: fchname
 logical :: uhf

 i = iargc()
 if(i<1 .or. i>2) then
  write(iout,'(/,A)') ' ERROR in subroutine bas_fch2py: wrong command line arguments!'
  write(iout,'(/,A)') ' Example 1(RHF): bas_fch2py a.fch'
  write(iout,'(/,A,/)') ' Example 2(UHF): bas_fch2py a.fch -uhf'
  stop
 end if

 call getarg(1, fchname)
 uhf = .false.

 if(i == 2) then
  str = ' '
  call getarg(2, str)
  if(str /= '-uhf') then
   write(iout,'(A)') 'ERROR in subroutine bas_fch2py: wrong command line arguments!'
   write(iout,'(A)') "The 2nd argument(if provided) must be '-uhf'."
   stop
  end if
  uhf = .true.
 end if

 call bas_fch2py(fchname, uhf)
 stop
end program main

! generate PySCF format basis set (.py file) from Gaussian .fch(k) file
subroutine bas_fch2py(fchname, uhf)
 implicit none
 integer :: i, system
 character(len=240) :: inpname
 character(len=240), intent(in) :: fchname
 logical, intent(in) :: uhf
 logical :: cart

 call determine_sph_or_cart(fchname, cart) 

 i = index(fchname, '.fch', back=.true.)
 inpname = fchname(1:i-1)//'.inp'

 if(uhf) then
  i = system('fch2inp '//TRIM(fchname)//' -uhf >/dev/null')
 else ! RHF
  i = system('fch2inp '//TRIM(fchname)//' >/dev/null')
 end if

 if(cart) then ! Cartesian functions
  i = system('bas_gms2py '//TRIM(inpname))
 else          ! sperical harmonic functions
  i = system('bas_gms2py '//TRIM(inpname)//' -sph')
 end if

 ! delete the inpname
 open(newunit=i,file=TRIM(inpname),status='old')
 close(unit=i,status='delete')

 return
end subroutine bas_fch2py

