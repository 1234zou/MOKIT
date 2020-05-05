! written by jxzou at 20200420: automatically do multireference calculations

! If the input file is .gjf,
!  1) if singlet, RHF, UHF(guess=mix stable=opt) will be successively performed;
!  2) if not singlet, UHF(stable=opt) will be performed.
! Otherwise the input file must be .fch(k) (which is recommended for difficult SCF case),
!  the .fch(k) file must contain a UHF job.

program automr
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 character(len=240) :: fname = ' '

 i = iargc()
 if(i /= 1) then
  write(iout,'(/,A)') ' ERROR in program automr: wrong command line arguments!'
  write(iout,'(A)')   ' Example 1: automr a.gjf'
  write(iout,'(A,/)') ' Example 2: automr a.fchk'
  stop
 end if

 call getarg(1, fname)

 stop
end program automr

