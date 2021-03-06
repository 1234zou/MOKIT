! written by jxzou at 20220119: convert one or more .chk file(s) to .gbw file(s)

program chk2gbw_main
 use util_wrapper, only: iout, chk2gbw
 implicit none
 integer :: i, n
 character(len=240) :: chkname = ' '

 n = iargc()
 if(n == 0) then
  write(iout,'(/,A)')' ERROR in program chk2gbw_main: wrong command line arguments!'
  write(iout,'(A,/)')' Format: chk2gbw 1.chk [2.chk ...]'
  stop
 end if

 do i = 1, n, 1
  call getarg(i, chkname)
  call chk2gbw(chkname)
 end do ! for i
 stop
end program chk2gbw_main

