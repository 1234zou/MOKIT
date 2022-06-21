! writteen by jxzou at 20220621

program main
 use util_wrapper, only: fch2gbw
 implicit none
 integer :: i
 character(len=240) :: fchname, gbwname

 i = iargc()
 if(i<1 .or. i>2) then
  write(6,'(/,A)') ' ERROR in subroutine fch2gbw: wrong command line arguments!'
  write(6,'(A)')   ' Example 1: fch2gbw a.fch'
  write(6,'(A,/)') ' Example 2: fch2gbw a.fch b.gbw'
  stop
 end if

 call getarg(1, fchname)
 call require_file_exist(fchname)

 gbwname = ' '
 if(i == 1) then
  call fch2gbw(fchname)
 else
  call getarg(2, gbwname)
  call fch2gbw(fchname, gbwname)
 end if
end program main

