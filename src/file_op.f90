! written by jxzou at 20210113: file operations

subroutine require_file_exist(fname)
 implicit none
 integer, parameter :: iout = 6
 character(len=240), intent(in) :: fname
 logical :: alive

 inquire(file=TRIM(fname),exist=alive)
 if(.not. alive) then
  write(iout,'(A)') 'ERROR in subroutine require_file_exist: file does not exist!'
  write(iout,'(A)') 'Filename='//TRIM(fname)
  stop
 end if
 return
end subroutine require_file_exist

