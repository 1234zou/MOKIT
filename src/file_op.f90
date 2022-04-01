! written by jxzou at 20210113: file operations

subroutine require_file_exist(fname)
 implicit none
 integer, parameter :: iout = 6
 character(len=240), intent(in) :: fname
 logical :: alive

 inquire(file=TRIM(fname),exist=alive)
 if(.not. alive) then
  write(iout,'(/,A)') 'ERROR in subroutine require_file_exist: file does not exist!'
  write(iout,'(A)') 'Filename='//TRIM(fname)
  stop
 end if
 return
end subroutine require_file_exist

subroutine open_file(fname, is_rewind, fid)
 implicit none
 integer :: i
 integer, parameter :: iout = 6
 integer, intent(out) :: fid
 character(len=240), intent(in) :: fname
 logical, intent(in) :: is_rewind

 if(is_rewind) then
  open(newunit=fid,file=TRIM(fname),status='old',position='rewind',iostat=i)
 else
  open(newunit=fid,file=TRIM(fname),status='old',position='append',iostat=i)
 end if

 if(i /= 0) then
  write(iout,'(/,A)') 'Failed to open file '//TRIM(fname)
  stop
 end if
 return
end subroutine open_file

! delete the specified file (if not exist, return)
subroutine delete_file(fname)
 implicit none
 integer :: fid
 character(len=*), intent(in) :: fname
 logical :: alive

 inquire(file=TRIM(fname),exist=alive)

 if(alive) then
  inquire(file=TRIM(fname),opened=alive,number=fid)
  if(.not. alive) open(newunit=fid,file=TRIM(fname),status='old')
  close(fid,status='delete')
 end if
end subroutine delete_file

! delete a set of files
subroutine delete_files(n, fname)
 implicit none
 integer :: i
 integer, intent(in) :: n
 character(len=240), intent(in) :: fname(n)

 do i = 1, n, 1
  call delete_file(fname(i))
 end do ! for i
end subroutine delete_files

subroutine delete_files_in_path(path, n, fname)
 implicit none
 integer :: i
 integer, intent(in) :: n
 character(len=240), intent(in) :: path, fname(n)

 do i = 1, n, 1
  call delete_file(TRIM(path)//'/'//TRIM(fname(i)))
 end do ! for i
end subroutine delete_files_in_path

! copy file fname1 to fname2 (if delete=.True., delete fname1)
subroutine copy_file(fname1, fname2, delete)
 implicit none
 integer :: i, fid1, fid2
 character(len=240) :: buf
 character(len=240), intent(in) :: fname1, fname2
 logical, intent(in) :: delete

 call require_file_exist(fname1)
 open(newunit=fid1,file=TRIM(fname1),status='old',position='rewind')
 open(newunit=fid2,file=TRIM(fname2),status='replace')

 do while(.true.)
  read(fid1,'(A)',iostat=i) buf
  if(i /= 0) exit
  write(fid2,'(A)') TRIM(buf)
 end do ! for while

 if(delete) then
  close(fid1, status='delete')
 else
  close(fid1)
 end if

 close(fid2)
 return
end subroutine copy_file

! copy binary file (if delete=.True., delete fname1)
subroutine copy_bin_file(fname1, fname2, delete)
 implicit none
 integer :: i, system
 integer, parameter :: iout = 6
 character(len=240), intent(in) :: fname1, fname2
 logical, intent(in) :: delete

#ifdef _WIN32
 i = system('copy /Y '//TRIM(fname1)//' '//TRIM(fname2)//' > NUL')
#else
 i = system('cp '//TRIM(fname1)//' '//TRIM(fname2))
#endif

 if(delete) then
#ifdef _WIN32
 i = system('del '//TRIM(fname1))
#else
 i = system('rm -f '//TRIM(fname1))
#endif
 end if

 if(i /= 0) then
  write(iout,'(A)') 'ERROR in subroutine copy_bin_file: fail to copy binary&
                   & file from'
  write(iout,'(A)') TRIM(fname1)//' to '//TRIM(fname2)
  stop
 end if
 return
end subroutine copy_bin_file

! create a directory
subroutine create_dir(dirname)
 implicit none
 integer :: i, system
 character(len=*), intent(in) :: dirname

#ifdef _WIN32
 i = system('CD '//dirname)
 if(i /= 0) i = system('MD '//dirname)
#else
 i = system('mkdir -p '//dirname)
#endif
 return
end subroutine create_dir

