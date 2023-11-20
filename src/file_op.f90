! written by jxzou at 20210113: file operations

subroutine require_file_exist(fname)
 implicit none
 character(len=240), intent(in) :: fname
 logical :: alive

 inquire(file=TRIM(fname),exist=alive)
 if(.not. alive) then
  write(6,'(/,A)') 'ERROR in subroutine require_file_exist: file does not exist!'
  write(6,'(A)') 'Filename='//TRIM(fname)
  stop
 end if
end subroutine require_file_exist

subroutine open_file(fname, is_rewind, fid)
 implicit none
 integer :: i
 integer, intent(out) :: fid
 character(len=240), intent(in) :: fname
 logical, intent(in) :: is_rewind

 call require_file_exist(fname)
 if(is_rewind) then
  open(newunit=fid,file=TRIM(fname),status='old',position='rewind',iostat=i)
 else
  open(newunit=fid,file=TRIM(fname),status='old',position='append',iostat=i)
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine open_file: failed to open file '//&
                    TRIM(fname)
  stop
 end if
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

subroutine delete_cfour_tmp_files(nproc)
 implicit none
 integer :: i
 integer, intent(in) :: nproc
 integer, parameter :: nfile = 25
 character(len=7) :: dirname
 character(len=12) :: filelist(nfile)
 data filelist /'BASINFO.DATA','DIPOL','FILES','fort.81','fort.82','fort.83',&
  'GAMLAM','HFNES1','HFNES2','HFNES3','HFNES4','HFNES5','HFNES6','HFNES7',&
  'IIII','JAINDX','JMOLplot','JOBARC','MOABCD','MOINTS','MOL','MOLDEN','ncpu',&
  'NTOTAL','OPTARC'/

 if(nproc > 0) then
  do i = 0, nproc-1, 1
   write(dirname,'(A4,I3.3)') 'rank', i
   call remove_dir(TRIM(dirname))
  end do ! for i
 end if

 do i = 1, nfile, 1
  call delete_file(TRIM(filelist(i)))
 end do ! for i
end subroutine delete_cfour_tmp_files

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

 close(fid2)
 if(delete) then
  close(fid1, status='delete')
 else
  close(fid1)
 end if
end subroutine copy_file

! Copy a file using system commands. If delete=.True., delete fname1.
subroutine sys_copy_file(fname1, fname2, delete)
 implicit none
 integer :: i, SYSTEM
 character(len=*), intent(in) :: fname1, fname2
 logical, intent(in) :: delete

#ifdef _WIN32
 i = SYSTEM('copy /Y '//TRIM(fname1)//' '//TRIM(fname2)//' > NUL')
#else
 i = SYSTEM('cp '//TRIM(fname1)//' '//TRIM(fname2))
#endif

 if(delete) then
#ifdef _WIN32
 i = SYSTEM('del '//TRIM(fname1))
#else
 i = SYSTEM('rm -f '//TRIM(fname1))
#endif
 end if

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine sys_copy_file: fail to copy file from'
  write(6,'(A)') TRIM(fname1)//' to '//TRIM(fname2)
  stop
 end if
end subroutine sys_copy_file

! copy a file from a path/directory to the current directory
subroutine copy_file_from_a_path(path, fname)
 implicit none
 character(len=1) :: s
 character(len=700) :: buf
 character(len=*), intent(in) :: path, fname

#ifdef _WIN32
 s = '\'
#else
 s = '/'
#endif

 buf = TRIM(path)//s//TRIM(fname)
 call sys_copy_file(TRIM(buf), TRIM(fname), .false.)
end subroutine copy_file_from_a_path

! create a directory
subroutine create_dir(dirname)
 implicit none
 integer :: i, SYSTEM
 character(len=*), intent(in) :: dirname

#ifdef _WIN32
 i = SYSTEM('CD '//dirname)
 if(i == 0) then
  i = SYSTEM('CD ..')
 else
  i = SYSTEM('MD '//dirname)
 end if
#else
 i = SYSTEM('mkdir -p '//dirname)
#endif
end subroutine create_dir

! remove/delete a directory. USE WITH CAUTION!!!
subroutine remove_dir(dirname)
 implicit none
 integer :: i, SYSTEM
 character(len=*), intent(in) :: dirname

#ifdef _WIN32
 i = SYSTEM('rd '//dirname)
#else
 ! in case that important directories are deleted, check the dirname firstly
 i = LEN_TRIM(dirname)
 if(i == 1) then
  write(6,'(A)') 'ERROR in subroutine remove_dir: wrong directory name: '//&
                  TRIM(dirname)
  stop
 end if

 if(i > 3) then
  select case(dirname(1:4))
  case('/usr','/bin','/lib','/etc')
   write(6,'(A)') 'ERROR in subroutine remove_dir: dangerous directory name: '&
                   //TRIM(dirname)
   stop
  case default
  end select
 end if

 if(i > 4) then
  if(dirname(1:5) == '/sbin') then
   write(6,'(A)') 'ERROR in subroutine remove_dir: dangerous directory name: '&
                   //TRIM(dirname)
   stop
  end if
 end if

 if(i > 5) then
  if(dirname(1:6) == '/lib64') then
   write(6,'(A)') 'ERROR in subroutine remove_dir: dangerous directory name: '&
                   //TRIM(dirname)
   stop
  end if
 end if

 i = SYSTEM('rm -rf '//dirname)
#endif
end subroutine remove_dir

! export a real(kind=8) array into a plain .txt file
subroutine export_rarray2txt(txtname, label, n, a)
 implicit none
 integer :: i, fid
 integer, intent(in) :: n
!f2py intent(in) :: n
 real(kind=8), intent(in) :: a(n)
!f2py intent(in) :: a
!f2py depend(n) :: a
 character(len=25), parameter :: star = '*************************'
 character(len=*), intent(in) :: label
!f2py intent(in) :: label
 character(len=240), intent(in) :: txtname
!f2py intent(in) :: txtname

 open(newunit=fid,file=TRIM(txtname),status='replace')
 write(fid,'(A,I0)') 'n=', n
 write(fid,'(A)') ' '//star//' '//label//' '//star

 do i = 1, n, 1
  write(fid,'(I5,1X,F11.5)') i, a(i)
 end do ! for i

 close(fid)
end subroutine export_rarray2txt

! export a square matrix into a plain .txt file
subroutine export_mat_into_txt(txtname, n, mat, lower, label)
 implicit none
 integer :: i, j, k, m, na, fid
 integer, intent(in) :: n
 integer, allocatable :: a(:)
 real(kind=8), intent(in) :: mat(n,n)
 character(len=25), parameter :: star = '*************************'
 character(len=25), intent(in) :: label
 character(len=240), intent(in) :: txtname
 logical, intent(in) :: lower ! True: lower triangle

 open(newunit=fid,file=TRIM(txtname),status='replace')
 write(fid,'(A,L1)') 'Lower_Triangle=',lower
 write(fid,'(A,I0)') 'n=', n
 write(fid,'(A)') ' '//star//' '//label//' '//star

 m = n/5
 if(n-5*m > 0) m = m + 1

 if(lower) then ! lower triangle
  do i = 1, m, 1
   if(i < m) then
    na = 5
   else
    na = n - 5*(m-1)
   end if
   allocate(a(na))
   forall(j=1:na) a(j) = 5*i - 5 + j
   write(fid,'(5I14)') a(1:na)
   deallocate(a)

   k = 5*i - 4
   do j = k, n, 1
    write(fid,'(I6,5(1X,ES15.8))') j, mat(k:min(k+4,j),j)
   end do ! for j
  end do ! for i

 else           ! full matrix
  do i = 1, m, 1
   if(i < m) then
    na = 5
   else
    na = n - 5*(m-1)
   end if
   allocate(a(na))
   forall(j=1:na) a(j) = 5*i - 5 + j
   write(fid,'(5I14)') a(1:na)
   deallocate(a)

   k = 5*i - 4
   do j = 1, n, 1
    write(fid,'(I6,5(1X,ES15.8))') j, mat(j,k:k+na-1)
   end do ! for j
  end do ! for i
 end if

 close(fid)
end subroutine export_mat_into_txt

! simplify a .fch(k) file
subroutine simplify_fch(fchname)
 implicit none
 integer :: i, fid, fid1, RENAME
 character(len=240) :: buf, fchname1
 character(len=240), intent(in) :: fchname

 i = INDEX(fchname, '.fch', back=.true.)
 fchname1 = fchname(1:i-1)//'.t'

 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')
 open(newunit=fid1,file=TRIM(fchname1),status='replace')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:12)=='Mulliken Cha' .or. buf(1:12)=='Number of ex') exit
  write(fid1,'(A)') TRIM(buf)
 end do ! for while

 if(i /= 0) then
  write(6,'(A)') 'ERROR in subroutine simplify_fch: keywords for termination no&
                 &t found in file '//TRIM(fchname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine simplify_fch

