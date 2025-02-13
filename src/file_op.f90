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
 if(LEN_TRIM(fname2) == 0) then
  write(6,'(/,A)') 'ERROR in subroutine copy_file: empty filename for fname2.'
  write(6,'(A)') 'For bug tracking: fname1='//TRIM(fname1)
  stop
 end if

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

 if(LEN_TRIM(fname1) == 0) then
  write(6,'(/,A)') 'ERROR in subroutine sys_copy_file: empty filename for fname1.'
  write(6,'(A)') 'For bug tracking: fname2='//TRIM(fname2)
  stop
 end if

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

! move a file into a specified directory
subroutine move_file(fname, dirname)
 implicit none
 integer :: i, SYSTEM
 character(len=*), intent(in) :: fname, dirname

#ifdef _WIN32
 i = SYSTEM('move '//TRIM(fname)//' '//TRIM(dirname)//'\')
#else
 i = SYSTEM('mv '//TRIM(fname)//' '//TRIM(dirname)//'/')
#endif

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine move_file: failed to move file to the s&
                   &pecified directory.'
  write(6,'(A)') 'fname='//TRIM(fname)
  write(6,'(A)') 'dirname='//TRIM(dirname)
  stop
 end if
end subroutine move_file

! move several files into a specified directory
subroutine move_files(n, fname, dirname)
 implicit none
 integer :: i
 integer, intent(in) :: n
 character(len=240), intent(in) :: fname(n), dirname

 do i = 1, n, 1
  call move_file(TRIM(fname(i)), TRIM(dirname))
 end do ! for i
end subroutine move_files

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

! check if MOKIT_ROOT present. Set it if not present.
subroutine check_mokit_root()
 implicit none
 integer :: i, SYSTEM
 character(len=240) :: mokit_root

 write(6,'(/,A)') 'Checking MOKIT_ROOT ... '
 mokit_root = ' '
 call getenv('MOKIT_ROOT', mokit_root)

 if (LEN_TRIM(mokit_root) == 0) then
  ! assume we are under conda install
  i = SYSTEM("echo `get_mokit_loc.py` > ~/.mokitrc ")
  !i = SYSTEM("echo `get_mokit_loc.py` ")
  !i = SYSTEM('echo $MOKIT_ROOT')
  if(i /= 0) then
   write(6,'(/,A)') 'ERROR in subroutine check_mokit_root: '
   write(6,'(A)') '    fail to set MOKIT_ROOT for conda-installed MOKIT.'
   write(6,'(A)') 'If MOKIT is installed via conda, please report this issue.'
   write(6,'(A)') 'Otherwise, it means your MOKIT_ROOT is not properly set.'
   stop
  else
   write(6,'(A)') 'reset MOKIT_ROOT for conda-installed version'
  end if
 endif
end subroutine check_mokit_root

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
!f2py intent(in) :: fchname

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
  write(6,'(/,A)') 'ERROR in subroutine simplify_fch: keywords for termination &
                   &not found in'
  write(6,'(A)') 'file '//TRIM(fchname)
  close(fid)
  close(fid1,status='delete')
  stop
 end if

 close(fid,status='delete')
 close(fid1)
 i = RENAME(TRIM(fchname1), TRIM(fchname))
end subroutine simplify_fch

! read the number of MOs from a Gaussian .fch(k) file
subroutine read_nif_from_fch(fchname, nif)
 implicit none
 integer :: i, fid
 integer, intent(out) :: nif
!f2py intent(out) :: nif
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 nif = 0
 open(newunit=fid,file=TRIM(fchname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:13) == 'Number of ind') exit
 end do ! for while

 close(fid)
 if(i /= 0) then
  write(6,'(/,A)') "ERROR in subroutine read_nif_from_fch: no 'Number of ind' f&
                   &ound in file "//TRIM(fchname)
  stop
 end if

 read(buf(45:),*) nif
end subroutine read_nif_from_fch

! read spin multipliticity from a given .fch(k) file
subroutine read_mult_from_fch(fchname, mult)
 implicit none
 integer :: i, fid
 integer, intent(out) :: mult
!f2py intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: fchname
!f2py intent(in) :: fchname

 call open_file(fchname, .true., fid)
 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:4) == 'Mult') exit
 end do ! for while
 close(fid)

 if(i /= 0) then
  write(6,'(A)') "ERROR in subroutine read_mult_from_fch: no 'Mult' found in&
                & file "//TRIM(fchname)
  stop
 end if
 read(buf(50:),*) mult
end subroutine read_mult_from_fch

! read spin multiplicity from a specified ORCA input file
subroutine read_mult_from_orca_inp(inpname, mult)
 implicit none
 integer :: i, j, itype, fid
 integer, intent(out) :: mult
!f2py intent(out) :: mult
 character(len=240) :: buf
 character(len=240), intent(in) :: inpname
!f2py intent(in) :: inpname

 mult = 1; itype = 1
 open(newunit=fid,file=TRIM(inpname),status='old',position='rewind')

 do while(.true.)
  read(fid,'(A)',iostat=i) buf
  if(i /= 0) exit
  if(buf(1:1) == '*') exit       ! * xyz 0 1
  if(buf(1:7) == '%coords') then ! Mult = 1
   itype = 2
   exit
  end if
 end do ! for while

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_mult_from_orca_inp: spin multiplic&
                   &ity cannot be'
  write(6,'(A)') 'found in file '//TRIM(inpname)
  close(fid)
  stop
 end if

 select case(itype)
 case(1)
  close(fid)
  i = LEN_TRIM(buf)
  j = INDEX(buf(1:i), ' ', back=.true.)
  read(buf(j+1:i),*) mult
 case(2)
  do i = 1, 4
   read(fid,'(A)') buf
   if(INDEX(buf(1:5),'Mult') > 0) exit
  end do ! for i
  close(fid)
  if(i == 5) then
   write(6,'(/,A,I0)') 'ERROR in subroutine read_mult_from_orca_inp: spin multi&
                       &plicity cannot be'
   write(6,'(A)') 'found in file '//TRIM(inpname)
   stop
  end if
  i = INDEX(buf, '=')
  read(buf(i+1:),*) mult
 case default
  write(6,'(/,A,I0)') 'ERROR in subroutine read_mult_from_orca_inp: invalid ity&
                      &pe=', itype
  write(6,'(A)') 'inpname='//TRIM(inpname)
  close(fid)
  stop
 end select
end subroutine read_mult_from_orca_inp

subroutine write_grad_into_fch(fchname, natom, grad)
 implicit none
 integer :: fid
 integer, intent(in) :: natom
 real(kind=8), intent(in) :: grad(3*natom)
 character(len=240), intent(in) :: fchname

 open(newunit=fid,file=TRIM(fchname),status='old',position='append')
 write(fid,'(A,4X,I8)') 'Cartesian Gradient                         R   N=',natom*3
 write(fid,'(5(1X,ES15.8))') grad
 close(fid)
end subroutine write_grad_into_fch

! read the size of a symmetric matrix from .npy file
subroutine read_sym_mat_size_from_npy(npyname, n)
 implicit none
 integer :: i, fid
 integer, intent(out) :: n
 character(len=128) :: str
 character(len=240), intent(in) :: npyname

 n = 0; str = ' '
 open(newunit=fid,file=TRIM(npyname),status='old',form='unformatted',&
      access='stream')
 read(fid) str
 close(fid)

 i = INDEX(str, '(')
 read(str(i+1:),*) n
end subroutine read_sym_mat_size_from_npy

! read a symmetric matrix from .npy file
subroutine read_sym_mat_from_npy(npyname, n, a)
 implicit none
 integer :: i, fid
 integer, intent(in) :: n
 real(kind=8), intent(out) :: a(n,n)
 character(len=128) :: str
 character(len=240), intent(in) :: npyname

 str = ' '
 open(newunit=fid,file=TRIM(npyname),status='old',form='unformatted',&
      access='stream')
 read(fid) str
 read(fid,iostat=i) a
 close(fid)

 if(i /= 0) then
  write(6,'(/,A)') 'ERROR in subroutine read_sym_mat_from_npy: failed to read s&
                   &ymmetric matrix from'
  write(6,'(A)') 'the file '//TRIM(npyname)
  stop
 end if
end subroutine read_sym_mat_from_npy

! find and delete the target .pyc file
subroutine find_and_del_pyc(proname, py_ver)
 implicit none
 integer :: i
 character(len=5) :: ver
 character(len=300) :: pycname, py_ver1
 character(len=240), intent(in) :: proname, py_ver
!f2py intent(in) :: proname, py_ver

 if(py_ver(1:6) == 'Python') then
  py_ver1 = ADJUSTL(TRIM(py_ver(7:)))
 else
  py_ver1 = TRIM(py_ver)
 end if

 i = INDEX(py_ver1, '.')
 ver = py_ver1(1:i-1)
 py_ver1 = ADJUSTL(py_ver1(i+1:))
 i = INDEX(py_ver1, '.')
 ver = TRIM(ver)//py_ver1(1:i-1)

 pycname = '__pycache__/'//TRIM(proname)//'.cpython-'//TRIM(ver)//'.pyc'
 call delete_file(TRIM(pycname))
end subroutine find_and_del_pyc

! modify the file xxx_uno.txt
subroutine modify_uno_out(uno_out, ndb, npair, nopen)
 implicit none
 integer :: k, fid, idx(3)
 integer, intent(in) :: ndb, npair, nopen
!f2py intent(in) :: ndb, npair, nopen
 character(len=240) :: buf
 character(len=240), intent(in) :: uno_out
!f2py intent(in) :: uno_out

 open(newunit=fid,file=TRIM(uno_out),status='old',position='append')
 do while(.true.)
  BACKSPACE(fid)
  BACKSPACE(fid)
  read(fid,'(A)') buf
  if(buf(1:3) == 'ndb') exit
 end do ! for while

 BACKSPACE(fid)
 k = ndb + 1
 idx = [k, k+nopen+2*npair, nopen]
 write(fid,'(A6,I5)') 'ndb  =', ndb
 write(fid,'(A6,I5)') 'nact =', npair+nopen
 write(fid,'(A6,I5)') 'nact0=', npair
 write(fid,'(A6,3I5)')'idx  =', idx
 close(fid)
end subroutine modify_uno_out

! print the Cartesian coordinates of LMO/MLWF centers
subroutine prt_mo_center2xyz(nmo, coor, xyzname)
 implicit none
 integer :: i, fid
 integer, intent(in) :: nmo
!f2py intent(in) :: nmo
 real(kind=8), intent(in) :: coor(3,nmo)
!f2py intent(in) :: coor
!f2py depend(nmo) :: coor
 character(len=240), intent(in) :: xyzname
!f2py intent(in) :: xyzname

 open(newunit=fid,file=TRIM(xyzname),status='replace')
 write(fid,'(I0,/)') nmo

 do i = 1, nmo, 1
  write(fid,'(A,3(1X,F18.8))') 'X', coor(:,i)
 end do ! for i
 close(fid)
end subroutine prt_mo_center2xyz

! Note: the parameter mem must be provided in unit MB.
subroutine gen_gau_opt_gjf(gjfname, mem, charge, mult, natom, elem, coor, &
                           numfreq)
 implicit none
 integer :: i, fid
 integer,intent(in) :: mem, charge, mult, natom
 real(kind=8), intent(in) :: coor(3,natom)
 character(len=2), intent(in) :: elem(natom)
 character(len=240), intent(in) :: gjfname
 logical, intent(in) :: numfreq

 open(newunit=fid,file=TRIM(gjfname),status='replace')
 write(fid,'(A,I0,A)') '%mem=', min(mem,6000), 'MB'
 write(fid,'(A)') '%nprocshared=1'
 write(fid,'(A)',advance='no') '# opt(nomicro,maxcycles=300)'
 if(numfreq) write(fid,'(A)',advance='no') ' freq=numer'
 write(fid,'(A)') " def2SVPP nosymm external='gau_external'"
 ! TODO: automatically switch to UGBS or mixed basis set when there is any
 ! element which is out of range of def2SVPP.

 write(fid,'(/,A)') 'Using Gaussian as the geometry optimizer'
 write(fid,'(/,I0,1X,I0)') charge, mult
 do i = 1, natom, 1
  write(fid,'(A2,3(1X,F18.8))') elem(i), coor(:,i)
 end do ! for i

 write(fid,'(/)',advance='no')
 close(fid)
end subroutine gen_gau_opt_gjf

